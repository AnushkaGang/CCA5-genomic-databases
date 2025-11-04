"""
fasta_processor.py — CCA5 · Part A · Q1
Robust FASTA handlers with streaming, metadata extraction, indexing, and writing.

Features
--------
1) Parse FASTA (multiple sequences) as a memory-efficient generator
2) Extract headers & basic metadata
3) Handle large genomic files efficiently (streaming, no full read)
4) Write FASTA with proper 60-char wrapping
5) Optional on-disk random-access index (.fai-like) for fast sub-sequence fetches
6) CLI with common tasks (inspect, stats, rewrap/write, index, fetch range)

No external libraries. Python 3.10+.

How to run
----------
$ python fasta_processor.py

Menu options will guide you.

Complexity notes
----------------
- Streaming parse: O(total_bytes) time, O(1) extra RAM (besides current sequence buffer).
- Writing: O(n) time, O(1) extra space.
- Index build: O(file_size) time (one pass), index ~ O(#records).
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Generator, Iterable, Tuple, Optional, Dict


FASTA_WRAP = 60  # standard wrapping width when writing


@dataclass(slots=True)
class FastaRecord:
    """Simple container for one FASTA record."""
    header: str      # full header line without '>'
    seq: str         # sequence (A/T/G/C/N etc.), uppercased


# ---------------------- FASTA PARSER (streaming) ---------------------- #

def parse_fasta(path: str) -> Generator[FastaRecord, None, None]:
    """
    Stream-parse a FASTA file yielding FastaRecord objects.

    - Does NOT load whole file into memory.
    - Uppercases sequence.
    - Ignores empty lines and whitespace.
    """
    header: Optional[str] = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    yield FastaRecord(header=header, seq="".join(chunks).upper())
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        # final flush
        if header is not None:
            yield FastaRecord(header=header, seq="".join(chunks).upper())


def fasta_iter_string(text: str) -> Generator[FastaRecord, None, None]:
    """
    Alternative: parse from a string (for testing). Same behavior as parse_fasta.
    """
    header: Optional[str] = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield FastaRecord(header=header, seq="".join(chunks).upper())
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        yield FastaRecord(header=header, seq="".join(chunks).upper())


# ---------------------- METADATA & STATS ---------------------- #

def parse_primary_id(header: str) -> str:
    """
    Extract a primary ID from header (token up to first space).
    e.g., 'chr1 length=248956422' -> 'chr1'
    """
    return header.split()[0] if header else ""


def seq_stats(seq: str) -> Dict[str, float]:
    """
    Return basic sequence statistics:
    length, GC%, N_count, A/T/G/C counts.
    """
    n = len(seq)
    if n == 0:
        return {
            "length": 0,
            "GC%": 0.0,
            "N": 0,
            "A": 0, "T": 0, "G": 0, "C": 0
        }
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    n_count = seq.count("N")
    gc = (g + c) * 100.0 / n
    return {
        "length": n,
        "GC%": gc,
        "N": n_count,
        "A": a, "T": t, "G": g, "C": c
    }


def summarize_fasta(path: str, limit: int = 10) -> Dict[str, float]:
    """
    Print and return global summary:
    - number of sequences
    - total length
    - longest id/length
    - average GC%
    """
    total_len = 0
    total_gc_bases = 0
    count = 0
    longest_id = ""
    longest_len = -1

    print("\nFirst few records (up to", limit, "):")
    for i, rec in enumerate(parse_fasta(path), 1):
        sid = parse_primary_id(rec.header)
        st = seq_stats(rec.seq)
        if i <= limit:
            print(f"  {i:>3}. {sid} | len={st['length']} | GC%={st['GC%']:.2f} | N={st['N']}")
        # global
        count += 1
        total_len += st["length"]
        total_gc_bases += int(round(st["GC%"] * st["length"] / 100.0))
        if st["length"] > longest_len:
            longest_len = st["length"]
            longest_id = sid

    avg_gc = (total_gc_bases * 100.0 / total_len) if total_len else 0.0
    summary = {
        "num_sequences": count,
        "total_length": total_len,
        "avg_GC%": avg_gc,
        "longest_id": longest_id,
        "longest_len": longest_len if longest_len >= 0 else 0,
    }
    print("\nGlobal summary:", summary, "\n")
    return summary


# ---------------------- WRITE FASTA ---------------------- #

def wrap_seq(seq: str, width: int = FASTA_WRAP) -> Iterable[str]:
    """Yield wrapped lines of the sequence."""
    n = len(seq)
    for i in range(0, n, width):
        yield seq[i:i+width]


def write_fasta(records: Iterable[FastaRecord], path: str, width: int = FASTA_WRAP) -> None:
    """
    Write records to FASTA with standard wrapping.
    """
    with open(path, "w", encoding="utf-8") as out:
        for rec in records:
            out.write(">" + rec.header + "\n")
            for line in wrap_seq(rec.seq, width):
                out.write(line + "\n")


# ---------------------- FASTA INDEX (.fai-like) ---------------------- #
# We build a simple index file "<path>.fai" with one line per record:
# primary_id, header, seq_start_byte, seq_len
# seq_start_byte is the byte offset where sequence letters (no '>') begin,
# seq_len is total sequence length (letters only, no newlines).

@dataclass(slots=True)
class FastaIndexEntry:
    primary_id: str
    header: str
    seq_start: int  # byte offset in file where first base begins
    seq_len: int    # total bases (no newlines)


def build_fasta_index(path: str, index_path: Optional[str] = None) -> str:
    """
    Build a minimal index for random access of sequences (no per-line length info).
    Supports direct fetch of whole contig quickly; subrange fetch uses a linear scan
    from seq_start (still faster than full-file scan for large FASTA with few contigs).

    Returns the index path.
    """
    if index_path is None:
        index_path = path + ".fai"

    entries: list[FastaIndexEntry] = []
    with open(path, "rb") as fh:
        pos = 0
        current_header: Optional[str] = None
        seq_start = -1
        seq_len = 0

        def flush():
            nonlocal current_header, seq_start, seq_len
            if current_header is not None:
                primary = parse_primary_id(current_header)
                entries.append(FastaIndexEntry(primary, current_header, seq_start, seq_len))
            current_header = None
            seq_start = -1
            seq_len = 0

        while True:
            line = fh.readline()
            if not line:
                # EOF
                flush()
                break
            s = line.strip()
            if s.startswith(b">"):
                # header line
                flush()
                current_header = s[1:].decode("utf-8", errors="replace")
                seq_start = fh.tell()  # next byte after newline
            else:
                # sequence bytes (exclude newlines)
                seq_len += sum(1 for b in s if b not in (ord("\n"), ord("\r")))
            pos = fh.tell()

    # write index
    with open(index_path, "w", encoding="utf-8") as idx:
        for e in entries:
            idx.write(f"{e.primary_id}\t{e.header}\t{e.seq_start}\t{e.seq_len}\n")

    print(f"Index written: {index_path} ({len(entries)} entries)")
    return index_path


def read_fasta_index(index_path: str) -> Dict[str, FastaIndexEntry]:
    """Load index into memory."""
    table: Dict[str, FastaIndexEntry] = {}
    with open(index_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            raw = raw.strip()
            if not raw:
                continue
            primary, header, start_s, len_s = raw.split("\t")
            table[primary] = FastaIndexEntry(
                primary_id=primary,
                header=header,
                seq_start=int(start_s),
                seq_len=int(len_s),
            )
    return table


def fetch_full_sequence(path: str, entry: FastaIndexEntry) -> str:
    """
    Return the full sequence for an entry using the index.
    Reads line by line from seq_start, strips newlines, stops when length reached.
    """
    bases: list[str] = []
    need = entry.seq_len
    with open(path, "rb") as fh:
        fh.seek(entry.seq_start)
        while need > 0:
            chunk = fh.readline()
            if not chunk:
                break
            # remove whitespace/newlines
            letters = [chr(b) for b in chunk if b not in (ord("\n"), ord("\r"))]
            take = min(need, len(letters))
            if take:
                bases.extend(letters[:take])
                need -= take
    return "".join(bases)


def fetch_range(path: str, entry: FastaIndexEntry, start: int, end: int) -> str:
    """
    Fetch subrange [start, end) (0-based, end-exclusive) of a sequence using index.
    We still scan from seq_start but we stop early to avoid reading entire contig.
    For very large contigs, this is much faster than naive full-file parse.
    """
    if start < 0 or end > entry.seq_len or start >= end:
        return ""
    want_len = end - start
    got = 0
    pos = 0
    out: list[str] = []
    with open(path, "rb") as fh:
        fh.seek(entry.seq_start)
        for raw in fh:
            # strip newlines, whitespace
            letters = [chr(b) for b in raw if b not in (ord("\n"), ord("\r"))]
            if not letters:
                continue
            next_pos = pos + len(letters)
            # overlap with [start, end)?
            if next_pos > start and pos < end:
                s = max(0, start - pos)
                e = min(len(letters), end - pos)
                out.extend(letters[s:e])
                got += (e - s)
                if got >= want_len:
                    break
            pos = next_pos
            if pos >= end:
                break
    return "".join(out)


# ---------------------- CLI ---------------------- #

def _menu() -> None:
    print("\nFASTA Processor — CCA5 Q1")
    print("Options:")
    print(" 1) Inspect FASTA (first few records + global summary)")
    print(" 2) Write (rewrap) FASTA to a new file")
    print(" 3) Build index (.fai) for random access")
    print(" 4) Fetch FULL sequence by ID using index")
    print(" 5) Fetch RANGE by ID using index (start,end)")
    print(" 0) Exit\n")


def _cli():
    while True:
        _menu()
        choice = input("Enter 1/2/3/4/5/0: ").strip()

        if choice == "0":
            break

        elif choice == "1":
            path = input("FASTA path: ").strip()
            if not os.path.isfile(path):
                print("File not found.\n"); continue
            summarize_fasta(path)

        elif choice == "2":
            path = input("Source FASTA path: ").strip()
            if not os.path.isfile(path):
                print("File not found.\n"); continue
            outp = input("Output FASTA path (e.g., out.fa): ").strip() or "out.fa"
            width = int(input(f"Wrap width [{FASTA_WRAP}]: ").strip() or str(FASTA_WRAP))
            # stream parse → write
            recs = parse_fasta(path)
            write_fasta(recs, outp, width=width)
            print(f"Written: {outp}\n")

        elif choice == "3":
            path = input("FASTA path: ").strip()
            if not os.path.isfile(path):
                print("File not found.\n"); continue
            index_path = build_fasta_index(path)
            print("Index built:", index_path, "\n")

        elif choice == "4":
            path = input("FASTA path: ").strip()
            idx_path = input("Index path (.fai) [auto]: ").strip() or (path + ".fai")
            if not (os.path.isfile(path) and os.path.isfile(idx_path)):
                print("Paths not found.\n"); continue
            table = read_fasta_index(idx_path)
            key = input("Primary ID to fetch (e.g., chr1): ").strip()
            if key not in table:
                print("ID not in index.\n"); continue
            seq = fetch_full_sequence(path, table[key])
            print(f"Fetched length: {len(seq)}")
            show = min(120, len(seq))
            print("Preview:", seq[:show] + ("..." if len(seq) > show else ""), "\n")

        elif choice == "5":
            path = input("FASTA path: ").strip()
            idx_path = input("Index path (.fai) [auto]: ").strip() or (path + ".fai")
            if not (os.path.isfile(path) and os.path.isfile(idx_path)):
                print("Paths not found.\n"); continue
            table = read_fasta_index(idx_path)
            key = input("Primary ID to fetch (e.g., chr1): ").strip()
            if key not in table:
                print("ID not in index.\n"); continue
            try:
                start = int(input("Start (0-based): ").strip())
                end = int(input("End (exclusive): ").strip())
            except ValueError:
                print("Invalid start/end.\n"); continue
            sub = fetch_range(path, table[key], start, end)
            print(f"Fetched subrange length: {len(sub)}")
            show = min(120, len(sub))
            print("Preview:", sub[:show] + ("..." if len(sub) > show else ""), "\n")

        else:
            print("Please enter a valid option.\n")


if __name__ == "__main__":
    _cli()
