"""
data_integration.py — CCA5 · Part A · Q2
- Tiny SQLite DB for sequences (id, header, seq, source)
- Import from FASTA (using fasta_processor)
- FASTQ basics: parse identifiers, sequences, qualities; compute mean Q
- Simple QC rules + search/retrieve APIs
"""

from __future__ import annotations
import os, sqlite3
from typing import Iterator, Tuple, Dict, List
from fasta_processor import parse_fasta, FastaRecord, parse_primary_id

DB_PATH = "sequences.db"

# ---------- FASTQ basics ----------
def parse_fastq(path: str) -> Iterator[Tuple[str,str,List[int]]]:
    """
    Yield (id, sequence, qualities[]) assuming Sanger FASTQ with ASCII+33.
    """
    with open(path, "r", encoding="utf-8") as fh:
        while True:
            h = fh.readline()
            if not h: break
            if not h.startswith("@"): continue
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()
            qs = [ord(ch)-33 for ch in qual]
            yield (h[1:].strip(), seq.upper(), qs)

def mean_quality(quals: List[int]) -> float:
    return sum(quals) / (len(quals) or 1)

# ---------- SQLite schema ----------
SCHEMA = """
CREATE TABLE IF NOT EXISTS sequences (
  id TEXT PRIMARY KEY,
  header TEXT,
  seq TEXT NOT NULL,
  source TEXT
);
CREATE INDEX IF NOT EXISTS idx_source ON sequences(source);
"""

def db_open(path: str = DB_PATH) -> sqlite3.Connection:
    con = sqlite3.connect(path)
    con.executescript(SCHEMA)
    return con

def db_insert_fasta(con: sqlite3.Connection, path: str, source: str = "FASTA") -> int:
    cur = con.cursor()
    n = 0
    for rec in parse_fasta(path):
        sid = parse_primary_id(rec.header)
        cur.execute("INSERT OR REPLACE INTO sequences(id, header, seq, source) VALUES (?,?,?,?)",
                    (sid, rec.header, rec.seq, source))
        n += 1
    con.commit()
    return n

def db_insert_fastq(con: sqlite3.Connection, path: str, qmin: int = 10, qmean: float = 20.0,
                    source: str = "FASTQ") -> int:
    """
    Import FASTQ reads that pass a simple QC: all bases ATGCN, mean Q >= qmean, each qual >= qmin.
    """
    cur = con.cursor()
    n = 0
    for rid, seq, qs in parse_fastq(path):
        if any(ch not in "ATGCN" for ch in seq): 
            continue
        if min(qs or [0]) < qmin: 
            continue
        if mean_quality(qs) < qmean:
            continue
        cur.execute("INSERT OR REPLACE INTO sequences(id, header, seq, source) VALUES (?,?,?,?)",
                    (rid, rid, seq, source))
        n += 1
    con.commit()
    return n

def db_get(con: sqlite3.Connection, seq_id: str) -> Dict[str,str] | None:
    cur = con.execute("SELECT id, header, seq, source FROM sequences WHERE id=?", (seq_id,))
    row = cur.fetchone()
    if not row: return None
    return {"id": row[0], "header": row[1], "seq": row[2], "source": row[3]}

def db_search_prefix(con: sqlite3.Connection, prefix: str, limit: int = 10) -> List[Dict[str,str]]:
    cur = con.execute("SELECT id, header, length(seq) FROM sequences WHERE id LIKE ? ORDER BY id LIMIT ?",
                      (prefix + "%", limit))
    return [{"id": r[0], "header": r[1], "length": r[2]} for r in cur.fetchall()]

# --------- CLI ----------
def _menu():
    print("\nData Integration — Q2")
    print(" 1) Import FASTA to SQLite")
    print(" 2) Import FASTQ (with QC) to SQLite")
    print(" 3) Fetch by ID")
    print(" 4) Search by ID prefix")
    print(" 0) Exit\n")

if __name__ == "__main__":
    con = db_open(DB_PATH)
    try:
        while True:
            _menu()
            ch = input("Enter 1/2/3/4/0: ").strip()
            if ch == "0": break
            elif ch == "1":
                p = input("FASTA path: ").strip()
                if not os.path.isfile(p): print("File not found.\n"); continue
                n = db_insert_fasta(con, p)
                print(f"Imported {n} sequences.\n")
            elif ch == "2":
                p = input("FASTQ path: ").strip()
                if not os.path.isfile(p): print("File not found.\n"); continue
                qmin = int(input("min per-base Q [10]: ") or "10")
                qmean = float(input("min mean Q [20]: ") or "20")
                n = db_insert_fastq(con, p, qmin=qmin, qmean=qmean)
                print(f"Imported {n} reads passing QC.\n")
            elif ch == "3":
                sid = input("ID: ").strip()
                print(db_get(con, sid), "\n")
            elif ch == "4":
                pref = input("Prefix: ").strip()
                print(db_search_prefix(con, pref), "\n")
            else:
                print("Invalid.\n")
    finally:
        con.close()
