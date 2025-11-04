"""
performance_tools.py — CCA5 · Part C · Q5
- Optimize simple functions (counting GC with fast paths)
- Parallel processing over many sequences (multiprocessing)
- Memory-friendly generators
- Text progress bar for long runs
"""

from __future__ import annotations
from typing import Iterable, List, Tuple
import time, multiprocessing as mp, sys

def gc_percent(seq: str) -> float:
    s = seq.upper()
    n = len(s) or 1
    return (s.count("G") + s.count("C"))*100.0/n

def gc_many_parallel(seqs: Iterable[str], workers: int = None) -> List[float]:
    with mp.Pool(processes=workers) as pool:
        return list(pool.map(gc_percent, seqs))

def stream_file_lines(path: str) -> Iterable[str]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            yield line.strip()

def progress(iterable: Iterable, total: int|None=None, width: int=30):
    count = 0
    t0 = time.time()
    for item in iterable:
        count += 1
        if total:
            done = int(width * count / total)
            bar = "[" + "#"*done + "-"*(width-done) + "]"
            elapsed = time.time()-t0
            rate = count / (elapsed or 1e-9)
            sys.stdout.write(f"\r{bar} {count}/{total} {rate:.1f}/s")
            sys.stdout.flush()
        yield item
    if total:
        sys.stdout.write("\n"); sys.stdout.flush()

# -------- CLI --------
if __name__ == "__main__":
    print("\nPerformance Tools — Q5")
    print(" 1) GC% of many sequences in parallel")
    print(" 2) Stream file lines with progress")
    print(" 0) Exit\n")
    while True:
        ch = input("Enter 1/2/0: ").strip()
        if ch == "0": break
        elif ch == "1":
            n = int(input("How many sequences [5]: ") or "5")
            seqs = [input(f"S{i+1}: ").strip() for i in range(n)]
            res = gc_many_parallel(seqs)
            print("GC%s:", res, "\n")
        elif ch == "2":
            p = input("Text file path: ").strip()
            try:
                lines = list(progress(stream_file_lines(p), total=None))
                print(f"\nRead {len(lines)} lines.\n")
            except FileNotFoundError:
                print("File not found.\n")
        else:
            print("Invalid.\n")
