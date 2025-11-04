"""
alignment_basics.py — CCA5 · Part B · Q3
- Needleman–Wunsch global alignment (simple scoring)
- Hamming / p-distance / Jukes-Cantor distance
- LCS for multiple sequences (pairwise-iterative)
- Consensus sequence from equal-length alignment
"""

from __future__ import annotations
from typing import Tuple, List, Dict
import math
from collections import Counter

def needleman_wunsch(a: str, b: str, match: int=1, mismatch: int=-1, gap: int=-1) -> Tuple[str,str,int]:
    n, m = len(a), len(b)
    dp = [[0]*(m+1) for _ in range(n+1)]
    bt = [[None]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0] = i*gap; bt[i][0] = "U"
    for j in range(1,m+1): dp[0][j] = j*gap; bt[0][j] = "L"
    for i in range(1,n+1):
        for j in range(1,m+1):
            diag = dp[i-1][j-1] + (match if a[i-1]==b[j-1] else mismatch)
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            best = max(diag, up, left)
            dp[i][j] = best
            bt[i][j] = "D" if best==diag else ("U" if best==up else "L")
    # traceback
    i, j = n, m
    al_a, al_b = [], []
    while i>0 or j>0:
        if bt[i][j] == "D":
            al_a.append(a[i-1]); al_b.append(b[j-1]); i-=1; j-=1
        elif bt[i][j] == "U":
            al_a.append(a[i-1]); al_b.append("-"); i-=1
        else:
            al_a.append("-"); al_b.append(b[j-1]); j-=1
    return "".join(reversed(al_a)), "".join(reversed(al_b)), dp[n][m]

def hamming(a: str, b: str) -> int:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i]!=b[i]) + abs(len(a)-len(b))

def p_distance(a: str, b: str) -> float:
    n = max(len(a), len(b)) or 1
    return hamming(a,b)/n

def jukes_cantor(a: str, b: str) -> float:
    p = p_distance(a,b)
    if p >= 0.75:  # avoid log of negative
        return float("inf")
    try:
        return -3/4 * math.log(1 - (4/3)*p)
    except ValueError:
        return float("inf")

def lcs_two(a: str, b: str) -> str:
    n, m = len(a), len(b)
    dp = [[""]*(m+1) for _ in range(n+1)]
    for i in range(n):
        for j in range(m):
            dp[i+1][j+1] = (dp[i][j] + a[i]) if a[i]==b[j] else (dp[i][j+1] if len(dp[i][j+1])>=len(dp[i+1][j]) else dp[i+1][j])
    return dp[n][m]

def lcs_multi(seqs: List[str]) -> str:
    if not seqs: return ""
    cur = seqs[0]
    for s in seqs[1:]:
        cur = lcs_two(cur, s)
        if not cur: break
    return cur

def consensus(aligned: List[str]) -> str:
    if not aligned: return ""
    L = len(aligned[0])
    for s in aligned:
        if len(s)!=L:
            raise ValueError("All aligned strings must be same length.")
    out = []
    for col in zip(*aligned):
        base = Counter(col).most_common(1)[0][0]
        out.append(base)
    return "".join(out)

# -------- CLI --------
def _menu():
    print("\nAlignment Basics — Q3")
    print(" 1) Needleman–Wunsch (global)")
    print(" 2) Distances (Hamming / p / JC)")
    print(" 3) LCS for multiple sequences")
    print(" 4) Consensus from aligned strings")
    print(" 0) Exit\n")

if __name__ == "__main__":
    while True:
        _menu()
        ch = input("Enter 1/2/3/4/0: ").strip()
        if ch == "0": break
        elif ch == "1":
            a = input("Seq A: ").strip().upper()
            b = input("Seq B: ").strip().upper()
            al_a, al_b, score = needleman_wunsch(a,b)
            print("A:", al_a); print("B:", al_b); print("Score:", score, "\n")
        elif ch == "2":
            a = input("Seq A: ").strip().upper()
            b = input("Seq B: ").strip().upper()
            print("Hamming:", hamming(a,b))
            print("p-distance:", p_distance(a,b))
            print("Jukes-Cantor:", jukes_cantor(a,b), "\n")
        elif ch == "3":
            n = int(input("How many sequences [3]: ") or "3")
            seqs = [input(f"S{i+1}: ").strip().upper() for i in range(n)]
            print("LCS:", lcs_multi(seqs), "\n")
        elif ch == "4":
            n = int(input("How many aligned strings [3]: ") or "3")
            arr = [input(f"A{i+1}: ").strip() for i in range(n)]
            print("Consensus:", consensus(arr), "\n")
        else:
            print("Invalid.\n")
