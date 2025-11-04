"""
advanced_patterns.py — CCA5 · Part B · Q4
- Suffix array (naive O(n log n)) + substring search
- Find repeated substrings (max repeat)
- Palindromic substrings (simple expand-around-center)
- Tiny assembly simulation with k-mer overlaps (greedy)
- UPGMA clustering from a distance matrix (phylogeny demo)
"""

from __future__ import annotations
from typing import List, Tuple, Dict
from collections import defaultdict
import math

# --- suffix array ---
def suffix_array(s: str) -> List[int]:
    return sorted(range(len(s)), key=lambda i: s[i:])

def sa_search(text: str, pattern: str) -> List[int]:
    sa = suffix_array(text)
    lo, hi = 0, len(sa)
    out = []
    # binary search lower bound
    while lo < hi:
        mid = (lo+hi)//2
        if text[sa[mid]:].startswith(pattern) or text[sa[mid]:] >= pattern:
            hi = mid
        else:
            lo = mid+1
    start = lo
    while start < len(sa) and text[sa[start]:].startswith(pattern):
        out.append(sa[start]); start += 1
    return sorted(out)

# --- repeats ---
def longest_repeat(s: str) -> str:
    sa = suffix_array(s)
    best = ""
    for i in range(1, len(sa)):
        a, b = sa[i-1], sa[i]
        l = 0
        while a+l < len(s) and b+l < len(s) and s[a+l]==s[b+l]:
            l += 1
        if l > len(best):
            best = s[a:a+l]
    return best

# --- palindromes ---
def palindromes(s: str, min_len: int = 3) -> List[Tuple[int,int,str]]:
    out: List[Tuple[int,int,str]] = []
    def expand(l: int, r: int):
        while l>=0 and r<len(s) and s[l]==s[r]:
            if r-l+1>=min_len:
                out.append((l,r+1,s[l:r+1]))
            l-=1; r+=1
    for i in range(len(s)):
        expand(i,i)      # odd
        expand(i,i+1)    # even
    return out

# --- assembly simulation (greedy overlap) ---
def assemble_reads(reads: List[str], k: int = 20) -> str:
    reads = reads[:]
    if not reads: return ""
    while len(reads) > 1:
        best_i = best_j = -1
        best_ov = -1
        # find best overlap suffix(p) with prefix(q) of length >= k
        for i, p in enumerate(reads):
            for j, q in enumerate(reads):
                if i==j: continue
                max_ov = min(len(p), len(q))
                ov = 0
                for L in range(max_ov, k-1, -1):
                    if p.endswith(q[:L]):
                        ov = L; break
                if ov > best_ov:
                    best_ov = ov; best_i, best_j = i, j
        if best_ov < k:
            # no strong overlap; just concatenate
            reads[0] += reads.pop()
        else:
            reads[best_i] = reads[best_i] + reads[best_j][best_ov:]
            reads.pop(best_j)
    return reads[0]

# --- UPGMA clustering ---
def upgma(names: List[str], dist: List[List[float]]) -> Tuple[str, Dict[str,float]]:
    """
    Very small UPGMA implementation returning a Newick string and heights.
    dist: symmetric matrix (names x names), zeros on diagonal.
    """
    clusters = {i: names[i] for i in range(len(names))}
    heights: Dict[int,float] = {i:0.0 for i in range(len(names))}
    active = set(range(len(names)))
    def d(i,j): return dist[i][j]
    while len(active) > 1:
        pair, best = None, float("inf")
        al = list(active)
        for a_i in range(len(al)):
            for b_i in range(a_i+1, len(al)):
                i, j = al[a_i], al[b_i]
                if d(i,j) < best:
                    best, pair = d(i,j), (i,j)
        i, j = pair
        # new cluster index
        k = max(clusters.keys()) + 1
        hi = best/2
        newick = f"({clusters[i]}:{hi - heights[i]:.3f},{clusters[j]}:{hi - heights[j]:.3f})"
        clusters[k] = newick
        heights[k] = hi
        # update distances (average)
        newrow = {}
        for t in list(active):
            if t in (i,j): continue
            newd = (d(i,t) + d(j,t))/2
            newrow[t] = newd
        # expand dist matrix
        for t in newrow:
            # ensure dist lists large enough
            while len(dist) <= k: dist.append([0.0]*len(dist))
            for row in dist:
                while len(row) <= k: row.append(0.0)
            dist[k][t] = dist[t][k] = newrow[t]
        active.remove(i); active.remove(j); active.add(k)
    root = clusters[list(active)[0]]
    return root+";", {names[i]:heights[i] for i in range(len(names))}

# -------- CLI --------
def _menu():
    print("\nAdvanced Patterns — Q4")
    print(" 1) Suffix-array substring search")
    print(" 2) Longest repeated substring")
    print(" 3) Palindromic substrings")
    print(" 4) Greedy assembly simulation (overlaps)")
    print(" 5) UPGMA clustering (Newick)")
    print(" 0) Exit\n")

if __name__ == "__main__":
    while True:
        _menu()
        ch = input("Enter 1/2/3/4/5/0: ").strip()
        if ch == "0": break
        elif ch == "1":
            text = input("Text: ").strip().upper()
            pat  = input("Pattern: ").strip().upper()
            print("Indices:", sa_search(text, pat), "\n")
        elif ch == "2":
            s = input("Sequence: ").strip().upper()
            print("Longest repeat:", longest_repeat(s), "\n")
        elif ch == "3":
            s = input("Sequence: ").strip().upper()
            m = int(input("Min palindrome length [3]: ") or "3")
            pals = palindromes(s, m)
            print("Count:", len(pals)); 
            for p in pals[:10]: print(p)
            print()
        elif ch == "4":
            n = int(input("How many reads [4]: ") or "4")
            reads = [input(f"read{i+1}: ").strip().upper() for i in range(n)]
            k = int(input("min overlap k [3]: ") or "3")
            print("Assembly:", assemble_reads(reads, k), "\n")
        elif ch == "5":
            n = int(input("How many taxa [4]: ") or "4")
            names = [input(f"name{i+1}: ").strip() for i in range(n)]
            print("Enter upper triangle distances row-wise (excluding diagonal).")
            dist = [[0.0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1,n):
                    val = float(input(f"d({names[i]},{names[j]}): "))
                    dist[i][j]=dist[j][i]=val
            newick, heights = upgma(names, dist)
            print("Newick:", newick)
            print("Heights:", heights, "\n")
        else:
            print("Invalid.\n")
