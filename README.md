# CCA5 – Genomic Databases & Advanced Applications

Assignment 3 implementation in Python. Each module has a simple **user-input CLI**.

## Student
**Name:** Anushka Gangwar  
**PRN:** 1032233324  
**Course:** TY BTech CSE Panel A — Roll No 46

---

## Modules / Questions

| File | Question | Topic |
|------|----------|-------|
| fasta_processor.py | Part A Q1 | FASTA parsing, metadata, writing, indexing |
| data_integration.py | Part A Q2 | Mini SQLite DB, FASTQ basics & QC, search |
| alignment_basics.py | Part B Q3 | Global alignment, LCS, distances, consensus |
| advanced_patterns.py | Part B Q4 | Suffix array search, repeats, palindromes, assembly, UPGMA |
| performance_tools.py | Part C Q5 | Parallel processing, streaming, progress bar |
| cli_tools.py | Part C Q6 | Unified CLI with argparse + logging |
| tests/test_all.py | — | Quick functional tests |

---

## How to Run (examples)

```bash
# Individual modules (interactive):
python fasta_processor.py
python data_integration.py
python alignment_basics.py
python advanced_patterns.py
python performance_tools.py

# Unified CLI:
python cli_tools.py gc --seq ATGCGC
python cli_tools.py fasta-summary --path genome.fa

# Tests:
python tests/test_all.py


Performance Notes

Streaming FASTA parsing and range fetch (index) are O(n) time, O(1) extra memory.

Global alignment is classic DP (O(nm)).

Suffix array implemented with naive sort (O(n log n) build) — OK for coursework scale.

Parallel map uses Python multiprocessing.

Documentation

All modules include docstrings and clear prompts with example input.


