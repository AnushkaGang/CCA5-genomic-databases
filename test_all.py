"""
tests/test_all.py — CCA5 quick functional tests

Run:
  python tests/test_all.py
"""

from fasta_processor import FastaRecord, write_fasta, parse_fasta
from data_integration import mean_quality
from alignment_basics import needleman_wunsch, lcs_multi, jukes_cantor, consensus
from advanced_patterns import sa_search, longest_repeat, palindromes, assemble_reads
from performance_tools import gc_percent

def main():
    print("\n--- CCA5 BASIC TESTS ---\n")

    # FASTA write+parse roundtrip
    recs = [FastaRecord("seq1", "ATGCGC"), FastaRecord("seq2", "TTTTAAAACCC")]
    write_fasta(recs, "tmp.fa", 5)
    parsed = list(parse_fasta("tmp.fa"))
    assert parsed[0].seq == "ATGCGC"

    # FASTQ mean quality
    print("meanQ:", mean_quality([40,40,30]))

    # Alignment + LCS + JC distance
    a,b,score = needleman_wunsch("GATTACA", "GCATGCU")
    print("NW score:", score)
    print("LCS:", lcs_multi(["GATTACA","GCATGCU","GCTACA"]))
    print("JC:", jukes_cantor("AAAA","AAAT"))

    # Consensus (aligned)
    print("Consensus:", consensus(["A-CG","ATCG","AT-G"]))

    # Patterns
    print("SA find:", sa_search("ATGCATGCAT", "CAT"))
    print("Longest repeat:", longest_repeat("banana"))
    print("Palindromes:", palindromes("ATGCGAATTC", 4)[:3])
    print("Assembly:", assemble_reads(["ATGCG","GCGTT","TTAAA"], k=2))

    # GC%
    print("GC%:", gc_percent("ATGCGC"))

    print("\n--- TESTS DONE ---\n")

if __name__ == "__main__":
    main()
