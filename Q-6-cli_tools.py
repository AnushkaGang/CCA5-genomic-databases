"""
cli_tools.py — CCA5 · Part C · Q6
Unified command-line entry point + logging & error handling.

Examples:
  python cli_tools.py gc --seq ATGCGC
  python cli_tools.py fasta-summary --path genome.fa
"""

from __future__ import annotations
import argparse, logging, sys
from fasta_processor import summarize_fasta
from performance_tools import gc_percent

def make_logger(level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("cca5")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    if not logger.handlers:
        h = logging.StreamHandler(sys.stdout)
        fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger

def main():
    parser = argparse.ArgumentParser(description="CCA5 unified CLI")
    parser.add_argument("--log", default="INFO", help="log level (INFO/DEBUG/WARN)")
    sub = parser.add_subparsers(dest="cmd")

    p1 = sub.add_parser("gc", help="compute GC% of one sequence")
    p1.add_argument("--seq", required=True)

    p2 = sub.add_parser("fasta-summary", help="summarize a FASTA file")
    p2.add_argument("--path", required=True)

    args = parser.parse_args()
    log = make_logger(args.log)

    try:
        if args.cmd == "gc":
            val = gc_percent(args.seq)
            log.info("GC%% = %.2f", val)
        elif args.cmd == "fasta-summary":
            summarize_fasta(args.path)
        else:
            parser.print_help()
    except Exception as e:
        log.error("Error: %s", e)

if __name__ == "__main__":
    main()
