"""Usage: python -m nanosim.run_case case.txt out.csv"""
import sys
from .core.case import read_case, run_case, write_csv

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python -m nanosim.run_case case.txt out.csv")
    write_csv(sys.argv[2], run_case(read_case(sys.argv[1])))
