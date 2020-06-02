"""
Utility script for easy access to most of the VEnCode scripts and tests.
Example on how to use in command-line:
>python process.py run_tests
"""

import sys
sys.path.append("..")

import argparse
from difflib import get_close_matches

from VEnCode import tests


def run_tests():
    tests.run_all_tests()


def main(process, possibilities):
    close_matches = get_close_matches(process, possibilities, n=3, cutoff=0.6)
    if process not in close_matches:
        if close_matches:
            print(f"Process {process} does ot exist. Did you mean one of the following:\n", "\n ".join(close_matches))
        else:
            print(f"Process {process} does ot exist.")
        return
    eval(f"{process}()")


if __name__ == "__main__":
    list_of_processes = ["run_tests"]
    parser = argparse.ArgumentParser(description="""Utility script for easy access to most of the VEnCode scripts and
    tests.
    Example on how to use in command-line:
    >python process.py run_tests
    """)
    parser.add_argument("process", help="The process/script to run.")

    args = parser.parse_args()
    main(args.process, list_of_processes)
