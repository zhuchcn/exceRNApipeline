from .unlock import parse_args as parse_args_unlock
from .clean_scratch import parse_args as parse_args_clean_scratch
from .run import parse_args as parse_args_run
from .init import parse_args as parse_args_init
import subprocess as sp
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="commands", dest="command")
    parse_args_unlock(subparsers)
    parse_args_clean_scratch(subparsers)
    parse_args_init(subparsers)
    parse_args_run(subparsers)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()