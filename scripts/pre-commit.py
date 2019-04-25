#!/usr/bin/env python
"""Lokal pre-commit hook

- no trailing whitespace in any file
- no lines with a "DEBUG" comment in a python file

This can be extended with project-specific checks. You may also consider
third-party hooks available from https://pre-commit.com/hooks.html
"""
import argparse
import re
import sys


def no_trailing_whitespace(filename):
    "Check that file has not trailing whitespace"
    success = True
    with open(filename) as in_fh:
        for (line_index, line) in enumerate(in_fh):
            if line.endswith(" \n"):
                print(
                    "%s:%d has trailing whitespace"
                    % (filename, line_index + 1)
                )
                success = False
    return success


def no_debug_comment(filename):
    "Check that file has not trailing whitespace"
    success = True
    if filename.endswith(".py"):
        rx_debug = re.compile(r'#\s*DEBUG')
        with open(filename) as in_fh:
            for (line_index, line) in enumerate(in_fh):
                if rx_debug.search(line):
                    print(
                        "%s:%d has a DEBUG marker comment"
                        % (filename, line_index + 1)
                    )
                    success = False
    return success


def main(argv=None):
    """Main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='*', help='Filenames to fix')
    args = parser.parse_args(argv)

    return_code = 0
    success = True
    for filename in args.filenames:
        success &= no_trailing_whitespace(filename)
        success &= no_debug_comment(filename)
    if not success:
        return_code = 1
    return return_code


if __name__ == '__main__':
    sys.exit(main())
