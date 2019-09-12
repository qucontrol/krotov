#!/usr/bin/env python
"""Stand-in for pip that processes packages sequentially.

`pip install <packages>` compiles *all* the given packages before installing
them. This can be a problem if the compilation of one package depends on other
packages being installed (most likely, cython/numpy). This script provides an
ad-hoc solution by translating `pip install <packages` into a sequential `pip
install <package>` for every package in <packages>. It can be used in tox.ini
as

    install_command=
        python scripts/pip_sequential.py install {opts} -- {packages}
"""

import subprocess
import sys


def main(argv=None):
    """Main function"""
    if argv is None:
        argv = sys.argv
    command = 'help'
    options = []
    args = []
    if len(argv) > 1 and not argv[1].startswith('-'):
        command = argv[1]
        bucket = options
        for arg in argv[2:]:
            if arg == '--':
                # everything before '--' is definitely an option, everything
                # afterwards *may* be an arg
                bucket = args
            else:
                if arg.startswith('-'):
                    options.append(arg)
                else:
                    bucket.append(arg)
    if len(args) == 0:
        print("Usage: %s command [options] -- <specs>" % __file__)
        return 1
    try:
        for arg in args:
            cmd = [sys.executable, '-m', 'pip', command, *options, arg]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc_info:
        print("ERROR: %s" % exc_info)
        return 1
    else:
        return 0


if __name__ == '__main__':
    sys.exit(main())
