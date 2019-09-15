#!/usr/bin/env python
"""Bootstrap script for setting up tox and pre-commit hooks

This scripts is called either by any invocation of the Makefile, or by the tox
configuration in setup.cvg, if there is no tox.ini.

It sets up the main tox.ini file and the pre-commit hooks.
"""
import os
import pathlib
import shutil
import sys


_TOX_ERR_MSG = r'''
tox is not available. See https://tox.readthedocs.io for installation
instructions.
'''

_TOX_OUTDATED_WARN = r'''
WARNING: tox.ini is older than tox-pyenv.ini and tox-conda.ini. Touch or modify
tox.ini to confirm.
'''


def main(argv=None):
    """Main function"""
    if argv is None:
        argv = sys.argv
    root = pathlib.Path(__file__).parent.parent

    # 1. Bootstrap the tox.ini file
    if (root / 'tox.ini').is_file():
        # if there are updates to the tracked tox-pyenv.ini and tox-conda.ini,
        # people might forget that this doesn't update their local,
        # bootstrapped tox.ini. This remindes them.
        tox_ini_mtime = (root / 'tox.ini').stat().st_mtime
        tox_pyenv_mtime = (root / 'tox-pyenv.ini').stat().st_mtime
        tox_conda_mtime = (root / 'tox-conda.ini').stat().st_mtime
        if tox_ini_mtime < max(tox_pyenv_mtime, tox_conda_mtime):
            print(_TOX_OUTDATED_WARN)
    else:
        tox_ini = os.environ.get('TOXINI', 'tox-pyenv.ini')
        if not (root / tox_ini).is_file():
            tox_ini = 'tox-pyenv.ini'
        tox_src = root / tox_ini
        tox_dst = root / "tox.ini"
        print("Copying %s to %s" % (tox_src, tox_dst))
        shutil.copyfile(tox_src, tox_dst)

    # 2. Ensure tox is installed
    try:
        import tox
    except ImportError:
        print(_TOX_ERR_MSG)
        sys.exit(1)

    # 3. Ensure pre-commit hooks are installed
    if not (root / ".git" / "hooks" / "pre-commit").is_file():
        # we're using the tox.ini environments that we got from step 1
        print("bootstrapping pre-commit hook")
        cmdline = ['-e', 'run-cmd', '--', 'pre-commit', 'install']
        print("tox " + " ".join(cmdline))
        tox.cmdline(cmdline)


if __name__ == '__main__':
    sys.exit(main())
