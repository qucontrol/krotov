#!/usr/bin/env python
"""Local pre-commit hooks"""
import os
import re
import subprocess
import sys
from argparse import ArgumentParser, RawTextHelpFormatter


try:
    from tox.session import load_config
except ImportError:
    print("tox must be available for pre-commit hooks.")
    print("See https://tox.readthedocs.io for installation instructions.")
    sys.exit(1)


def no_trailing_whitespace(filenames):
    """Check that files have no trailing whitespace."""
    success = True
    for filename in filenames:
        with open(filename) as in_fh:
            for (line_index, line) in enumerate(in_fh):
                if line.endswith(" \n"):
                    print(
                        "%s:%d has trailing whitespace"
                        % (filename, line_index + 1)
                    )
                    success = False
    return success


def no_debug_comment(filenames):
    """Check that files have to DEBUG comments."""
    success = True
    for filename in filenames:
        print("Checking %s for debug comments" % filename)
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


def run_tox_env_commands(envname, *args):
    """Run tox-environment commands as a pre-commit hook.

    Roughly equivalent to checking that ``tox -e <envname> -- <args>`` finishes
    with exit status 0. However, the ``tox`` executable is not invoked.
    Instead, the environment directory for the given `envname` is extracted
    from the default tox configuration file, as well as the list of commands
    for that environment. The commands are then executed directly with a PATH
    variable set to the bin-folder of the environment directory. This means
    that the environment directory must exist (from a previous manual
    invocation of tox). This avoids potential "hangs" and other unforeseen
    problems from having to initialize a tox environment during the execution
    of a pre-commit hook (which should be relatively fast).

    Returns:
        bool: whether all the commands defined in the default tox configuration
        file for the given `envname` complete with exit code 0.

    Raises:
        OSError: If the tox environment directory for `envname` does not
        already exist.
    """
    # Besides speed, another reason we are not invoking `tox` directly is to
    # allow for the possibility that a tox environment was previously created
    # from an alternative tox.ini file
    tox_cmdline = ['-e', envname, '--', *args]
    config = load_config(tox_cmdline)
    try:
        env_config = config.envconfigs[envname]
    except KeyError:
        print("tox is not set up correctly. Run `tox -e bootstrap`")
        sys.exit(1)
    if env_config.envdir.isdir():
        success = True
        shell_env = {'PATH': str(env_config.envbindir)}
        if 'PATH' in os.environ:
            shell_env['PATH'] += os.pathsep + os.environ['PATH']
        for varname in env_config.passenv:
            if varname == 'PATH':
                continue
            else:
                if varname in os.environ:
                    shell_env[varname] = os.environ[varname]
        for cmd in env_config.commands:
            proc = subprocess.run(cmd, env=shell_env, cwd=env_config.changedir)
            success &= proc.returncode == 0
        return success
    else:
        raise OSError(
            "The tox environment is not initialized. Run 'tox -e %s' manually "
            "once before trying to commit!" % envname
        )


def black(filenames):
    """Check that files conform to Black formatting.

    This uses the `run-blackcheck` tox environment to ensure that the results
    of the pre-commit hook and a manual invocation of `tox -e run-blackcheck`
    are identical.
    """
    return run_tox_env_commands('run-blackcheck', *filenames)


def isort(filenames):
    """Check that files conform to isort formatting.

    This uses the `run-isortcheck` tox environment to ensure that the results
    of the pre-commit hook and a manual invocation of `tox -e run-isortcheck`
    are identical.
    """
    return run_tox_env_commands('run-isortcheck', *filenames)


CHECKS = {
    'whitespace': no_trailing_whitespace,
    'debug-comments': no_debug_comment,
    'black': black,
    'isort': isort,
}


def main(argv=None):
    """Main function"""
    description = "Perform the given CHECK on the given FILENAMES.\n\n"
    description += "The following CHECKs are available:\n\n"
    for (name, check) in CHECKS.items():
        description += "    " + name + ":\n"
        description += "        " + check.__doc__.splitlines()[0] + "\n"
    parser = ArgumentParser(
        description=description, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument('CHECK', help='Name of check')
    parser.add_argument('FILENAMES', nargs='*', help='Filenames to check')
    args = parser.parse_args(argv)

    return_code = 0
    success = True

    check = CHECKS[args.CHECK]
    try:
        success = check(args.FILENAMES)
    except (ValueError, OSError) as exc_info:
        print(
            "Cannot check %s with %s: %s"
            % (args.FILENAMES, args.CHECK, exc_info)
        )
        return_code = 1
    if not success:
        return_code = 1
    return return_code


if __name__ == '__main__':
    sys.exit(main())
