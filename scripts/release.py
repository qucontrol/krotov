#!/usr/bin/env python
"""Automation script for making a release. Must be run from the root for the
repository"""
# Note: Version scheme according to https://www.python.org/dev/peps/pep-0440
import os
from os.path import join
import sys
from subprocess import run, DEVNULL, CalledProcessError
from pkg_resources import parse_version
import json
import urllib.request
import urllib.error
import urllib.parse
import shutil

import git
import click


def make_release(package_name):
    """Interactively create and publish a new release for the package"""
    click.confirm("Do you want to make a release?", abort=True)
    check_git_clean(git.Repo(os.getcwd()))
    run_tests()
    new_version = ask_for_release_version(package_name)
    set_version(join('.', 'src', package_name, '__init__.py'), new_version)
    edit_history(new_version)
    while not check_dist():
        click.confirm(
            "Fix errors manually! Continue?", default=True, abort=True
        )
    check_docs()
    make_release_commit(new_version)
    make_upload(test=True)
    push_release_commit()
    make_upload(test=False)
    make_and_push_tag(new_version)
    next_dev_version = get_next_dev_version(new_version)
    set_version(
        join('.', 'src', package_name, '__init__.py'), next_dev_version
    )
    make_next_dev_version_commit(next_dev_version)


class ReleaseError(ValueError):
    pass


def pypi_versions(package_name):
    """Return list of versions for the given package on PyPI"""
    url = "https://pypi.python.org/pypi/%s/json" % (package_name,)
    data = json.load(urllib.request.urlopen(urllib.request.Request(url)))
    versions = list(data["releases"].keys())
    versions.sort(key=parse_version)
    return versions


def get_version(filename):
    """Extract the package version"""
    with open(filename) as in_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                return parse_version(line.split('=')[1].strip()[1:-1])
    raise ReleaseError("Cannot extract version from %s" % filename)


def edit(filename):
    """Open filename in EDITOR"""
    editor = os.getenv('EDITOR', 'vi')
    if click.confirm("Open %s in %s?" % (filename, editor), default=True):
        run([editor, filename])


def check_git_clean(repo):
    """Ensure that a given git.Repo is clean"""
    if repo.active_branch.name != 'master':
        raise ReleaseError("Repository must be on the master branch")
    if repo.is_dirty():
        run(['git', 'status'])
        raise ReleaseError("Repository must be in a clean state")
    if repo.untracked_files:
        click.echo("WARNING: there are untracked files:")
        for filename in repo.untracked_files:
            click.echo("\t%s" % filename)
        click.confirm("Continue?", default=False, abort=True)


def run_tests():
    """Run 'make test'"""
    run(['make', 'test'], check=True)


def ask_for_release_version(package_name):
    """Ask for the version number of the release.

    The version number must be greater than the last PyPI release
    """
    current_version = get_version(
        join('.', 'src', package_name, '__init__.py')
    )
    try:
        pypi_version = pypi_versions(package_name)[-1]
    except:
        pypi_version = '0.0'
    proposed_version = current_version.base_version
    if parse_version(proposed_version) <= parse_version(pypi_version):
        proposed_version = ''
    new_version = '0.0'
    while parse_version(new_version) <= parse_version(pypi_version):
        new_version = click.prompt(
            "What version > %s would you like to release?" % pypi_version,
            default=proposed_version,
        )
    click.confirm("Confirm version %s?" % new_version, abort=True)
    return str(new_version)


def set_version(filename, version):
    """Set the package version (in main __init__.py)"""
    shutil.copyfile(filename, filename + '.bak')
    click.echo("Modifying %s to set version %s" % (filename, version))
    with open(filename + '.bak') as in_fh, open(filename, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                line = line.split('=')[0].rstrip() + " = '" + version + "'\n"
            out_fh.write(line)
    if get_version(filename) == parse_version(version):
        os.remove(filename + ".bak")
    else:
        # roll back
        shutil.copyfile(filename + ".bak", filename)
        raise ReleaseError(
            "Failed to set version in %s (restored original)" % filename
        )


def edit_history(version):
    """Interactively edit HISTORY.rst"""
    click.echo(
        "Edit HISTORY.rst to add changelog and release date for %s" % version
    )
    edit('HISTORY.rst')
    click.confirm("Is HISTORY.rst up to date?", default=True, abort=True)


def check_dist():
    """Quietly make dist and check it. This is mainly to ensure that the README
    and HISTORY metadata are well-formed"""
    click.echo("Making and verifying dist and metadata...")
    try:
        run(['make', 'dist'], check=True, stdout=DEVNULL)
        run(['make', 'dist-check'], check=True)
        return True
    except CalledProcessError as exc_info:
        click.echo("ERROR: %s" % str(exc_info))
        return False


def check_docs():
    """Verify the documentation (interactively)"""
    click.echo("Making the documentation....")
    run(['make', 'docs'], check=True, stdout=DEVNULL)
    click.echo(
        "Check documentation in file://"
        + os.getcwd()
        + "/docs/_build/html/index.html"
    )
    click.confirm(
        "Does the documentation look correct?", default=True, abort=True
    )


def make_release_commit(version):
    """Commit 'Bump version to xxx and update HISTORY'"""
    click.confirm("Make release commit?", default=True, abort=True)
    run(
        [
            'git',
            'commit',
            '-a',
            '-m',
            "Bump version to %s and update HISTORY" % version,
        ],
        check=True,
    )


def make_upload(test=True):
    """Upload to PyPI or test.pypi"""
    if test:
        cmd = ['make', 'test-upload']
        url = 'https://test.pypi.org'
    else:
        url = 'https://pypi.org'
        cmd = ['make', 'upload']
    click.confirm(
        "Ready to upload release to %s?" % url, default=True, abort=True
    )
    success = False
    while not success:
        try:
            run(cmd, check=True)
        except CalledProcessError as exc_info:
            click.confirm(
                "Failed to upload: %s. Try again?" % str(exc_info),
                default=True,
                abort=(not test),
            )
            success = False
        else:
            success = True
            click.confirm(
                "Please check release on %s. Continue?" % url,
                default=True,
                abort=True,
            )


def push_release_commit():
    """Push local commits to origin"""
    click.confirm("Push release commit to origin?", default=True, abort=True)
    run(['git', 'push', 'origin', 'master'], check=True)
    click.confirm(
        "Please check Continuous Integration success. Continue?",
        default=True,
        abort=True,
    )


def make_and_push_tag(version):
    """Tag the current commit and push that tag to origin"""
    click.confirm(
        "Push tag '%s' to origin?" % version, default=True, abort=True
    )
    run(['git', 'tag', "v%s" % version], check=True)
    run(['git', 'push', '--tags', 'origin'], check=True)


def get_next_dev_version(released_version):
    """Ask for the post-release version number"""
    released_version = parse_version(str(released_version))
    if released_version.is_prerelease:
        return str(released_version.base_version) + "-dev"
    else:
        return str(released_version.base_version) + "+dev"


def make_next_dev_version_commit(version):
    """Commit 'Bump version to xxx'"""
    click.confirm(
        "Make commit for bumping to %s?" % version, default=True, abort=True
    )
    run(
        ['git', 'commit', '-a', '-m', "Bump version to %s" % version],
        check=True,
    )


@click.command(help=__doc__)
@click.help_option('--help', '-h')
@click.argument('package_name')
def main(package_name):
    try:
        make_release(package_name)
    except Exception as exc_info:
        click.echo(str(exc_info))
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
