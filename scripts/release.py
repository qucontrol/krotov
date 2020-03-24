#!/usr/bin/env python
"""Automation script for making a release. Must be run from the root for the
repository"""
# Note: Version scheme according to https://www.python.org/dev/peps/pep-0440
import json
import os
import re
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request
from os.path import join
from subprocess import DEVNULL, CalledProcessError, run
from textwrap import dedent

import click
import git
import pytest
from pkg_resources import parse_version


RX_VERSION = re.compile(
    r'^(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+)'
    r'(?P<prepost>\.post\d+|-(dev|a|b|rc)\d+)?'
    r'(?P<devsuffix>[+-]dev)?$'
)

RX_NBVIEWER_URL = re.compile(
    r'(?P<baseurl>https?://nbviewer\.jupyter\.org/[\w/-]+/blob/)'
    r'(?P<branch>[\w.-]+)(?P<path>[\w%/.-]+)'
)
RX_BINDER_URL = re.compile(
    r'(?P<baseurl>https?://mybinder\.org/v\d/[\w/-]+/)'
    r'(?P<branch>[\w.-]+)(?P<path>\?filepath=[\w%/.-]+)'
)

RX_STABLE_DOCS_LINK = re.compile(
    r'''
    (?P<baseurl>
      https://(?P<orgname>[\w-]+)\.github\.io/(?P<projname>[\w-]+)/)
    (?P<version>v\d[\w.-]+)
    (?P<path>/[\w%/.#-]+)
    ''',
    re.X,
)


def make_release(package_name, fast_test=False):
    """Interactively create and publish a new release for the package.

    If `fast_test` is passed as True, do a fast test-run of the release script.
    This skips testing, and verification, does not push any
    commits, and make no releases on PyPI. It does make local commits, which
    should be reset after the script completes.
    """
    click.confirm("Do you want to make a release?", abort=True)
    check_git_clean()
    new_version = ask_for_release_version(package_name)
    set_version(join('.', 'src', package_name, '__init__.py'), new_version)
    edit_history(new_version)
    while not check_dist():
        click.confirm(
            "Fix errors manually! Continue?", default=True, abort=True
        )
    files_with_binder_links = [
        'README.rst',
        join('docs', '09_examples.rst'),
        join('docs', 'index.rst'),
    ]
    for filename in files_with_binder_links:
        set_binder_branch(filename, "v" + str(new_version))
    set_binder_package_version(version=new_version)
    set_stable_docs_links(new_version, 'README.rst', 'docs/_README.patch')
    make_release_commit(new_version)
    make_notebooks(fast_test=fast_test)
    make_notebooks_commit(new_version)
    squash_commits(n=2, commit_msg="Release %s" % new_version)
    if not fast_test:
        make_artifacts()
        run_tests()
        make_upload(test=True)
        push_release_commits()
        make_upload(test=False)
        make_and_push_tag(new_version)
    upload_artifacts(new_version)
    # release is finished; go back to a dev state
    next_dev_version = new_version + '+dev'
    set_version(
        join('.', 'src', package_name, '__init__.py'), next_dev_version
    )
    files_with_released_binder_links = [
        'README.rst',
        join('docs', 'index.rst'),
    ]
    for filename in files_with_binder_links:
        if filename in files_with_released_binder_links:
            # README and index always link to the latest released version
            continue
        set_binder_branch(filename, 'master')
    set_binder_package_version(branch='master')
    make_next_dev_version_commit(next_dev_version)


###############################################################################


class ReleaseError(ValueError):
    pass


def get_package_name():
    """Find and return the package name from src"""
    for name in os.listdir('src'):
        if 'egg-info' in name:
            continue
        if os.path.isdir(os.path.join('src', name)):
            return name
    raise ReleaseError("Cannot find package name")


def get_pypi_versions(package_name):
    """Return list of versions for the given package on PyPI"""
    url = "https://pypi.python.org/pypi/%s/json" % (package_name,)
    data = json.load(urllib.request.urlopen(urllib.request.Request(url)))
    versions = list(data["releases"].keys())
    versions.sort(key=parse_version)
    return versions


def get_local_versions():
    """Return list of versions based on local tags

    For every version, there must be a tag "v<version>"
    """
    repo = git.Repo(os.getcwd())
    return [tag.name[1:] for tag in repo.tags if tag.name.startswith('v')]


def get_version(filename):
    """Extract the package version, as a str"""
    with open(filename) as in_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                return line.split('=')[1].strip()[1:-1]
    raise ReleaseError("Cannot extract version from %s" % filename)


def edit(filename):
    """Open filename in EDITOR."""
    editor = os.getenv('EDITOR', 'vi')
    if click.confirm("Open %s in %s?" % (filename, editor), default=True):
        run([editor, filename])


def check_git_clean():
    """Ensure that a given git.Repo is clean."""
    repo = git.Repo(os.getcwd())
    if repo.is_dirty():
        run(['git', 'status'])
        raise ReleaseError("Repository must be in a clean state")
    if repo.untracked_files:
        click.echo("WARNING: there are untracked files:")
        for filename in repo.untracked_files:
            click.echo("\t%s" % filename)
        click.confirm("Continue?", default=False, abort=True)


def run_tests():
    """Run 'make test'."""
    success = False
    while not success:
        try:
            run(['make', 'test'], check=True)
        except CalledProcessError as exc_info:
            print("Failed tests: %s\n" % exc_info)
            print("Fix the tests and ammend the release commit.")
            print("Then continue.\n")
            click.confirm("Continue?", default=True, abort=True)
            if not click.confirm("Retry?", default=True):
                break
        else:
            success = True


def split_version(version, base=True):
    """Split `version` into a tuple

    If `base` is True, only return (<major>, <minor>, <patch>) as a tuple of
    ints, stripping out pre/post/dev release tags. Otherwise, the are included
    as a possible fourth and fifth element in the tuple (as strings)
    """
    version = str(version)
    if not RX_VERSION.match(version):
        raise ValueError("Invalid version: %s" % version)
    if base:
        return tuple(
            [
                int(v)
                for v in str(parse_version(version).base_version).split(".")
            ]
        )
    else:
        m = RX_VERSION.match(version)
        if m:
            res = [
                int(m.group('major')),
                int(m.group('minor')),
                int(m.group('patch')),
            ]
            if m.group('prepost') is not None:
                res.append(m.group('prepost'))
            if m.group('devsuffix') is not None:
                res.append(m.group('devsuffix'))
            return tuple(res)
        else:
            raise ValueError("Invalid version string: %s" % version)


def list_versions(package_name):
    """List previously released versions

    This prints each released version on a new line, and returns the list of
    all released versions (based on PyPI and local tags)
    """
    try:
        pypi_versions = get_pypi_versions(package_name)
    except OSError:
        click.echo("PyPI versions no available")
        pypi_versions = []
    local_versions = get_local_versions()
    versions = sorted(
        set(pypi_versions).union(local_versions), key=parse_version
    )
    for version in versions:
        if version in pypi_versions and version in local_versions:
            status = 'PyPI/local'
        elif version in pypi_versions:
            status = 'PyPI only!'
        elif version in local_versions:
            status = 'local only!'
        click.echo("%-20s %s" % (version, status))
    return versions


def version_ok(version, dev_version, released_versions=None):
    """Check that `version` is a valid version for an upcoming release

    The `version` must be newer than the `dev_version` (from __version__, which
    should end in '-dev' or '+dev')
    """
    if released_versions is None:
        released_versions = []
    m = RX_VERSION.match(version)
    if m:
        if m.group('devsuffix') is not None:
            click.echo("Version %s contains a development suffix" % version)
            return False
        if version in released_versions:
            click.echo("Version %s is already released" % version)
            return False
        if parse_version(version) > parse_version(dev_version):
            return True
        else:
            click.echo("Version %s not newer than %s" % (version, dev_version))
            return False
    else:
        click.echo("Invalid version: %s" % version)
        return False


def propose_next_version(dev_version):
    """Return the most likely release version based on the current
    __version__"""
    dev_version = str(dev_version)
    if parse_version(dev_version).is_prerelease:
        return parse_version(dev_version).base_version
    else:
        base_version = parse_version(dev_version).base_version
        v = split_version(base_version)
        return "%d.%d.%d" % (v[0], v[1], v[2] + 1)


def ask_for_release_version(package_name):
    """Ask for the version number of the release.

    The version number is checked to be a valid next release
    """
    dev_version = get_version(join('.', 'src', package_name, '__init__.py'))
    proposed_version = propose_next_version(dev_version)
    released_versions = list_versions(package_name)
    new_version = click.prompt(
        "What version would you like to release?", default=proposed_version
    )
    while not version_ok(new_version, dev_version, released_versions):
        new_version = click.prompt(
            "What version would you like to release?", default=proposed_version
        )
    click.confirm("Confirm version %s?" % new_version, abort=True)
    return str(new_version)


def set_version(filename, version):
    """Set the package version (in main __init__.py)"""
    shutil.copyfile(filename, filename + '.bak')
    click.echo("Modifying %s to set version %s" % (filename, version))
    with open(filename + '.bak') as in_fh, open(filename, 'w') as out_fh:
        found_version_line = False
        for line in in_fh:
            if line.startswith('__version__'):
                found_version_line = True
                line = line.split('=')[0].rstrip() + " = '" + version + "'\n"
            out_fh.write(line)
    if get_version(filename) == version:
        os.remove(filename + ".bak")
    else:
        # roll back
        shutil.copyfile(filename + ".bak", filename)
        msg = "Failed to set version in %s (restored original)" % filename
        if not found_version_line:
            msg += ". Does not contain a line starting with '__version__'."
        raise ReleaseError(msg)


def set_binder_branch(filename, branch='master'):
    """Search for nbviewer/binder URLs and replace their branch name"""
    shutil.copyfile(filename, filename + '.bak')
    click.echo(
        "Modifying %s to set branch %s in binder/nbviewer URLs"
        % (filename, branch)
    )
    with open(filename + '.bak') as in_fh, open(filename, 'w') as out_fh:
        for line in in_fh:
            line = RX_NBVIEWER_URL.sub(r'\g<1>%s\g<3>' % branch, line)
            line = RX_BINDER_URL.sub(r'\g<1>%s\g<3>' % branch, line)
            out_fh.write(line)
    os.remove(filename + ".bak")


def set_binder_package_version(version=None, branch=None):
    """Set the version of the krotov package in binder/environment.yml.

    Either a `version` must be given, in which case Binder dependes on
    ``krotov==<version>``, or `branch`, in which case Binder depens on
    ``git+https://github.com/qucontrol/krotov.git@<branch>#egg=krotov``.
    """
    filename = 'binder/environment.yml'
    shutil.copyfile(filename, filename + '.bak')
    # fmt: off
    with open(filename + '.bak') as in_fh, open(filename, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('    - ') and 'krotov' in line:
                if version is not None:
                    click.echo("Modifying %s to set krotov version %s" % (filename, version))
                    line = '    - krotov==%s' % version
                elif branch is not None:
                    click.echo("Modifying %s to set krotov branch %s" % (filename, branch))
                    line = '    - git+https://github.com/qucontrol/krotov.git@%s#egg=krotov' % branch
                else:
                    raise ValueError("You must give either a version or a branch")
            out_fh.write(line)
    # fmt: on
    os.remove(filename + ".bak")


def set_stable_docs_links(version, *files):
    """Replace documentation-links in the given files to point to stable.

    Replace all strings that match RX_STABLE_DOCS_LINK with the correct
    projname with a URL that points to `version`.
    """
    projname = get_package_name()
    for filename in files:
        shutil.copyfile(filename, filename + '.bak')
        with open(filename + '.bak') as in_fh, open(filename, 'w') as out_fh:
            for line in in_fh:
                replacements = {}
                for m in RX_STABLE_DOCS_LINK.finditer(line):
                    if m.group('projname') == projname:
                        replacements[m.group(0)] = (
                            m.group('baseurl')
                            + 'v%s' % version
                            + m.group('path')
                        )
                for url, replacement in replacements.items():
                    click.echo(
                        "In %s, replace '%s' → '%s'"
                        % (filename, url, replacement)
                    )
                    line = line.replace(url, replacement)
                out_fh.write(line)
        os.remove(filename + ".bak")


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


def make_notebooks(fast_test=False):
    """Re-run all notebooks in the documentation.

    If `fast_test` is given as True, only the first notebook is run.
    """
    click.echo("Re-creating notebooks...")
    target = 'notebooks'
    if fast_test:
        # only make the first notebook
        target = 'docs/notebooks/01_example_simple_state_to_state.ipynb.log'
    try:
        run(['make', target], check=True)
        return True
    except CalledProcessError as exc_info:
        click.echo("ERROR: %s" % str(exc_info))
        return False


def make_artifacts():
    """Verify the documentation (interactively)."""
    click.echo("Making the documentation....")
    run(['make', 'docs-artifacts'], check=True, stdout=DEVNULL)
    click.echo(
        "Check documentation in file://"
        + os.getcwd()
        + "/docs/_build/html/index.html"
        + " and the artifacts in "
        + os.getcwd()
        + "/docs/_build/artifacts"
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


def make_notebooks_commit(version):
    """Commit 'Reevaluate notebooks for version xxx'"""
    click.confirm("Make commit for notebooks?", default=True, abort=True)
    run(
        [
            'git',
            'commit',
            '-a',
            '-m',
            "Reevaluate notebooks for version %s" % version,
        ],
        check=True,
    )


def squash_commits(n, commit_msg=None):
    """Squash release commits"""
    msg = (
        "Ready to squash %d release commits (change 'pick' to 's' for all "
        "but first commit)?" % n
    )
    if commit_msg:
        msg = (
            "Ready to squash %d release commits (change 'pick' to 's' for all "
            "but first commit, use '%s' as the commmit message')?"
            % (n, commit_msg)
        )
    click.confirm(
        msg, default=True, abort=True,
    )
    run(['git', 'rebase', '--interactive', 'HEAD~%d' % n], check=True)


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


def push_release_commits():
    """Push local commits to origin."""
    click.confirm("Push release commit to origin?", default=True, abort=True)
    run(['git', 'push', 'origin', 'master'], check=True)
    click.confirm(
        "Please check Continuous Integration success. Continue?",
        default=True,
        abort=True,
    )


def make_and_push_tag(version):
    """Tag the current commit and push that tag to origin."""
    click.confirm(
        "Push tag '%s' to origin?" % version, default=True, abort=True
    )
    run(['git', 'tag', "-s", "v%s" % version], check=True)
    run(['git', 'push', '--tags', 'origin'], check=True)


def upload_artifacts(version):
    """Prompt for manually uploading documentation artifacts."""
    # fmt: off
    click.echo(dedent(r'''
    You must now manually create a release on Github and upload the
    documentation artifacs from docs/_build/artifacts.

    * Go to https://github.com/qucontrol/krotov/releases
    * Find and click on the tag 'v{version}'
    * Click on "Edit tag"
    * Use "Release {version}" as the Release title
    * Add release notes in markdown format
    * Attach the documentation artifacs as binary files
    * Click "Publish release"

    '''.format(version=version)))
    # fmt: on
    click.confirm("Release published on Github?", default=True)


def make_next_dev_version_commit(version):
    """Commit 'Bump version to xxx'."""
    click.confirm(
        "Make commit for bumping to %s?" % version, default=True, abort=True
    )
    run(
        ['git', 'commit', '-a', '-m', "Bump version to %s" % version],
        check=True,
    )


###############################################################################


# run tests with `pytest -s scripts/release.py`


def test_list_versions():
    print("")
    versions = list_versions(get_package_name())
    print(versions)
    assert isinstance(versions, list)


def test_split_version():
    # fmt: off
    assert split_version('0.1.0') == (0, 1, 0)
    assert split_version('0.1.0', base=False) == (0, 1, 0)
    assert split_version('0.1.0-dev1', base=True) == (0, 1, 0)
    assert split_version('0.1.0-dev1', base=False) == (0, 1, 0, '-dev1')
    assert split_version('0.1.0.post1', base=True) == (0, 1, 0)
    assert split_version('0.1.0.post1', base=False) == (0, 1, 0, '.post1')
    assert split_version('0.1.0-rc1', base=True) == (0, 1, 0)
    assert split_version('0.1.0-rc1', base=False) == (0, 1, 0, '-rc1')
    assert split_version('0.1.0-rc1-dev', base=True) == (0, 1, 0)
    assert split_version('0.1.0-rc1-dev', base=False) == (0, 1, 0, '-rc1', '-dev')
    assert split_version('0.1.0-rc1+dev', base=True) == (0, 1, 0)
    assert split_version('0.1.0-rc1+dev', base=False) == (0, 1, 0, '-rc1', '+dev')
    assert split_version('0.1.0-dev', base=True) == (0, 1, 0)
    assert split_version('0.1.0-dev', base=False) == (0, 1, 0, '-dev')
    assert split_version('0.1.0+dev', base=True) == (0, 1, 0)
    assert split_version('0.1.0+dev', base=False) == (0, 1, 0, '+dev')
    with pytest.raises(ValueError):
        split_version('0.1.0.rc1')
    with pytest.raises(ValueError):
        split_version('0.1.0rc1')
    with pytest.raises(ValueError):
        split_version('0.1.0.1')
    with pytest.raises(ValueError):
        split_version('0.1')
    with pytest.raises(ValueError):
        split_version('0.1.0+dev1')
    # fmt: on


def test_version_ok():
    assert version_ok('0.1.0', '0.1.0-dev')
    assert version_ok('0.1.0-a1', '0.1.0-dev')
    assert version_ok('0.1.0-b1', '0.1.0-dev')
    assert version_ok('0.1.0-rc1', '0.1.0-dev')
    assert version_ok('0.2.0', '0.1.0+dev')
    assert version_ok('0.2.0-a1', '0.1.0+dev')
    assert version_ok('0.2.0-b1', '0.1.0+dev')
    assert version_ok('0.2.0-rc1', '0.1.0+dev')
    assert version_ok('0.2.0-dev1', '0.1.0+dev')
    assert version_ok('0.1.0.post1', '0.1.0+dev')
    assert version_ok('0.1.0.post1', '0.1.0')
    assert version_ok('0.2.0', '0.1.0')
    assert version_ok('0.2.0', '0.1.0+dev', ['0.1.0', '0.1.0.post1', '0.1.1'])
    print("")
    assert not version_ok('0.0.1-dev', '0.1.0-dev')
    assert not version_ok('0.1.0', '0.1.0')
    assert not version_ok('0.1.0', '0.1.0+dev')
    assert not version_ok('0.1.0+dev', '0.1.0')
    assert not version_ok('0.2.0-dev', '0.1.0+dev')
    assert not version_ok('0.1.0.1', '0.1.0-dev')
    assert not version_ok('0.1.0a1', '0.1.0-dev')
    assert not version_ok('0.1.0b1', '0.1.0-dev')
    assert not version_ok('0.1.0rc1', '0.1.0-dev')
    assert not version_ok('0.1.0dev1', '0.1.0-dev')
    assert not version_ok('0.1.0-post1', '0.1.0+dev')
    assert not version_ok('0.2.0', '0.1.0+dev', ['0.1.0', '0.2.0'])


def test_propose_next_version():
    assert propose_next_version('0.1.0') == '0.1.1'
    assert propose_next_version('0.1.0-dev') == '0.1.0'
    assert propose_next_version('0.1.0-rc1') == '0.1.0'
    assert propose_next_version('0.1.0-rc1+dev') == '0.1.0'
    assert propose_next_version('0.1.0+dev') == '0.1.1'
    assert propose_next_version('0.1.0.post1') == '0.1.1'
    assert propose_next_version('0.1.0.post1+dev') == '0.1.1'


def test_nbviewer_binder_regexes():
    url_nbviewer = (
        r'http://nbviewer.jupyter.org/github/qucontrol/krotov/'
        r'blob/master/docs/notebooks/'
        r'03_example_lambda_system_rwa_non_hermitian.ipynb'
    )
    url_binder = (
        r'https://mybinder.org/v2/gh/qucontrol/krotov/master'
        r'?filepath=docs/notebooks/'
        r'03_example_lambda_system_rwa_non_hermitian.ipynb'
    )

    match = RX_NBVIEWER_URL.match(url_nbviewer)
    assert match
    assert (
        match.group('baseurl')
        == 'http://nbviewer.jupyter.org/github/qucontrol/krotov/blob/'
    )
    assert match.group('branch') == 'master'
    assert (
        match.group('path')
        == '/docs/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb'
    )
    new_url_nbviewer = RX_NBVIEWER_URL.sub(
        r'\g<1>%s\g<3>' % 'v0.1.0', url_nbviewer
    )
    assert new_url_nbviewer == (
        r'http://nbviewer.jupyter.org/github/qucontrol/krotov/blob/v0.1.0/'
        r'docs/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb'
    )
    assert RX_NBVIEWER_URL.match(new_url_nbviewer).group('branch') == 'v0.1.0'

    match = RX_BINDER_URL.match(url_binder)
    assert match
    assert (
        match.group('baseurl')
        == 'https://mybinder.org/v2/gh/qucontrol/krotov/'
    )
    assert match.group('branch') == 'master'
    assert match.group('path') == (
        r'?filepath=docs/notebooks/'
        r'03_example_lambda_system_rwa_non_hermitian.ipynb'
    )
    new_url_binder = RX_BINDER_URL.sub(r'\g<1>%s\g<3>' % 'v0.1.0', url_binder)
    assert new_url_binder == (
        r'https://mybinder.org/v2/gh/qucontrol/krotov/v0.1.0'
        r'?filepath=docs/notebooks/'
        r'03_example_lambda_system_rwa_non_hermitian.ipynb'
    )
    assert RX_BINDER_URL.match(new_url_binder).group('branch') == 'v0.1.0'


###############################################################################


@click.command(help=__doc__)
@click.help_option('--help', '-h')
@click.option(
    '--test', is_flag=True, help="Do a fast test-run of the release script"
)
def main(test):
    try:
        make_release(get_package_name(), fast_test=bool(test))
    except Exception as exc_info:
        click.echo(str(exc_info))
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
