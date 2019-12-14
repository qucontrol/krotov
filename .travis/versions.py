"""Utilities for generating versions.json."""
import json
import subprocess
from pathlib import Path

from pkg_resources import parse_version
from pkg_resources.extern.packaging.version import LegacyVersion


def get_stable(versions):
    """Select the latest stable release from `versions`.

    If there is no stable version, return None.
    """
    try:
        return sorted(
            [v for v in versions if not parse_version(v).is_prerelease],
            key=parse_version,
        )[-1]
    except IndexError:
        return None


def _clear_broken_symlinks():
    """Clear all broken symlinks in the current directory."""
    for name in Path().iterdir():
        if name.is_symlink():
            try:
                name.resolve(strict=True)
            except FileNotFoundError:
                print("Remove broken symlink %s" % name)
                subprocess.run(['git', 'rm', name], check=True)


def write_versions_json(versions_data, outfile):
    """Write the versions data to a json file.

    This json file will be processed by the javascript that generates the
    version-selector.
    """
    with open(outfile, 'w') as out_fh:
        json.dump(versions_data, out_fh)
    print(json.dumps(versions_data, indent=2))
    subprocess.run(['git', 'add', outfile], check=True)


def _is_unreleased(folder):
    """Default `is_unreleased` function.

    The following are considered "unreleased":

    * Anything that doesn't look like a proper version according to PEP440.
      Proper versions are e.g. "1.0.0", "v1.0.0", "1.0.0.post1". Specifically,
      any branch names like "master", "develop" are considered unreleased.
    * Anything that PEP440 considers a pre-release, e.g. "1.0.0-dev", "1.0-rc1"
    * Anything that includes a "local version identifier" according to PEP440,
      e.g. "1.0.0+dev"
    """
    version = parse_version(folder)
    if isinstance(version, LegacyVersion):
        return True
    if version.is_prerelease:
        return True
    if version.local is not None:
        return True
    return False


def _find_latest_release(folders):
    try:
        return sorted(folders, key=parse_version)[-1]
    except IndexError:
        return None


def _find_downloads(folder):
    return []


def get_versions_data(
    hidden=None,
    is_unreleased=None,
    find_latest_release=None,
    sort_key=None,
    labels=None,
    suffix_latest_release=' (latest release)',
    suffix_unreleased=' (dev)',
    find_downloads=None,
):
    """Get the versions data, to be serialized to json."""
    if hidden is None:
        hidden = []
    if is_unreleased is None:
        is_unreleased = _is_unreleased
    if find_latest_release is None:
        find_latest_release = _find_latest_release
    if find_downloads is None:
        find_downloads = _find_downloads
    if sort_key is None:
        sort_key = parse_version
    if labels is None:
        labels = {}
    folders = sorted(
        [
            str(f)
            for f in Path().iterdir()
            if (
                f.is_dir()
                and not f.is_symlink()
                and not str(f).startswith('.')
                and not str(f).startswith('_')
            )
        ],
        key=sort_key,
    )
    labels = {folder: labels.get(folder, str(folder)) for folder in folders}
    versions = []
    unreleased = []
    for folder in folders:
        if folder not in hidden:
            versions.append(folder)
        if is_unreleased(folder):
            unreleased.append(folder)
            labels[folder] += suffix_unreleased
    latest_release = find_latest_release(
        [f for f in versions if f not in unreleased]
    )
    outdated = []
    if latest_release is not None:
        labels[latest_release] += suffix_latest_release
        outdated = [
            folder
            for folder in versions
            if (folder != latest_release and folder not in unreleased)
        ]
    versions_data = {
        # list of *all* folders
        'folders': folders,
        #
        # folder => labels for every folder in "Versions"
        'labels': labels,
        #
        # list folders that appear in "Versions"
        'versions': versions,
        #
        # list of folders that do not appear in "Versions"
        'hidden': hidden,
        #
        # list of folders that should warn & point to latest release
        'outdated': outdated,
        #
        # list of dev-folders that should warn & point to latest release
        'unreleased': unreleased,
        #
        # the latest stable release folder
        'latest_release': latest_release,
        #
        # folder => list of (label, file)
        'downloads': {folder: find_downloads(folder) for folder in folders},
    }

    return versions_data
