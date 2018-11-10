"""High-level tests for `krotov` package."""

from pkg_resources import parse_version

import krotov


def test_valid_version():
    """Check that the package defines a valid __version__"""
    assert parse_version(krotov.__version__) >= parse_version("0.0.1")
