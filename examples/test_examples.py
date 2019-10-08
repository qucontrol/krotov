"""Verify the example scripts."""

import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pytest

PYTHON = sys.executable


def get_examples():
    """Return a list of all .py/.out pairs in the current folder."""
    res = []
    for expected_output_file in Path(__file__).parent.glob('*.out'):
        script_file = expected_output_file.with_suffix('.py')
        if script_file.is_file:
            res.append((script_file, expected_output_file))
    return res


# fmt: off
def sanitize(output):
    """Sanitize the script output, before comparing."""
    replacements = [
        (r'\d{4}-\d{2}-\d{2}', 'DATE-STAMP'),
        (r'\d{1,2}:\d{2}:\d{2}', 'TIME-STAMP'),
        (r'(?<=[\s\d]{5}[\s\d.e+-]{9}[\s\d.e+-]{12}([\s\d.e+-]{11}){3})[\s\d]{6}\n', "\n"),
        (r'(?<=[\s\d]{5}[\s\d.e+-]{9}[\s\d.e+-]{12}([\s\d.e+-]{11}){1}(        n/a){2})[\s\d]{6}\n', "\n"),
    ]
    for (rx, replacement) in replacements:
        output = re.sub(rx, replacement, output)
    return output
# fmt: on


@contextmanager
def protect(*files):
    """Protect the given files from modification."""
    with tempfile.TemporaryDirectory() as tmpdir:
        for file in files:
            tmpfile = Path(tmpdir) / file.name
            shutil.copyfile(file, tmpfile)
        yield
        for file in files:
            tmpfile = Path(tmpdir) / file.name
            shutil.copyfile(tmpfile, file)


@pytest.mark.parametrize("script_file,expected_output_file", get_examples())
def test_example(script_file, expected_output_file):
    """Test that running the `script_file` generates the expected output."""
    print(f"Checking {script_file}")
    root = Path(__file__).parent
    with protect(root / 'tls_state_to_state.dump'):
        completed = subprocess.run(
            [PYTHON, script_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding='utf8',
            check=True,
        )
    output = sanitize(completed.stdout)
    expected_output = sanitize(
        Path(expected_output_file).read_text(encoding='utf8')
    )
    assert output == expected_output
