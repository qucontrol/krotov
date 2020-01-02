#!/usr/bin/env docs
"""Clean development artifacts.

Usage: clean.py <selection>

with <selection> being one of the following:

    tests: remove coverage and pytest caches

    build: remove build, dist, egg, __pycache__ folders

    docs: remove built documentation

    venv: remove tox environments (.tox folder)

    all:  all of the above
"""
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).parent.parent

DOCSDIR = ROOT / 'docs'

DOCSBUILDDIR = DOCSDIR / '_build'

# fmt: off
FILES_TO_DELETE = {
    'tests': [
        (ROOT, '.coverage*'),
        ROOT / 'htmlcov',
        ROOT / '.pytest_cache'
    ],
    'build': [
        ROOT / 'build',
        ROOT / 'dist',
        (ROOT / 'src', '*.egg-info'),
        (ROOT, '[!.]*/**/__pycache__'),
        (ROOT, '**/.DS_Store'),
    ],
    'docs': [
        DOCSBUILDDIR,
        (DOCSDIR / 'API', '*.rst'),
    ],
    'venv': [
        ROOT / '.tox',
        ROOT / '.venv'
    ],
}
# fmt: on


def main(selection):
    """Main function"""
    if selection == 'all':
        for key in FILES_TO_DELETE:
            main(key)
    else:
        for entry in FILES_TO_DELETE[selection]:
            if isinstance(entry, tuple):
                path, pattern = entry
                for file_or_folder in path.glob(pattern):
                    if file_or_folder.is_file():
                        print("Remove file %s" % file_or_folder)
                        file_or_folder.unlink()
                    elif file_or_folder.is_dir():
                        print("Remove folder %s" % file_or_folder)
                        shutil.rmtree(file_or_folder)
            elif entry.is_file():
                print("Remove file %s" % entry)
                entry.unlink()
            elif entry.is_dir():
                print("Remove folder %s" % entry)
                shutil.rmtree(entry)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[-1])
