#!/usr/bin/env python
import os
import subprocess
from pathlib import Path

from versions import get_versions_data, write_versions_json


INDEX_HTML = r'''<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="Refresh" content="0; url={default_branch}" />
  </head>
  <body>
    <p>Got to <a href="{default_branch}">default documentation</a>.</p>
  </body>
</html>
'''


def write_index_html(default_branch):
    """Write an index.html that redirects to the DEFAULT_BRANCH."""
    with open("index.html", "w") as out_fh:
        out_fh.write(INDEX_HTML.format(default_branch=default_branch))
    subprocess.run(['git', 'add', 'index.html'], check=True)


def find_downloads(folder):
    """Find files in the 'download' subfolder of the given `folder`."""
    downloads = []
    for filename in Path(folder).glob(r'download/*'):
        label = "".join(filename.suffixes).replace('.', '').lower()
        if len(label) > 0:
            downloads.append((label, str(filename)))
    return downloads


def main():
    """Main function."""
    print("Post-processing documentation on gh-pages")
    print("Gather versions info")
    versions_data = get_versions_data(find_downloads=find_downloads)
    latest_release = versions_data['latest_release']
    if latest_release is None:
        latest_release = 'master'
    print("Write index.html")
    write_index_html(latest_release)
    print("Write versions.json")
    write_versions_json(versions_data, outfile='versions.json')
    print("DONE post-processing")


if __name__ == "__main__":
    main()
