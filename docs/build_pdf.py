"""Script for building the documentation in PDF format.

This script is invoked by `make docs-pdf` (Unix systems only!)

It assumes that you have lualatex installed (https://www.tug.org/texlive/).
You must also have the DejaVu fonts installed on your system
(https://dejavu-fonts.github.io).
"""

import re
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).parent


def get_version(filename):
    """Extract the package version"""
    with open(filename, encoding='utf8') as in_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                return line.split('=')[1].strip()[1:-1]
    raise ValueError("Cannot extract version from %s" % filename)


def _patch_line(line):
    if line.startswith('\sphinxhref{') and line.endswith(".svg}}\n"):
        return None
    if line == r'\chapter{References}' + "\n":
        return None
    if line.endswith(r'\label{\detokenize{99_bibliography::doc}}' + "\n"):
        return None
    line = line.replace('{{krotovscheme}.svg}', '{{krotovscheme}.pdf}')
    line = line.replace(
        '{{oct_decision_tree}.svg}', '{{oct_decision_tree}.pdf}'
    )
    line = line.replace(
        r'\sphinxurl{krotov\_pseudocode.pdf}',
        r'\url{https://qucontrol.github.io/krotov/master/krotov_pseudocode.pdf}',
    )
    if line.startswith(r'\(\newcommand'):
        return None
    if line == r'\begin{split}\begin{equation}' + "\n":
        return None
    if line == r'\end{equation}\end{split}' + "\n":
        return None
    if line == r'\begin{split}\begin{split}' + "\n":
        return r'\begin{split}' + "\n"
    if line == r'\end{split}\end{split}' + "\n":
        return r'\end{split}' + "\n"
    if line == r'  \end{split}\end{split}' + "\n":
        return r'\end{split}' + "\n"
    if line.startswith(r'\author{'):
        line = r'\author{Michael Goerz \textit{et. al.}}' + "\n"
    if line.startswith(r'\section{'):
        # don't put section numbers in the HISTORY, in front of version numbers
        match = re.match(
            r'\\section\{(\d+\.\d+\.\d+\s+\([\d-]+\)|[\s(]*next version[)\s]*)\}',
            line.strip(),
        )
        if match:
            line = r'\section*{' + match.group(1) + "}\n"
    return line


def patch_krotov_tex_lines(texfile):
    """Fix errors line-by-line in the given texfile."""
    with tempfile.TemporaryDirectory() as tmpdir:
        orig = Path(tmpdir) / texfile.name
        shutil.copyfile(texfile, orig)
        with orig.open() as in_fh, texfile.open("w") as out_fh:
            for line in in_fh:
                line = _patch_line(line)
                if line is None:
                    continue
                out_fh.write(line)


def _multiline_str(*lines):
    return "\n".join(lines)


def patch_krotov_tex(texfile):
    """Fix errors in the given texfile, acting on the whole text."""
    tex = texfile.read_text(encoding='utf8')
    tex = tex.replace(
        r'\begin{equation*}' + "\n" + r'\begin{split}\begin{align}',
        r'\begin{align*}',
    )
    tex = tex.replace(
        r'\end{align}\end{split}' + "\n" + r'\end{equation*}', r'\end{align*}'
    )
    tex = tex.replace(
        _multiline_str(
            r'\chapter{Indices and tables}',
            r'\label{\detokenize{index:indices-and-tables}}\begin{itemize}',
            r'\item {} ',
            r'\DUrole{xref,std,std-ref}{genindex}',
            r'',
            r'\item {} ',
            r'\DUrole{xref,std,std-ref}{modindex}',
            r'',
            r'\end{itemize}',
        ),
        '',
    )
    texfile.write_text(tex, encoding='utf8')


def latex(
    texfile,
    executable='lualatex',
    texliveonfly=None,
    stdout=None,
    log_filename=None,
):
    """Run lualatex to compile the given texfile.

    Args:
        texfile: path-like object pointing to tex file to compile
        exectuable (str): one of pdflatex, lualatex
        texliveonfly (bool or None): whether to use texliveonfly. If not given,
            defaults to True when running on Travis, False otherwise
    """
    if texliveonfly is None:
        if 'TRAVIS' in os.environ:
            print("Running on Travis: using texliveonfly (%s)" % executable)
            texliveonfly = True
        else:
            print("Not running on Travis: using %s" % executable)
            texliveonfly = False
    cmd = [
        executable,
        '--interaction=nonstopmode',
        '--halt-on-error',
        texfile.name,
    ]
    if texliveonfly:
        cmd = [
            'texliveonfly',
            '--compiler=%s' % executable,
            texfile.name,
        ]
    print(" ".join(cmd))
    res = subprocess.run(
        cmd, cwd=texfile.parent, check=True, capture_output=True,
    )
    if log_filename is not None:
        with (texfile.parent / log_filename).open("wb") as out_fh:
            out_fh.write(res.stdout)


def main():
    """Main function."""
    texfile = ROOT / '_build/tex/krotov.tex'
    if not texfile.is_file():
        print("%s does not exist" % texfile)
        sys.exit(1)
    print("Patching %s..." % texfile)
    patch_krotov_tex_lines(texfile)
    patch_krotov_tex(texfile)
    print("Compiling %s..." % texfile)
    latex(texfile)
    latex(texfile)
    latex(texfile, texliveonfly=False, log_filename='krotov_final.log')  # DEBUG
    print("Done compiling %s" % texfile)
    sys.exit(0)


if __name__ == "__main__":
    main()
