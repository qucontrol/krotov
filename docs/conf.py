# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest import mock

import git
import pybtex.style.formatting
import pybtex.style.formatting.unsrt
import pybtex.style.template
from pybtex.plugin import register_plugin as pybtex_register_plugin
from sphinx.ext.autodoc import (
    ClassLevelDocumenter,
    InstanceAttributeDocumenter,
)

import krotov


DOCS = Path(__file__).parent
ROOT = DOCS / ".."

sys.path.insert(0, str((DOCS / "_extensions").resolve()))

exclude_patterns = [
    '_build',
]

# -- Generate API documentation ------------------------------------------------
def run_apidoc(app):
    """Generage API documentation"""
    import better_apidoc

    better_apidoc.APP = app
    better_apidoc.main(
        [
            "better-apidoc",
            "-t",
            str(DOCS / "_templates"),
            "--force",
            "--no-toc",
            "--separate",
            "-o",
            str(DOCS / "API"),
            str(DOCS / ".." / "src" / "krotov"),
        ]
    )


# -- Generate patched README documentation ------------------------------------
def generate_patched_readme(_):
    shutil.copyfile(DOCS / ".." / "README.rst", DOCS / "_README.rst")
    cmd = ["patch", str(DOCS / "_README.rst"), str(DOCS / "_README.patch")]
    subprocess.run(cmd, check=True)
    assert (DOCS / "_README.rst").is_file()


# -- General configuration -----------------------------------------------------

# Report broken links as warnings
nitpicky = True
nitpick_ignore = [("py:class", "callable")]

extensions = [
    "doctr_versions_menu",
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
    "dollarmath",  # must be loaded after sphinx.ext.autodoc
]
if os.getenv("SPELLCHECK"):
    extensions += ("sphinxcontrib.spelling",)
    spelling_show_suggestions = True
    spelling_lang = os.getenv("SPELLCHECK")
    spelling_word_list_filename = "spelling_wordlist.txt"
    spelling_ignore_pypi_package_names = True

intersphinx_mapping = {
    # Note: the inventory "python37.inv" was patched to include a reference
    # for :class:`multiprocessing.context.Process`. The patched files was
    # created with the help of https://sphobjinv.readthedocs.io/
    "python": ("https://docs.python.org/3.7", "python37.inv"),
    "sympy": ("https://docs.sympy.org/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "qutip": ("http://qutip.org/docs/latest/", None),
    "glom": ("https://glom.readthedocs.io/en/latest/", None),
    "weylchamber": ("https://weylchamber.readthedocs.io/en/latest/", None),
    "loky": ("https://loky.readthedocs.io/en/stable/", None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

source_suffix = ".rst"
master_doc = "index"
project = "Krotov"
year = str(datetime.datetime.now().year)
author = "Michael Goerz"
copyright = "{0}, {1}".format(year, "Michael Goerz et al.")
version = krotov.__version__
release = version
git_tag = "v%s" % version
if version.endswith("dev"):
    try:
        last_commit = str(git.Repo(ROOT).head.commit)[:7]
        release = "%s (%s)" % (version, last_commit)
        git_tag = str(git.Repo(ROOT).head.commit)
    except git.exc.InvalidGitRepositoryError:
        git_tag = "master"
numfig = True

html_extra_path = ["./pseudocode/krotov_pseudocode.pdf"]

pygments_style = "friendly"
extlinks = {
    "issue": ("https://github.com/qucontrol/krotov/issues/%s", "#"),
    "pr": ("https://github.com/agkoch/krotov/pull/%s", "PR #"),
}

# autodoc settings
autoclass_content = "both"
autodoc_member_order = "bysource"
autodoc_mock_imports = [
    "numpy",
    "numpy.linalg",
    "scipy",
    "scipy.sparse",
    "matplotlib",
    "matplotlib.pyplot",
]


html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_sidebars = {"**": ["searchbox.html", "globaltoc.html", "sourcelink.html"]}
html_short_title = "%s-%s" % (project, version)

# Mathjax settings
mathjax_path = (
    "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js"
)
mathjax_config = {
    "extensions": ["tex2jax.js"],
    "jax": ["input/TeX", "output/SVG"],
    "TeX": {
        "extensions": ["AMSmath.js", "AMSsymbols.js"],
        "Macros": {
            "tr": ["{\\operatorname{tr}}", 0],
            "diag": ["{\\operatorname{diag}}", 0],
            "abs": ["{\\operatorname{abs}}", 0],
            "pop": ["{\\operatorname{pop}}", 0],
            "ee": ["{\\text{e}}", 0],
            "ii": ["{\\text{i}}", 0],
            "aux": ["{\\text{aux}}", 0],
            "opt": ["{\\text{opt}}", 0],
            "tgt": ["{\\text{tgt}}", 0],
            "init": ["{\\text{init}}", 0],
            "lab": ["{\\text{lab}}", 0],
            "rwa": ["{\\text{rwa}}", 0],
            "bra": ["{\\langle#1\\vert}", 1],
            "ket": ["{\\vert#1\\rangle}", 1],
            "Bra": ["{\\left\\langle#1\\right\\vert}", 1],
            "Braket": [
                "{\\left\\langle #1\\vphantom{#2} \\mid #2\\vphantom{#1}\\right\\rangle}",
                2,
            ],
            "ketbra": ["{\\vert#1\\rangle\\!\\langle#2\\vert}", 2],
            "Ket": ["{\\left\\vert#1\\right\\rangle}", 1],
            "mat": ["{\\mathbf{#1}}", 1],
            "op": ["{\\hat{#1}}", 1],
            "Op": ["{\\hat{#1}}", 1],
            "dd": ["{\\,\\text{d}}", 0],
            "daggered": ["{^{\\dagger}}", 0],
            "transposed": ["{^{\\text{T}}}", 0],
            "Liouville": ["{\\mathcal{L}}", 0],
            "DynMap": ["{\\mathcal{E}}", 0],
            "identity": ["{\\mathbf{1}}", 0],
            "Norm": ["{\\left\\lVert#1\\right\\rVert}", 1],
            "norm": ["{\\lVert#1\\rVert}", 1],
            "Abs": ["{\\left\\vert#1\\right\\vert}", 1],
            "avg": ["{\\langle#1\\rangle}", 1],
            "Avg": ["{\\left\langle#1\\right\\rangle}", 1],
            "AbsSq": ["{\\left\\vert#1\\right\\vert^2}", 1],
            "Re": ["{\\operatorname{Re}}", 0],
            "Im": ["{\\operatorname{Im}}", 0],
            "Real": ["{\\mathbb{R}}", 0],
            "Complex": ["{\\mathbb{C}}", 0],
            "Integer": ["{\\mathbb{N}}", 0],
        },
    },
}

# LaTeX settings
latex_engine = "lualatex"
latex_fontpkg = r"""
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
"""
latex_preamble = r"""
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
\usepackage{emptypage}
\usepackage{braket}
\newcommand{\tr}[0]{\operatorname{tr}}
\newcommand{\diag}[0]{\operatorname{diag}}
\newcommand{\abs}[0]{\operatorname{abs}}
\newcommand{\pop}[0]{\operatorname{pop}}
\newcommand{\aux}[0]{\text{aux}}
\newcommand{\opt}[0]{\text{opt}}
\newcommand{\tgt}[0]{\text{tgt}}
\newcommand{\init}[0]{\text{init}}
\newcommand{\lab}[0]{\text{lab}}
\newcommand{\rwa}[0]{\text{rwa}}
\renewcommand{\Braket}[2]{\left\langle{}#1\vphantom{#2}\mid{}#2\vphantom{#1}\right\rangle}
\newcommand{\ketbra}[2]{\vert#1\rangle\!\langle#2\vert}
\newcommand{\op}[1]{\hat{#1}}
\newcommand{\Op}[1]{\hat{#1}}
\newcommand{\dd}[0]{\,\text{d}}
\newcommand{\Liouville}[0]{\mathcal{L}}
\newcommand{\DynMap}[0]{\mathcal{E}}
\newcommand{\identity}[0]{\mathbf{1}}
\newcommand{\norm}[1]{\lVert#1\rVert}
\newcommand{\Norm}[1]{\left\lVert#1\right\rVert}
\newcommand{\Abs}[1]{\left\vert#1\right\vert}
\newcommand{\avg}[1]{\langle#1\rangle}
\newcommand{\Avg}[1]{\left\langle#1\right\rangle}
\newcommand{\AbsSq}[1]{\left\vert#1\right\vert^2}
\renewcommand{\Re}[0]{\operatorname{Re}}
\renewcommand{\Im}[0]{\operatorname{Im}}
"""
latex_printindex = r"\footnotesize\raggedright\printindex"
latex_fncychap = r"\usepackage[Bjornstrup]{fncychap}"
latex_elements = {
    "fontpkg": latex_fontpkg,
    "preamble": latex_preamble,
    "fncychap": latex_fncychap,
    "printindex": latex_printindex,
    "babel": "",
}
latex_show_urls = "no"


# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True


# Sphinxcontrib-bibtex
bibtex_bibfiles = ['refs.bib']
pybtex.style.formatting.unsrt.date = pybtex.style.template.words(sep="")[
    "(", pybtex.style.template.field("year"), ")"
]


class ApsStyle(pybtex.style.formatting.unsrt.Style):
    """Style that mimicks APS journals."""

    def __init__(
        self,
        label_style=None,
        name_style=None,
        sorting_style=None,
        abbreviate_names=True,
        min_crossrefs=2,
        **kwargs
    ):
        super().__init__(
            label_style=label_style,
            name_style=name_style,
            sorting_style=sorting_style,
            abbreviate_names=abbreviate_names,
            min_crossrefs=min_crossrefs,
            **kwargs
        )

    def format_title(self, e, which_field, as_sentence=True):
        """Set titles in italics."""
        formatted_title = pybtex.style.template.field(
            which_field, apply_func=lambda text: text.capitalize()
        )
        formatted_title = pybtex.style.template.tag("em")[formatted_title]
        if as_sentence:
            return pybtex.style.template.sentence[formatted_title]
        else:
            return formatted_title

    def get_article_template(self, e):
        volume_and_pages = pybtex.style.template.first_of[
            # volume and pages
            pybtex.style.template.optional[
                pybtex.style.template.join[
                    " ",
                    pybtex.style.template.tag("strong")[
                        pybtex.style.template.field("volume")
                    ],
                    ", ",
                    pybtex.style.template.field(
                        "pages",
                        apply_func=pybtex.style.formatting.unsrt.dashify,
                    ),
                ],
            ],
            # pages only
            pybtex.style.template.words[
                "pages",
                pybtex.style.template.field(
                    "pages", apply_func=pybtex.style.formatting.unsrt.dashify
                ),
            ],
        ]
        template = pybtex.style.formatting.toplevel[
            self.format_names("author"),
            self.format_title(e, "title"),
            pybtex.style.template.sentence(sep=" ")[
                pybtex.style.template.field("journal"),
                pybtex.style.template.optional[volume_and_pages],
                pybtex.style.formatting.unsrt.date,
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_book_template(self, e):
        template = pybtex.style.formatting.toplevel[
            self.format_author_or_editor(e),
            self.format_btitle(e, "title"),
            self.format_volume_and_series(e),
            pybtex.style.template.sentence(sep=" ")[
                pybtex.style.template.sentence(add_period=False)[
                    pybtex.style.template.field("publisher"),
                    pybtex.style.template.optional_field("address"),
                    self.format_edition(e),
                ],
                pybtex.style.formatting.unsrt.date,
            ],
            pybtex.style.template.optional[
                pybtex.style.template.sentence[self.format_isbn(e)]
            ],
            pybtex.style.template.sentence[
                pybtex.style.template.optional_field("note")
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_incollection_template(self, e):
        template = pybtex.style.formatting.toplevel[
            pybtex.style.template.sentence[self.format_names("author")],
            self.format_title(e, "title"),
            pybtex.style.template.words[
                "In",
                pybtex.style.template.sentence[
                    pybtex.style.template.optional[
                        self.format_editor(e, as_sentence=False)
                    ],
                    self.format_btitle(e, "booktitle", as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    self.format_chapter_and_pages(e),
                ],
            ],
            pybtex.style.template.sentence(sep=" ")[
                pybtex.style.template.sentence(add_period=False)[
                    pybtex.style.template.optional_field("publisher"),
                    pybtex.style.template.optional_field("address"),
                    self.format_edition(e),
                ],
                pybtex.style.formatting.unsrt.date,
            ],
            self.format_web_refs(e),
        ]
        return template


pybtex_register_plugin("pybtex.style.formatting", "apsstyle", ApsStyle)


# -- Monkeypatch for instance attribs (sphinx bug #2044) -----------------------


def iad_add_directive_header(self, sig):
    ClassLevelDocumenter.add_directive_header(self, sig)


InstanceAttributeDocumenter.add_directive_header = iad_add_directive_header

# -- Options for HTML output ---------------------------------------------------

# on_rtd is whether we are on readthedocs.org, this line of code grabbed from
# docs.readthedocs.org
on_rtd = os.environ.get("READTHEDOCS", None) == "True"

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme

    html_theme = "sphinx_rtd_theme"
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# html_theme = 'sphinxdoc'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "collapse_navigation": True,
    "display_version": True,
    "navigation_depth": 4,
}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = 'favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# JavaScript filenames, relative to html_static_path
html_js_files = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

doctr_versions_menu_conf = {'menu_title': 'Docs'}

# -- Options for nbsphinx -0---------------------------------------------------

nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='docs') %}

.. only:: html

    .. role:: raw-html(raw)
        :format: html

    :raw-html:`<a href="http://nbviewer.jupyter.org/github/qucontrol/krotov/blob/<<GIT_TAG>>/{{ docname }}" target="_blank"><img alt="Render on nbviewer" src="https://img.shields.io/badge/render%20on-nbviewer-orange.svg" style="vertical-align:text-bottom"></a>&nbsp;<a href="https://mybinder.org/v2/gh/qucontrol/krotov/<<GIT_TAG>>?filepath={{ docname }}" target="_blank"><img alt="Launch Binder" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>`
""".replace(
    "<<GIT_TAG>>", git_tag
)


# -----------------------------------------------------------------------------
def setup(app):
    app.connect("builder-inited", run_apidoc)
    app.connect("builder-inited", generate_patched_readme)
