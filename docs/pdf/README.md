# PDF documentation for the `krotov` package


The documentation for different versions of the
[krotov Python package][package] is available in PDF format through the
versions menu in the [online documentation][] or on the [Releases page][].

For non-released development versions, a PDF of the documentation can be
generated from a git checkout of [`krotov`][package]. Run

    tox -e bootstrap
    tox -e docs -- _build/tex -b latex
    cp docs/*.pdf docs/_build/tex/
    tox -e run-cmd -- python docs/build_pdf.py

or simply `make docs-pdf` to create a file `docs/_build/tex/krotov.pdf`.

This assumes that you have `lualatex` (e.g. through [texlive][] on Linux/macOS
or [MikTeX][] on Windows). You must also have the [DejaVu fonts][] installed on
your system.


[package]: https://github.com/qucontrol/krotov
[online documentation]: https://qucontrol.github.io/krotov/
[Releases page]: https://github.com/qucontrol/krotov/releases
[texlive]: https://www.tug.org/texlive/
[MikTex]: https://miktex.org
[DejaVu fonts]: https://dejavu-fonts.github.io
