# Krotov: A Python implementation of Krotov's method for quantum optimal control

[![Github](https://img.shields.io/badge/github-qucontrol/krotov-blue.svg)][package]
[![DOI](https://img.shields.io/badge/DOI-10.21468/SciPostPhys.7.6.080-blue.svg)][DOI]
[![arXiv](https://img.shields.io/badge/arXiv-1902.11284-red.svg)][arxiv]

This is the LaTeX source of the paper

> M. H. Goerz et al., *Krotov: A Python implementation of Krotov's method for quantum optimal control*, [SciPost Phys. 7, 080][DOI] (2019)

that describes `v1.0.0` of the [krotov Python package][package].

[arxiv]: https://arxiv.org/abs/1902.11284
[package]: https://github.com/qucontrol/krotov
[DOI]: https://doi.org/10.21468/SciPostPhys.7.6.080

## Compilation

Run `make` to generate the paper from source.

Run `make arxiv.tgz` to create an arXiv submission.

Run `make clean`/`make distclean` to remove all compiled files and the virtual environment that was used for compilation.


## Examples

The example script used in Section 2.4 of the paper can be found in the [examples folder](examples).

The folder also contains example snippets used in Section 3 of the paper.

You can verify that the examples produce the expected output by running `make test`.
