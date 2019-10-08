# Examples

This directory contains examples scripts and snippets from the paper

> M. H. Goerz et al., *Krotov: A Python implementation of Krotov's method for quantum optimal control*, [SciPost Phys. 7, 080][DOI] (2019)

The following files are included:

* [`tls_state_to_state.py`][tls]: Complete example script for a simple state-to-state transition in a two level system, see Section 2.4 of the paper. A Jupyter notebook version of this example is available as part of the documentation as [`01_example_simple_state_to_state.ipynb`][tlsnb].
* `quantum_gate_hilbert.py` and `quantum_gate_liouville.py`: Snippets for the optimization towards a quantum gate, see Section 3.2 of the paper. Full examples are available as Jupyter notebooks [`05_example_transmon_xgate.ipynb`][gatenb] and [` 06_example_3states.ipynb`][lvgatenb].
* `ensemble.py`: Snippet for an ensemble optimization, see Section 3.3 of the paper. A full example is available as [`08_example_ensemble.ipynb`][ensemblenb].
* `second_order.py`: Snippet for a second-order Krotov optimization, see Section 3.4 of the paper. A full example is available as [`07_example_PE.ipynb`][penb].

The `.out` files contain the expected output of each Python script.

Only the first example script, [`tls_state_to_state.py`][tls], is intended to be used for study.

The remaining files only contain snippets used in the corresponding discussion in the paper, together with minimal setup and test code that verifies the snippet's correctness. This additional code should not be used as an example for the proper use of the `krotov` package. Instead, refer to the Jupyter notebook files listed above.

For Section 3.1 of the paper (optimization with complex-valued control fields), which does not contain any code snippets, there is a Jupyter notebook [`02_example_lambda_system_rwa_complex_pulse.ipynb`][complexnb] in the online documentation.

These examples are tested against the `v1.0.0` release of the `krotov` package.


[arxiv]: https://arxiv.org/abs/1902.11284
[DOI]: https://doi.org/10.21468/SciPostPhys.7.6.080
[tls]: tls_state_to_state.py
[tlsnb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/01_example_simple_state_to_state.ipynb
[complexnb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb
[gatenb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/05_example_transmon_xgate.ipynb
[lvgatenb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/06_example_3states.ipynb
[ensemblenb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/08_example_ensemble.ipynb
[penb]: https://mybinder.org/v2/gh/qucontrol/krotov/v1.0.0?filepath=docs/notebooks/07_example_PE.ipynb
