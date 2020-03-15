.. include:: _README.rst

.. The file _README.rst is a patched copy of ../README.rst. The patch ./_README.patch is applied in conf.py.

.. _CitingKrotov:

Citing the Krotov Package
-------------------------


.. attention::

   Please cite the :mod:`krotov` package as

   * \M. H. Goerz *et al*., *Krotov: A Python implementation of Krotov's method for quantum optimal control*, `SciPost Phys. 7, 080 <https://scipost.org/SciPostPhys.7.6.080>`_ (2019)

You can also print this from ``krotov.__citation__``:

.. doctest::

   >>> print(krotov.__citation__)
   M. H. Goerz et al., Krotov: A Python implementation of Krotov's method for quantum optimal control, SciPost Phys. 7, 080 (2019)

The corresponding BibTeX entry is available in ``krotov.__bibtex__``:

.. doctest::

   >>> print(krotov.__bibtex__)
   @article{GoerzSPP2019,
       author = {Michael H. Goerz and Daniel Basilewitsch and Fernando Gago-Encinas and Matthias G. Krauss and Karl P. Horn and Daniel M. Reich and Christiane P. Koch},
       title = {Krotov: A {Python} implementation of {Krotov's} method for quantum optimal control},
       journal={SciPost Phys.},
       volume={7},
       pages={80},
       year={2019},
       doi={10.21468/SciPostPhys.7.6.080},
   }
