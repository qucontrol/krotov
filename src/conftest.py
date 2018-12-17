"""Set up the environment for doctests

This file is automatically evaluated by py.test. It ensures that we can write
doctests without distracting import statements in the doctest.
"""
import numpy
import krotov
from collections import OrderedDict

# noinspection PyPackageRequirements
import pytest


@pytest.fixture(autouse=True)
def set_doctest_env(doctest_namespace):
    doctest_namespace['numpy'] = numpy
    doctest_namespace['krotov'] = krotov
    doctest_namespace['OrderedDict'] = OrderedDict
