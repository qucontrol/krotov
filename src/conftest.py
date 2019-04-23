"""Set up the environment for doctests

This file is automatically evaluated by py.test. It ensures that we can write
doctests without distracting import statements in the doctest.
"""
import inspect
from collections import OrderedDict

import numpy
import pytest

import krotov


@pytest.fixture(autouse=True)
def set_doctest_env(doctest_namespace):
    doctest_namespace['numpy'] = numpy
    doctest_namespace['krotov'] = krotov
    doctest_namespace['inspect'] = inspect
    doctest_namespace['OrderedDict'] = OrderedDict
