"""Test that we can serialize a Resuls object"""

import logging
import os
import pickle

import numpy as np

from krotov.result import Result


def test_serialization_roundtrip(request, tmpdir, caplog):
    """Test load/dump of Result"""
    testdir = os.path.splitext(request.module.__file__)[0]
    with caplog.at_level(logging.WARNING):
        result = Result.load(os.path.join(testdir, 'oct_result.dump'))
    assert 'Result.objectives contains control placeholders' in caplog.text
    assert isinstance(result, Result)
    dumpfile = str(tmpdir.join('oct_result.dump'))
    result.dump(dumpfile)
    with open(dumpfile, 'rb') as dump_fh:
        result2 = pickle.load(dump_fh)
    for name in result.__dict__:
        val1 = getattr(result, name)
        val2 = getattr(result2, name)
        _check_recursive_equality(val1, val2)


def test_serialization_finalize(request, caplog):
    testdir = os.path.splitext(request.module.__file__)[0]
    dumpfile = os.path.join(testdir, 'oct_result_incomplete.dump')
    with caplog.at_level(logging.WARNING):
        result = Result.load(dumpfile)
    assert "not finalized" in caplog.text
    nt = len(result.tlist)
    assert len(result.optimized_controls[0]) == nt - 1
    result = Result.load(dumpfile, finalize=True)
    assert len(result.optimized_controls[0]) == nt


def test_serialization_broken(request, caplog):
    testdir = os.path.splitext(request.module.__file__)[0]
    dumpfile = os.path.join(testdir, 'oct_result_broken.dump')
    with caplog.at_level(logging.ERROR):
        Result.load(dumpfile)
    assert "incongruent" in caplog.text


def _check_recursive_equality(val1, val2):
    if isinstance(val1, list):
        for (v1, v2) in zip(val1, val2):
            _check_recursive_equality(v1, v2)
    elif isinstance(val1, np.ndarray):
        assert np.all(val1 == val2)
    else:
        assert val1 == val2
