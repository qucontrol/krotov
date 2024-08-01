"""Test the dump_result convergence routine"""

import copy
import os

import krotov
import pytest


oct_result_name = "oct_result.dump"


def incl_range(a, b, step=1):
    e = 1 if step > 0 else -1
    return range(a, b + e, step)


# Note: The `oct_result.dump` file was manually created from
# docs/notebooks/01_example_simple_state_to_state.ipynb (the `opt_result`
# variable)


@pytest.fixture
def oct_result_18(request):
    testdir = os.path.splitext(request.module.__file__)[0]
    filename = os.path.join(testdir, oct_result_name)
    return krotov.result.Result.load(filename)


@pytest.fixture
def oct_result_5(oct_result_18):
    result = copy.deepcopy(oct_result_18)
    result.iters = list(incl_range(0, 5))
    return result


@pytest.fixture
def oct_result_6(oct_result_18):
    result = copy.deepcopy(oct_result_18)
    result.iters = list(incl_range(0, 6))
    return result


@pytest.fixture
def oct_result_10(oct_result_18):
    result = copy.deepcopy(oct_result_18)
    result.iters = list(incl_range(0, 10))
    return result


def test_invalid_dump_result(oct_result_18):
    with pytest.raises(ValueError):
        krotov.convergence.dump_result(oct_result_name, every=0)
    filename = "non_existent_folder/oct_result.dump"
    # Use every=18 because this result only has 18 iterations.
    # The default every=10 returns None.
    dump_result = krotov.convergence.dump_result(filename, every=18)
    response = dump_result(oct_result_18)
    assert "Could not store" in response


def test_dump_result_overwrite(oct_result_5, oct_result_10, tmpdir):
    filename = str(tmpdir.join("oct_result_a.dump"))
    assert not os.path.isfile(filename)
    dump_result = krotov.convergence.dump_result(filename, every=5)

    response = dump_result(oct_result_5)
    assert response is None
    assert os.path.isfile(filename)
    result = krotov.result.Result.load(filename)
    assert result.iters[-1] == 5

    response = dump_result(oct_result_10)
    assert response is None
    assert os.path.isfile(filename)
    result = krotov.result.Result.load(filename)
    assert result.iters[-1] == 10


def test_dump_result_keep(oct_result_5, oct_result_6, oct_result_10, tmpdir):
    filename = str(tmpdir.join("oct_result_{iter:06d}.dump"))
    assert not os.path.isfile(filename)
    dump_result = krotov.convergence.dump_result(filename, every=5)

    response = dump_result(oct_result_5)
    assert response is None
    assert os.path.isfile(str(tmpdir.join("oct_result_000005.dump")))

    response = dump_result(oct_result_6)
    assert response is None
    assert not os.path.isfile(str(tmpdir.join("oct_result_000006.dump")))

    response = dump_result(oct_result_10)
    assert response is None
    assert os.path.isfile(str(tmpdir.join("oct_result_000010.dump")))
