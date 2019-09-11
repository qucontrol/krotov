"""This file is automatically executed by pytest when testing anything in the
docs folder"""


def pytest_collectstart(collector):
    """Ignore stderr and javascript output when verifying notebooks.

    Works around a test failure on Travis/AppVeyor (See Travis Build 56).

    https://nbval.readthedocs.io/en/latest/#Skipping-certain-output-types
    """
    if collector.__class__.__name__ == 'IPyNbFile':
        collector.skip_compare += (
            'application/javascript',
            'stderr',
            'application/vnd.jupyter.widget-view+json',
        )
