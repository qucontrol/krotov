#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages


def get_version(filename):
    """Extract the package version"""
    with open(filename) as in_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                return line.split('=')[1].strip()[1:-1]
    raise ValueError("Cannot extract version from %s" % filename)


with open('README.rst') as readme_file:
    readme = readme_file.read()

try:
    with open('HISTORY.rst') as history_file:
        history = history_file.read()
except OSError:
    history = ''

# requirements for use
requirements = ['attrs', 'glom', 'numpy', 'scipy', 'qutip']

# requirements for development (testing, generating docs)
dev_requirements = [
    'jupyter', 'coverage', 'pytest', 'pytest-cov', 'pytest-xdist', 'nbval',
    'twine', 'pep8', 'flake8', 'wheel', 'sphinx', 'sphinx-autobuild',
    'sphinx_rtd_theme', 'nbsphinx', 'matplotlib', 'gitpython', 'watermark',
    'sphinxcontrib-bibtex', 'weylchamber', 'click']
dev_requirements.append('better-apidoc')

# some recommended packages that make development nicer
dev_extras = [
    'jupyterlab', 'pdbpp']

version = get_version('./src/krotov/__init__.py')

setup(
    author="Michael Goerz",
    author_email='mail@michaelgoerz.net',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Python implementation of Krotov's method for quantum optimal control",
    install_requires=requirements,
    extras_require={
        'dev': dev_requirements,
        'extras': dev_extras,
    },
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='krotov',
    name='krotov',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    url='https://github.com/qucontrol/krotov',
    version=version,
    zip_safe=False,
)
