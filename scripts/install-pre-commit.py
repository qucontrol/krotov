#!/usr/bin/env python
"""Install the pre-commit hooks

This should be run from the Makefile immediately after creating .venv/py36
"""
import os


if not os.path.isdir(".venv/py36"):
    os.system("make .venv/py36/bin/python")

if os.path.isdir(".git") and not os.path.isfile(".git/hooks/pre-commit"):
    print("##################################################################")
    print("Installing pre-commit hooks")
    print("##################################################################")
    os.system("./.venv/py36/bin/pre-commit install")
