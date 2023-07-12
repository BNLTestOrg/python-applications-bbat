#!/usr/bin/env python3

from setuptools import setup, find_packages


APP_NAME = "bbat"
PACKAGE_NAME = "bbat"

setup(
    name=APP_NAME,
    use_scm_version=dict(write_to=f"{PACKAGE_NAME}/version.py"),
    packages=find_packages(),
    setup_requires=["setuptools_scm"],
    entry_points={
        "console_scripts": [f"{APP_NAME}={PACKAGE_NAME}.__main__:main"],
    },
    include_package_data=True,
)
