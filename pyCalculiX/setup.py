#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import setuptools

long_description = """CalculiX-Python Utilities"""

setuptools.setup(
    name="pycalculix",
    version="0.0.1",
    author="Marin Lauber",
    author_email="marinlauber@gmail.com",
    description="CalculiX-Python Utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marinlauber/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)