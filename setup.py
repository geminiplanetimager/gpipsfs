#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from setuptools import setup, find_packages
setup(
    name = "gpipsfs",
    version = "0.1.0",
    packages = ["gpipsfs"],
    install_requires=['numpy>=1.8.0', 'astropy>=1.0.0', 'poppy>=0.3.4'],
    package_data = {
        'gpipsfs': ['data/*.txt'],
    },
    # metadata for upload to PyPI
    author = "Gemini Planet Imager Team",
    author_email = "geminiplanetimager@gmail.com",
    description = "GPI PSF Simulator",
    license = "BSD",
    url = "http://github.com/geminiplanetimager/gpipsfs/",
)
