#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from setuptools import setup, find_packages
setup(
    name = "gpipsfs",
    version = "0.2.0",
    packages = ["gpipsfs"],
    install_requires=['numpy>=1.10.0',
        'astropy>=2.0.0',
        'poppy>=0.7.0',
        'matplotlib>=1.3.0',
        'pysynphot>=0.9' ],
    package_data = {
        'gpipsfs': ['data/*.txt', 'data/*.fits'],
    },
    # metadata for upload to PyPI
    author = "Gemini Planet Imager Team",
    author_email = "geminiplanetimager@gmail.com",
    description = "GPI PSF Simulator",
    license = "BSD",
    url = "http://github.com/geminiplanetimager/gpipsfs/",
)
