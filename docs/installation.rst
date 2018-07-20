
Installation and Setup
=======================

Prerequisites
-------------

The easiest method is to use Conda to install the following:

 * numpy, scipy, matplotlib
 * astropy
 * poppy
 * pysynphot

Recent / current versions of all are recommended. 

Installing the GPI Data Pipeline
--------------------------------

``gpipsfs`` requires having the GPI data pipeline installed so that it can
access some of the configuration files, in particular reference data on the
filter transmission curves. The pipeline doesn't need to be operational on
your machine, in fact; you don't need to have IDL installed either. 
At minimum it just needs to be able to access the pipeline's ``config`` subdirectory
under the directory path indicated by the environment variable ``$GPI_DRP_DIR``. 

See the GPI Data Pipeline documentation for download links and installation instructions. 


Installing gpipsfs
------------------

To clone this repo from Github and install it::

  % git clone https://github.com/geminiplanetimager/gpipsfs
  % cd gpipsfs
  % python setup.py install

Equivalently you can do the same in one line using ``pip``::

  % pip install -e git+https://github.com/geminiplanetimager/gpipsfs.git#egg=gpipsfs

This will create a directory ``./src/gpipsfs`` in your current directory containing the cloned code. 

