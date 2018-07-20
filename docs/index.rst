
gpipsfs : A GPI PSF Simulation Toolkit
======================================================


.. _intro:

This Python package produces simulated PSFs for the Gemini Planet
Imager, a facility-class exoplanets imaging instrument at Gemini
South. 

This code provides a **toy model** of GPI intended primarily for understanding how the GPI coronagraph optics work; it is **not** an adaptive optics simulator, does not model wavefront errors, operates in a simplified optical approximation regime (Fraunhofer diffraction, not Fresnel), and is not intended to be a high fidelity model of GPI. 


Most of the documentation of this package is provided as Jupyter notebooks. These are included below in static forms, but
you are encourages to run them directly as notebooks in Jupyter yourself. 

You may also wish to check out the documentation for `POPPY <http://poppy-optics.readthedocs.org/>`_, 
the optics simulation framework that this is built on top of. 

.. caution:: 
   This software comes with even less warranty or guarantee of correctness than the data pipeline does... 
   

Contents
--------

.. toctree::
   :maxdepth: 1
   
   installation.rst
   usage.rst
   Getting Started with GPI PSFs.ipynb
   GPI Tutorial on Coronagraphs.ipynb
   Example Coronagraph Misalignments.ipynb
   Sat Spots and Sine Waves.ipynb

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`

Documentation last updated on |today|

