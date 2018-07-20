gpipsfs: GPI Point Spread Function (PSF) Simulation Toolkit
================================================================

This Python package produces simulated PSFs for the Gemini Planet
Imager, a facility-class exoplanets imaging instrument at Gemini
South.

This code provides a **toy model** of GPI intended primarily for understanding how the GPI coronagraph optics work; it is **not** an adaptive optics simulator, does not model wavefront errors, operates in a simplified optical approximation regime (Fraunhofer diffraction, not Fresnel), and is not intended to be a high fidelity model of GPI.


In lieu of more complete documentation, see
http://nbviewer.ipython.org/github/geminiplanetimager/gpipsfs/blob/master/notebooks/Getting%20Started%20with%20GPI%20PSFs.ipynb


Requirements & Installation
----------------------------------------

Prerequisites:
 * numpy, scipy, matplotlib, etc.
 * astropy
 * pysynphot
 * poppy (>= version 0.7.0)

And also
 * A copy of the GPI data reduction pipeline, with the environment variable $GPI_DRP_DIR configured to
   point to its location.


To install gpipsfs, clone this repo from Github, and then::

    > cd to that folder
    > python setup.py install


Getting Started
------------------

Check out the ipython notebooks, in particular the 'Getting Started" one.


Warnings, Caveats, and Disclaimers
---------------------------------------

Work in progress, incomplete, simplified code, no guarantees, etc! 

Currently no wavefront error terms included.

Cash prize available for the first person to find a sign error somewhere in this code.

To reiterate: This code provides a **toy model** of GPI intended primarily for understanding how the GPI coronagraph optics work; it is **not** an adaptive optics simulator, does not model wavefront errors, operates in a simplified optical approximation regime (Fraunhofer diffraction, not Fresnel), and is not intended to be a high fidelity model of GPI.
