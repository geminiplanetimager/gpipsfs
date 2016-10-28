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
 * poppy (you may need a development version. Install from Github)




To install gpipsfs, either clone this repo from Github, or
Zip files are available from https://github.com/geminiplanetimager/gpipsfs. see link at right to 'Download ZIP'
Download the zip file and unzip. And then::

    > cd to that folder

    > python setup.py install
    
    


Getting Started
------------------

Check out the ipython notebooks, in particular the 'Getting Started" one. 

Depending on which version of IPython (now a.k.a. Jupyter) you have installed, 
you make have trouble opening some of the notebooks. If you run into this problem, 
either: 

1. Install a recent version of IPython, 3.0 or greater. (The latest version of 
   the Anaconda scientific python distribution includes this, for instance.). Or else
2. Try the filenames with the "v3" in them, which are in "notebook format 3". 
   (Confusingly, ipython 3 uses notebook format 4, while ipython 2 uses file format 3!) 
   If that still doesn't work for you, ask Marshall for help.  


Warnings, Caveats, and Disclaimers
---------------------------------------

Work in progress, incomplete, simplified code, no guarantees, etc!  

Currently no wavefront error terms included. 

Cash prize available for the first person to find a sign error somewhere in this code. 

To reiterate: This code provides a **toy model** of GPI intended primarily for understanding how the GPI coronagraph optics work; it is **not** an adaptive optics simulator, does not model wavefront errors, operates in a simplified optical approximation regime (Fraunhofer diffraction, not Fresnel), and is not intended to be a high fidelity model of GPI. 
