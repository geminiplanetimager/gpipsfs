
Usage Instructions will Eventually Go Here
============================================


For now here is a basic example::


    >>> import gpipsfs
    >>> gpi = gpipsfs.GPI()
    >>> gpi.obsmode = "J_coron"
    >>> psf = gpi.calcPSF(display=True)
    >>> psf.writeto('my_J_psf.fits')


The GPI class has attributes for the filter, apodizer, occulter, and Lyot mask.
You can set them indivdually, or just set the obsmode to control the entire set
at once in the same manner as the actual instrument. 


Results are returned as FITS files. 

See the documentation for POPPY or WebbPSF for more general background on how
this code works. 
