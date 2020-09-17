
import matplotlib.pyplot as plt
import numpy as np
import poppy

def bilinear_interpolate(im, x, y):
    """ Simply and efficient (fast) bilinear interpolation
    From http://stackoverflow.com/questions/12729228/simple-efficient-bilinear-interpolation-of-images-in-numpy-and-python


    Parameters
    ------------
    im : ndarray
        2d image to interpolate
    x, y : ndarrays
        indices of coordinations for each point
        desired in the output array, in units of
        fraction pixels into the input image 
        to be interpolated. 
    """

    x = np.asarray(x)
    y = np.asarray(y)

    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, im.shape[1]-1);
    x1 = np.clip(x1, 0, im.shape[1]-1);
    y0 = np.clip(y0, 0, im.shape[0]-1);
    y1 = np.clip(y1, 0, im.shape[0]-1);

    Ia = im[ y0, x0 ]
    Ib = im[ y1, x0 ]
    Ic = im[ y0, x1 ]
    Id = im[ y1, x1 ]

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    return wa*Ia + wb*Ib + wc*Ic + wd*Id




def plot_contrast(coron_psf, direct_psf, overplot=False, label=None, ext=0, ax=None, **kwargs):
    """Plot a PSF, normalized into contrast units relative to a direct unocculted PSF
    """

    cor_radius, cor_profile = poppy.radial_profile(coron_psf, ext=ext, **kwargs)

    direct_peak = direct_psf[ext].data.max()

    cor_profile /= direct_peak

    if ax is None:
        ax = plt.gca()

    ax.semilogy(cor_radius, cor_profile, label=label, **kwargs)

    if not overplot:
        ax.set_xlabel("Radius [arcsec]")
        ax.set_ylabel("Contrast")



