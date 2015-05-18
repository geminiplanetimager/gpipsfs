import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import poppy
import scipy.interpolate

from .dms import GPITweeter


class AOTelemetryWFE(GPITweeter):
    """ GPI AO residuals from telemetry

    Parameters
    --------------
    telemetrypath : string
        File path to some _when_* GPI telemetry file. This can either
        be a full filename for a specific FITS file or just the base
        _When_<timestamp> part that is common to all of the files in a
        given set.
    index: int
        Index for which wavefront error slice to use

    """
    def __init__(self, telemetrypath, index=0):
        telemdir, telemfile = os.path.split(telemetrypath)
        telemwhen = "_".join(telemfile.split('_')[0:4])
        self.telem_base = os.path.join(telemdir,telemwhen)+"_"

        self._telem_header = fits.getheader(self.telembase+"phase.fits")
        self._subapmask = fits.getdata(self.telembase+"subapmask.fits")
        self._inten = fits.getdata(self.telembase+"inten.fits")
        self._phase = fits.getdata(self.telembase+"phase.fits")
        self._ttdata = fits.getdata(self.telembase+"ttdata.fits")
        self._hz = fits.getdata(self.telembase+"hz.fits")

        # Some derived quantities, based on gpilib load_telem.pro
        self._pingrid = (self._subapmask + np.roll(self._subapmask,0,1)
                         +np.roll(self._subapmask,1,0)+np.roll(self._subapmask,1,1))>1

        self.frame_rate = np.round(np.abs(hz).max()*2/500)*500  # round to 500 or 1000

        self._mean_inten = self.inten.mean(axis=0).shape
        #xphase =

