import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import poppy
import scipy.interpolate

from .main import GeminiPrimary

# Classes for dealing with AO Telemetry sets

class GPI_Globals(object):
    """ Container for same constants as gpilib's gpi_globals,
    with same variable names to ease porting of code. Plus some
    other variables as needed."""
    gpi_tweet_n = 48
    gpi_woof_n = 9
    gpi_numacross=43.2
    gpi_tweet_spacing = GeminiPrimary.primary_diameter/gpi_numacross
    gpi_woof_spacing  = GeminiPrimary.primary_diameter/gpi_numacross*5.5
    # below ones are not in gpi_globals
    pupil_center_subap = 23
    pupil_center_tweeter=23.5
    pupil_center_woofer=4


class DeformableMirror(poppy.AnalyticOpticalElement):
    """ Generic deformable mirror, of the continuous face sheet variety"""
    def __init__(self, shape=(10,10)):
        poppy.OpticalElement.__init__(self, planetype=poppy.poppy_core._PUPIL)
        self._shape = shape                  # number of actuators
        self.surface = np.zeros(shape)      # array for the DM surface WFE
        self.numacross = shape[0]           # number of actuators across diameter of
                                            # the optic's cleared aperture (may be
                                            # less than full diameter of array)
        self.actuator_spacing = 1.0/self.numacross  # distance between actuators,
                                                    # projected onto the primary
        self.pupil_center = (shape[0]-1.)/2 # center of clear aperture in actuator units
                                            # (may be offset from center of DM)
        self._stale_interpolator = True
    @property
    def shape(self):
        return self._shape

    def set_surface(self, new_surface, units='nm'):
        assert new_surface.shape == self.shape
        assert units=='nm'
        self.surface = np.asarray(new_surface, dtype=float)*1e-9

        self._stale_interpolator = True

    def _update_interpolator(self):
        """ Update the interpolator functions used for
        determining surface state at arbitrary samplings
        based on the DM surface
        """

        yc, xc = self.get_coordinates(one_d=True)
        self._stale_interpolator = False

        # TODO evaluate other interpolation options
        # stuff in scipy.ndimage.interpolation maybe?
        # This is clearly NOT the right way to do this! Given extreme slowness
        self._interp_surface = scipy.interpolate.RectBivariateSpline(xc, yc, self.surface)

    def set_actuator(self, actx, acty, new_value, units='nm'):
        # FIXME do something more comprehensive with units
        assert units=='nm'
        if actx < 0 or actx > self.shape[1]-1:
            raise ValueError("X axis coordinate is out of range")
        if acty < 0 or acty > self.shape[0]-1:
            raise ValueError("Y axis coordinate is out of range")
        self.surface[acty, actx] = new_value*1e-9
        self._stale_interpolators = True

    def get_coordinates(self, one_d=False):
        """ Y and X coordinates for the actuators

        Parameters
        ------------
        one_d : bool
            Return 1-dimensional arrays of coordinates per axis?
            Default is to return 2D arrays with same shape as full array.
        """

        y_act = (np.arange(self.shape[0])-self.pupil_center)*self.actuator_spacing
        x_act = (np.arange(self.shape[1])-self.pupil_center)*self.actuator_spacing

        if not one_d: # convert to 2D
            y_act.shape = (self.shape[0],1)
            y_act = y_act * np.ones( (1, self.shape[1]))

            x_act.shape = (1, self.shape[1])
            x_act = x_act * np.ones( (self.shape[0], 1))

        return y_act, x_act


    def annotate(self, marker='o', **kwargs):
        """ Overplot actuator coordinates on some already-existing pupil display

        """
        yc, xc = self.get_coordinates()
        ax = plt.gca()

        # jump through some hoops to avoid autoscaling the X,Y coords
        # of the prior plot here, but retain the autoscale state
        autoscale_state = (ax._autoscaleXon, ax._autoscaleYon)
        ax.autoscale(False)
        plt.scatter(xc, yc, marker=marker, **kwargs)
        ax._autoscaleXon, ax._autoscaleYon = autoscale_state

    def getPhasor(self,wave=None):
        """ Interpolate from the current optic surface state onto the
        desired coordinates.

        CAUTION: This right now uses fairly simple interpolation and does not
        take into account the actual actuator influence functions.
        That's a project for another day.


        if no wave object is supplied, will return the phasor on the 
        native coordinates of the DM, for a wavelength of 1 micron
        """

        if wave is None:
            wavelen = 1e-6
            phasor = np.exp(1.j * 2 * np.pi * self.surface/wavelen)
        else:
            if self._stale_interpolator: self._update_interpolator()

            y, x = wave.coordinates()
            interpolated_surface = self._interp_surface(x,y)

            phasor = np.exp(1.j * 2 * np.pi * interpolated_surface/wave.wavelength)
        return phasor


    def display(self, interpolated=False, annotate=False, what='phase', *args, **kwargs):

        if interpolated:
            # Display the interpolated optic on a finer sampling
            if self._stale_interpolator: self._update_interpolator()
            returnvalue = poppy.AnalyticOpticalElement.__display__(self, what=what, *args, **kwargs)
        else:
            # display in DM coordinates
            # temporarily set attributes appropriately as if this were a regular OpticalElement
            self.amplitude = np.ones_like(self.surface)
            self.opd = self.surface 
            self.pixelscale = self.actuator_spacing

            #then call parent class display
            returnvalue = poppy.OpticalElement.display(self,  what=what, **kwargs)

            # now un-set all the temporary attributes, since this is analytic and
            # these are unneeded
            del self.pixelscale
            del self.opd
            del self.amplitude


        if annotate: self.annotate()
        return returnvalue

class GPITweeter(DeformableMirror):
    def __init__(self, mems_print_through=True):
        DeformableMirror.__init__(self, shape=(GPI_Globals.gpi_tweet_n, GPI_Globals.gpi_tweet_n))
        self.name = "GPI Tweeter"
        self.numacross = GPI_Globals.gpi_numacross
        self.actuator_spacing = GPI_Globals.gpi_tweet_spacing
        self.pupil_center = GPI_Globals.pupil_center_tweeter
        self.pupil_diam = GPI_Globals.gpi_tweet_n*GPI_Globals.gpi_tweet_spacing   # for display, projected full area around 48x48 subaps

        self.mems_print_through = mems_print_through
        self._mems_print_through_amplitude = 15e-9
        # 15 nm, estimated from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf

        my_path=os.path.abspath(os.path.dirname(__file__))
        self._actuator_type_info = fits.open(os.path.join(my_path, 'data','GPI_tweeter_actuators.fits'))

    @property
    def bad_actuators(self):
        """Returns a list of coordinate indices for the actuators which are
        nonoperable """
        act_map = self._actuator_type_info
        wflagged = np.where(  ( act_map[0].data == act_map[0].header['DEAD'] |
                                act_map[0].data == act_map[0].header['WEAK']) )
        return wflagged

    def getPhasor(self, wave):
        phasor = DeformableMirror.getPhasor(self,wave)

        if self.mems_print_through:
            mems_print_through_phasor= self._getPhasor_MEMS_print_through(self,wave)
            phasor *= mems_print_through_phasor
        return phasor

    def _getPhasor_MEMS_print_through(self,wave):
        """ DM surface print through """

        # GPI tweeter actuators are reimaged to 18 cm subapertures

        # Boston DM print through info in:
        #  http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf
        #  ao4elt.lesia.obspm.fr/sites/ao4elt/IMG/ppt/Bifano.ppt

        # in horizontal direction, the print through is about 35/190 pixels = 18% of the width
        # in the vertical direction, closer to 31%, but it's more like 2 narrow bars each 10% wide
        # and there's a 10% wide dot in the middle of it too
        #printthrough properties:
        pt_col_width = 0.18
        # 15 nm, estimated from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf
        pt_col_value = self._mems_print_through_amplitude
        pt_row_width = 0.10
        pt_row_value = -1 * self._mems_print_through_amplitude

        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._PUPIL)

        opd = np.ones(wave.shape)
        y, x = wave.coordinates()


        pixscale = x[0,1] - x[0,0]
        opd[np.mod(x,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_row_value
        opd[np.mod(y,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_col_value

        phasor = np.exp(1.j * 2 * np.pi * opd/wave.wavelength)
        return phasor

    def annotate(self, markbad=True, badmarker='o', marker='+', **kwargs):
        # first plot all the normal ones
        DeformableMirror.annotate(self, marker=marker, **kwargs)

        if markbad:
            # now the less-than-good ones
            yc, xc = self.get_coordinates()
            ax = plt.gca()
            autoscale_state = (ax._autoscaleXon, ax._autoscaleYon)
            ax.autoscale(False)
            act_map = self._actuator_type_info
            for act_type, color in zip(['DEAD', 'COUPLED', 'WEAK','VARIABLE'],
                                ['red',  'orange', 'brown', 'magenta']):
                wflagged = np.where(act_map[0].data == act_map[0].header[act_type])
                plt.scatter(xc[wflagged], yc[wflagged], marker=badmarker, color=color)
            ax._autoscaleXon, ax._autoscaleYon = autoscale_state


class GPIWoofer(DeformableMirror):
    def __init__(self):
        DeformableMirror.__init__(self, shape=(GPI_Globals.gpi_woof_n, GPI_Globals.gpi_woof_n))
        self.name = "GPI Woofer"
        self.pupil_diam = 8.6   # for display, projected full area around 48x48 subaps of tweeter
        self.numacross = GPI_Globals.gpi_numacross
        self.actuator_spacing = GPI_Globals.gpi_woof_spacing
        self.pupil_center = GPI_Globals.pupil_center_woofer

    def annotate(self, marker='s', color='teal', s=50, alpha=0.4, **kwargs):
        """ Annotate the DM actuator coordinates.
        Applies some cosmetic defaults to distinguish Woofer from Tweeter actuators
        """
        DeformableMirror.annotate(self, marker=marker, color=color, s=s, alpha=alpha, **kwargs)


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

