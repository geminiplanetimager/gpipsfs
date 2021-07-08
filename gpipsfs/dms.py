import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import poppy

from .main import GeminiPrimary
my_path = os.path.abspath(os.path.dirname(__file__))
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
        self._surface = np.zeros(shape)      # array for the DM surface WFE
        self.numacross = shape[0]           # number of actuators across diameter of
                                            # the optic's cleared aperture (may be
                                            # less than full diameter of array)
        #self.actuator_spacing = 1.0/self.numacross  # distance between actuators,
                                                    # projected onto the primary
        self.pupil_center = (shape[0]-1.)/2 # center of clear aperture in actuator units
                                            # (may be offset from center of DM)
    
    @property
    def actuator_spacing(self):
        return  self.pupil_diam/self.numacross
    @property
    def shape(self):
        return self._shape

    @property
    def surface(self):
        """ The surface shape of the deformable mirror, in
        **meters** """
        return self._surface

    def set_surface(self, new_surface, units='nm'):
        """ Set the entire surface shape of the DM.

        Parameters
        -------------
        new_surface : 2d ndarray
            Desired DM surface shape
            (note that wavefront error will be 2x this)
        units : string
            Right now this *must* be 'nm' for nanometers,
            which is the default. Other units may be added later
            if needed.

        """

        assert new_surface.shape == self.shape
        if units!='nm':
            raise NotImplementedError("Units other than nanometers not yet implemented.")
        self._surface[:] = np.asarray(new_surface, dtype=float)*1e-9


    def set_actuator(self, actx, acty, new_value, units='nm'):
        """ Set the entire surface shape of the DM.

        Parameters
        -------------
        actx, acty : integers
            Coordinates of the actuator you wish to control
        new_value : float
            Desired surface height for that actuator
            (note that wavefront error will be 2x this)
        units : string
            Right now this *must* be 'nm' for nanometers,
            which is the default. Other units may be added later
            if needed.

        Example
        -----------
        dm.set_actuator(12,22, 123.4)

        """


        # FIXME do something more comprehensive with units
        assert units=='nm'
        if actx < 0 or actx > self.shape[1]-1:
            raise ValueError("X axis coordinate is out of range")
        if acty < 0 or acty > self.shape[0]-1:
            raise ValueError("Y axis coordinate is out of range")
        self._surface[acty, actx] = new_value*1e-9


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


    def get_opd(self,wave):
        """ Return the surface optical path delay for the optic.
        Interpolates from the current optic surface state onto the
        desired coordinates for the wave.

        CAUTION: This right now uses a fairly simple representation
        of the actuator influence function, which should not be
        taken too seriously just yet.
        """

        # the following could be replaced with a higher fidelity model if needed
        interpolated_surface = self._get_surface_via_gaussian_influence_functions(wave)
        return interpolated_surface


        #phasor = np.exp(1.j * 2 * np.pi * interpolated_surface/wave.wavelength)
        #return phasor

    def _get_surface_via_gaussian_influence_functions(self, wave):
        """ Infer a finely-sampled surface from simple Gaussian influence functions centered on
        each actuator.

        Work in progress, oversimplified, not a great representation of the true influence function
        """
        y, x = wave.coordinates()
        y_act, x_act = self.get_coordinates(one_d=True)

        interpolated_surface = np.zeros(wave.shape)

        crosstalk = 0.15 # amount of crosstalk on advancent actuator
        sigma = self.actuator_spacing/np.sqrt((-np.log(crosstalk)))

        pixelscale = x[0,1]-x[0,0]          # scale of x,y
        boxsize =  (3*sigma)/pixelscale     # half size for subarray (THIS IS NOT USED? --douglase)

        for yi, yc in enumerate(y_act):
            for xi, xc in enumerate(x_act):
                if self._surface[yi,xi] == 0: continue

                # 2d Gaussian
                r = ((x - xc)**2 + (y-yc)**2)/sigma**2

                interpolated_surface +=  np.nan_to_num(self._surface[yi,xi]) * np.exp(-r) #*np.sqrt(np.pi)

        return interpolated_surface


    def display(self, annotate=False, grid=False, what='opd', crosshairs=False, *args, **kwargs):
        """Display an Analytic optic by first computing it onto a grid.

        Parameters
        ----------
        wavelength : float
            Wavelength to evaluate this optic's properties at
        npix : int
            Number of pixels to use when sampling the optical element.
        what : str
            What to display: 'intensity', 'surface' or 'phase', or 'both'
        ax : matplotlib.Axes instance
            Axes to display into
        nrows, row : integers
            # of rows and row index for subplot display
        crosshairs : bool
            Display crosshairs indicating the center?
        colorbar : bool
            Show colorbar?
        colorbar_orientation : bool
            Desired orientation, horizontal or vertical?
            Default is horizontal if only 1 row of plots, else vertical
        opd_vmax : float
            Max value for OPD image display, in meters.
        title : string
            Plot label
        """


        if what=='both': raise NotImplementedError('still need to implement display both mode for display_actuators')

        kwargs['crosshairs']= crosshairs
        kwargs['what'] = what
        returnvalue = poppy.AnalyticOpticalElement.display(self, *args, **kwargs)

        if annotate: self.annotate()
        if grid: self.annotate_grid()
        return returnvalue


    def display_actuators(self, annotate=False, grid=True,  what='surface', crosshairs=False,  *args, **kwargs):
        """ Display the optical surface, viewed as discrete actuators

        Parameters
        ------------

        annotate : bool
            Annotate coordinates and types of actuators on the display? Default false.
        grid : bool
            Annotate grid of actuators on the display?  Default true.
        what : string
            What to display: 'intensity' transmission, 'surface' or 'phase', or 'both'
        """

        # display in DM coordinates
        # temporarily set attributes appropriately as if this were a regular OpticalElement
        self.amplitude = np.ones_like(self.surface)
        self.opd = self.surface
        self.pixelscale = self.actuator_spacing


        # back compatibility for older poppy syntax (which is confusing)
        if what=='surface': what='phase'

        #then call parent class display
        returnvalue = poppy.OpticalElement.display(self,  what=what, crosshairs=crosshairs, **kwargs)

        # now un-set all the temporary attributes, since this is analytic and
        # these are unneeded
        del self.pixelscale
        del self.opd
        del self.amplitude

        if annotate: self.annotate()
        if grid: self.annotate_grid()
        return returnvalue


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

    def annotate_grid(self, linestyle=":", color="black", **kwargs):

        y_act, x_act = self.get_coordinates(one_d=True)

        ax = plt.gca()
        for x in x_act:
            plt.axvline(x+ (self.actuator_spacing/2), linestyle=linestyle, color=color)
        for y in y_act:
            plt.axhline(y+ (self.actuator_spacing/2), linestyle=linestyle, color=color)

class ContinuousMEMS(DeformableMirror):
    def __init__(self,
                 name="ContinuousMEMS",
                 pupil_diam=8.6,
                 numacross=None,
                 mems_print_through_m=15e-9,#amplitude, 15 nm, estimated from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf
                 actuator_map_file= os.path.join(my_path, 'data','GPI_tweeter_actuators.fits'),
                 **kwargs
                 ):
        DeformableMirror.__init__(self, **kwargs)
        if mems_print_through_m != 0.0:
            self.mems_print_through = True
            self._mems_print_through_amplitude = mems_print_through_m
        else:
            self.mems_print_through = False
        self.name = name
        self.pupil_diam = pupil_diam   # for display, projected full area around 48x48 subaps of tweeter
        if numacross is None:
            self.numacross = self.shape[0]
        else:
            self.numacross = numacross

        self._actuator_type_info = fits.open(actuator_map_file)

    @property
    def bad_actuators(self):
        """Returns a list of coordinate indices for the actuators which are
        nonoperable """
        act_map = self._actuator_type_info
        wflagged = np.where(  ( act_map[0].data == act_map[0].header['DEAD']) |
                              ( act_map[0].data == act_map[0].header['WEAK']) )

        output = []
        for i in range(len(wflagged[0])):
            yc,xc = wflagged[0][i], wflagged[1][i]

            label = 'DEAD' if (act_map[0].data[yc,xc] == act_map[0].header['DEAD'] ) else 'WEAK'
            output.append([xc,yc,label])

        return output

    def get_opd(self, wave):
        opd = DeformableMirror.get_opd(self,wave)

        if self.mems_print_through:
            mems_print_through_opd = self._get_opd_MEMS_print_through(wave)
            opd += mems_print_through_opd
        return opd

    def _get_opd_MEMS_print_through(self,wave):
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

        opd = np.zeros(wave.shape)
        y, x = wave.coordinates()


        pixscale = x[0,1] - x[0,0]
        opd[np.mod(x,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_row_value
        opd[np.mod(y,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_col_value

        return opd
        #phasor = np.exp(1.j * 2 * np.pi * opd/wave.wavelength)
        #return phasor

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



class GPITweeter(DeformableMirror):
    def __init__(self, mems_print_through=True):
        DeformableMirror.__init__(self, shape=(GPI_Globals.gpi_tweet_n, GPI_Globals.gpi_tweet_n))
        self.name = "GPI Tweeter"
        self.numacross = GPI_Globals.gpi_numacross
        #self.actuator_spacing = GPI_Globals.gpi_tweet_spacing
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
        wflagged = np.where(  ( act_map[0].data == act_map[0].header['DEAD']) |
                              ( act_map[0].data == act_map[0].header['WEAK']) )

        output = []
        for i in range(len(wflagged[0])):
            yc,xc = wflagged[0][i], wflagged[1][i]

            label = 'DEAD' if (act_map[0].data[yc,xc] == act_map[0].header['DEAD'] ) else 'WEAK'
            output.append([xc,yc,label])

        return output

    def get_opd(self, wave):
        opd = DeformableMirror.get_opd(self,wave)

        if self.mems_print_through:
            mems_print_through_opd = self._get_opd_MEMS_print_through(wave)
            opd += mems_print_through_opd
        return opd

    def _get_opd_MEMS_print_through(self,wave):
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

        opd = np.zeros(wave.shape)
        y, x = wave.coordinates()


        pixscale = x[0,1] - x[0,0]
        opd[np.mod(x,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_row_value
        opd[np.mod(y,self.actuator_spacing) <= (self.actuator_spacing*pt_col_width)] += pt_col_value

        return opd
        #phasor = np.exp(1.j * 2 * np.pi * opd/wave.wavelength)
        #return phasor

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
        #self.actuator_spacing = GPI_Globals.gpi_woof_spacing
        self.pupil_center = GPI_Globals.pupil_center_woofer

    def annotate(self, marker='s', color='teal', s=50, alpha=0.4, **kwargs):
        """ Annotate the DM actuator coordinates.
        Applies some cosmetic defaults to distinguish Woofer from Tweeter actuators
        """
        DeformableMirror.annotate(self, marker=marker, color=color, s=s, alpha=alpha, **kwargs)


