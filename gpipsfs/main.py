import poppy
import numpy as np
import astropy.io.fits as fits
import astropy.units as u


_grayscale_pixels=4

def _grayscale_getphasor_decorator(getphasor):
    """ Decorator to add automatic resampling and greyscale pixels to
    any optical element's getPhasor function. 

    Works by calculating a more highly sampled phasor then binning down."""
    def decorated_getphasor(optic,wave):
        #newwave_shape = (wave.shape[0]*_grayscale_pixels, wave.shape[1]*_grayscale_pixels)

        # TODO FIXME extend this to work for both image and pupil planes?
        oversampled_wave = poppy.Wavefront( wavelength=wave.wavelength,
            npix=wave.shape[0]*_grayscale_pixels,
            diam=wave.diam,
            oversample=wave.oversample)

        oversampled_phasor = optic.getPhasor(oversampled_wave)

        rebinned_phasor = poppy.utils.rebin_array(oversampled_wave, (_grayscale_pixels,_grayscale_pixels))
        return rebinned_phasor






class GPI(poppy.Instrument):
    """ Class implementing GPI perfect PSF calculation """


    apodizer_list = ['Y', 'J', 'H', 'K1', 'K2', 'CLEAR','CLEAR_GPI','NRM','H_LIWA','H_STAR']
    occulter_list = ['Y', 'J', 'H', 'K1', 'SCIENCE']
    lyotmask_list = [ '080m12_02', '080m12_03', '080m12_04', '080m12_04_c', '080m12_06', 
                    '080m12_06_03', '080m12_07', '080m12_10', 'Open', 'Blank' ]

    # Table of 'modename': (filter, apodizer, occulter, lyotmask)
    # Consistent with TLC CONFIG.InstSeqAll file
    _obsmode_table = {  # OSBMODE      FILTER APODIZER  OCCULTER   LYOT MASK
                        'Y_coron':      ('Y',  'Y',     'Y',       '080m12_03'),
                        'J_coron':      ('J',  'J',     'J',       '080m12_04'),
                        'H_coron':      ('H',  'H',     'H',       '080m12_04'),
                        'K1_coron':     ('K1', 'K1',    'K1',      '080m12_06_03'),
                        'K2_coron':     ('K2', 'K2',    'K1',      '080m12_07'),
                        'H_starcor':    ('H',  'Hstar', 'H',       '080m12_03'),
                        'H_LIWAcor':    ('H',  'HLIWA', 'K1',      '080m12_04'),
                        'Y_unblocked':  ('Y',  'Y',     'SCIENCE', '080m12_03'),
                        'J_unblocked':  ('J',  'J',     'SCIENCE', '080m12_04'),
                        'H_unblocked':  ('H',  'H',     'SCIENCE', '080m12_04'),
                        'K1_unblocked': ('K1', 'K1',    'SCIENCE', '080m12_06_03'),
                        'K2_unblocked': ('K2', 'K2',    'SCIENCE', '080m12_07'),
                        'Y_direct':     ('Y',  'CLEAR', 'SCIENCE', 'Open'),
                        'J_direct':     ('J',  'CLEAR', 'SCIENCE', 'Open'),
                        'H_direct':     ('H',  'CLEAR', 'SCIENCE', 'Open'),
                        'K1_direct':    ('K1', 'CLEAR', 'SCIENCE', 'Open'),
                        'K2_direct':    ('K2', 'CLEAR', 'SCIENCE', 'Open'),
                        'NRM_Y':        ('Y',  'NRM',   'SCIENCE', 'Open'),
                        'NRM_J':        ('J',  'NRM',   'SCIENCE', 'Open'),
                        'NRM_H':        ('H',  'NRM',   'SCIENCE', 'Open'),
                        'NRM_K1':       ('K1', 'NRM',   'SCIENCE', 'Open'),
                        'NRM_K2':       ('K2', 'NRM',   'SCIENCE', 'Open'),
                        'DARK':         ('H',  'H',     'H',       'Blank')}


    pixelscale=0.0141  # IFS pixelscale
    """IFS lenslet scale, in arcseconds"""

    dms=False          # Include DMs?
    """ Include Woofer and Tweeter surfaces in model? """


    def __init__(self, obsmode='H_coron', npix=1024, lyot_tabs=True, satspots=True, undersized_secondary=False, 
            display_before=False, dms=False):
        """  GPI PSF Simulator

        Parameters
        --------------
        obsmode: string
            GPI Obsmode
        lyot_tabs : Bool
            Block out the bad actuators with tabs in the Lyot masks? Default is True
        satspots : Bool
            Include the sat spots grid? Default is True
        display_before : Bool
            Display the wavefront before the FPM and Lyot when using calcPSF(display=True)?
            This is mostly for pedagogical purposes. Default is False.
        npix : int
            Number of pixels to use across the pupil. Default is 1024.
        """

        super(GPI,self).__init__(name='GPI')
        self.lyot_tabs=lyot_tabs
        self.satspots=satspots
        self._display_before=display_before
        self.dms=dms
        self._undersized_secondary=undersized_secondary # use the (erroneous) value we used in GPI design
        self.obsmode = obsmode
        self.options['output_mode'] = 'detector'  # by default, always bin down to GPI actual pixels
        self.npix=npix


        from . import dms
        self.tweeter = dms.GPITweeter()
        self.woofer  = dms.GPIWoofer()

    def _getFilterList(self):
        """ Return filter list and dummy placeholder for synphot bandpasses, 
        in the manner expected by the poppy.Instrument class"""

        filter_list =   ['Y', 'J', 'H', 'K1', 'K2']
        return filter_list, None

    def _getSynphotBandpass(self,filtername):
        """ Return a pysynphot.ObsBandpass object for the given desired band.
        For GPI this piggybacks off the filter copes inside the GPI DRP directory
        """
        import os
        import astropy.table
        import pysynphot
        try:
            drpdir = os.environ['GPI_DRP_DIR']
        except:
            raise RuntimeError("Cannot find GPI DRP dir from $GPI_DRP_DIR environment variable")

        fn = os.path.join(drpdir, 'config', 'filters', 'GPI-filter-'+filtername+'.fits')
        filterdata = astropy.table.Table.read(fn)
        band = pysynphot.spectrum.ArraySpectralElement(
                throughput = filterdata['TRANSMISSION'].flatten(),
                wave=filterdata['WAVELENGTH'].flatten()*1e4, # convert from microns to angstroms
                waveunits='angstrom', name="GPI "+filtername)
        return band



    # create properties with error checking
    @property
    def apodizer(self):
        """Currently selected apodizer name"""
        return self._apodizer
    @apodizer.setter
    def apodizer(self, value):
        value = value.upper() # force to uppercase
        if value not in self.apodizer_list:
            raise ValueError("Instrument {0} doesn't have a apodizer called {1}.".format(self.name, value))
        self._apodizer = value

    @property
    def occulter(self):
        """Currently selected occulter name"""
        return self._occulter
    @occulter.setter
    def occulter(self, value):
        value = value.upper() # force to uppercase
        if value not in self.occulter_list:
            raise ValueError("Instrument {0} doesn't have a occulter called {1}.".format(self.name, value))
        self._occulter = value

    @property
    def lyotmask(self):
        """Currently selected lyotmask name"""
        return self._lyotmask
    @lyotmask.setter
    def lyotmask(self, value):
        # preserve case for this one since we're used to that with the lyot mask names
        if value not in self.lyotmask_list:
            raise ValueError("Instrument {0} doesn't have a Lyot mask called {1}.".format(self.name, value))
        self._lyotmask = value


    @property
    def obsmode_list(self):
        "Available Observation Modes"
        keys = self._obsmode_table.keys()
        keys.sort()
        return keys

    # Obsmode works differently since it's a meta-property that affects the other ones:
    @property
    def obsmode(self):
        """Currently selected obsmode name"""
        for modename, settings in self._obsmode_table.items():
            if (self.filter==settings[0].upper() and self.apodizer==settings[1].upper() and 
                self.occulter==settings[2].upper() and self.lyotmask==settings[3]):
                return modename
        return 'Custom'
    @obsmode.setter
    def obsmode(self, value):
        if value not in self.obsmode_list:
            raise ValueError("Instrument {0} doesn't have an Obsmode called {1}.".format(self.name, value))

        settings = self._obsmode_table[value]
        self.filter=settings[0]
        self.apodizer=settings[1]
        self.occulter=settings[2]
        self.lyotmask=settings[3]


    def _getDefaultFOV(self):
        """ return default FOV in arcseconds """
        return 2.8


    def calcPSF(self, wavelength=None, oversample=2, contrast_relative_to=None, verbose=False, *args, **kwargs):
        """ Compute a PSF

        Has all the same options as poppy.Instrument.calcPSF plus a few more

        Parameters
        ----------------
        nlambda : int
            How many wavelengths to model for broadband? 
            The default depends on how wide the filter is: (5,3,1) for types (W,M,N) respectively
        monochromatic : float, optional
            Setting this to a wavelength value (in meters) will compute a monochromatic PSF at that 
            wavelength, overriding filter and nlambda settings.
        wavelength : float, optional
            Synonym for `monochromatic`. 

        """

        # Wrapper to slightly apply cosmetic adjustment to output plots

        if verbose:
            original_return_intermediates=kwargs.get('return_intermediates',False)
            kwargs['return_intermediates'] = True


        if wavelength is not None and 'monochromatic' not in kwargs: kwargs['monochromatic'] = wavelength
        # we have adjusted above the default to 2, vs. the default=4 for poppy, so let's apply that.
        kwargs['oversample'] = oversample

        # The actual calculation:
        retval = super(GPI,self).calcPSF(*args, **kwargs)

        # Output and display
        if verbose:
            psf, intermediates = retval
            for plane in intermediates:
                print("{0:40s} total intensity = {1:.3g}".format(plane.location+",", plane.totalIntensity))
            if not original_return_intermediates: retval = psf

        if kwargs.get('display', False):
            import matplotlib.pyplot as plt
            plt.gcf().suptitle("GPI, "+self.obsmode, size='xx-large')
        return retval

    def _getOpticalSystem(self,fft_oversample=2, detector_oversample = None, fov_arcsec=2, fov_pixels=None, options=dict()):

        """ Return an OpticalSystem instance corresponding to the instrument as currently configured.

        When creating such an OpticalSystem, you must specify the parameters needed to define the
        desired sampling, specifically the oversampling and field of view.


        Parameters
        ----------

        fft_oversample : int
            Oversampling factor for intermediate plane calculations. Default is 2
        detector_oversample: int, optional
            By default the detector oversampling is equal to the intermediate calculation oversampling.
            If you wish to use a different value for the detector, set this parameter.
            Note that if you just want images at detector pixel resolution you will achieve higher fidelity
            by still using some oversampling (i.e. *not* setting `oversample_detector=1`) and instead rebinning
            down the oversampled data.
        fov_pixels : float
            Field of view in pixels. Overrides fov_arcsec if both set.
        fov_arcsec : float
            Field of view, in arcseconds. Default is 2
        options : dict
            Other arbitrary options for optical system creation


        Returns
        -------
        osys : poppy.OpticalSystem
            an optical system instance representing the desired configuration.

        """


        if detector_oversample is None: detector_oversample = fft_oversample

        #poppy_core._log.debug("Oversample: %d  %d " % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(name=self.name, oversample=fft_oversample)
        if 'source_offset_r' in options.keys(): optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options.keys(): optsys.source_offset_theta = options['source_offset_theta']

        optsys.npix = self.npix

        #---- set pupil intensity
        pupil_optic=GeminiPrimary(undersized=self._undersized_secondary)
        #if self._undersized_secondary:
            #pupil_optic.obscuration_diameter = 1.02375 # SM outer diameter (vs inner hole projected diameter)

        #---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            full_opd_path = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], basestring): # tuple with filename and slice
            full_opd_path =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else 
                    os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif isinstance(self.pupilopd, fits.HDUList): # OPD supplied as FITS HDUList object
            full_opd_path = self.pupilopd # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None:
            full_opd_path = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:"+str(type(self.pupilopd)))

        #---- apply pupil intensity and OPD to the optical model
        optsys.addPupil(name='Gemini Primary', optic=pupil_optic, opd=full_opd_path, opdunits='micron', rotation=self._rotation)

        if self.dms:
            optsys.addPupil(optic=self.woofer)
            optsys.addPupil(optic=self.tweeter)


        # GPI Apodizer
        apod = GPI_Apodizer(name=self.apodizer, satspots=self.satspots)
        optsys.addPupil(optic=apod)

        if self._display_before: optsys.addImage(optic=poppy.ScalarTransmission(name='Before FPM', transmission=1))

        # GPI FPM
        fpm = GPI_FPM(name=self.occulter)
        optsys.addImage(optic=fpm)

        if self._display_before: optsys.addPupil(optic=poppy.ScalarTransmission(name='Before Lyot', transmission=1))

        # GPI Lyot Mask
        lyot = GPI_LyotMask(name=self.lyotmask, tabs=self.lyot_tabs)
        optsys.addPupil(optic=lyot)

        #--- add the detector element.
        if fov_pixels is None:
            fov_pixels = np.round(fov_arcsec/self.pixelscale)
            if 'parity' in self.options.keys():
                if self.options['parity'].lower() == 'odd'  and np.remainder(fov_pixels,2)==0: fov_pixels +=1
                if self.options['parity'].lower() == 'even' and np.remainder(fov_pixels,2)==1: fov_pixels +=1

        optsys.addDetector(self.pixelscale, fov_pixels = fov_pixels, oversample = detector_oversample, name=self.name+" lenslet array")


        return optsys




class GeminiPrimary(poppy.CompoundAnalyticOptic):
    primary_diameter = 7.701      # projected diameter of baffle on M2
    obscuration_diameter = 1.2968 # projected diameter of M2 inner oversized hole
    support_angles = [90-43.10, 90+43.10, 270-43.10, 270+43.10]
    support_widths = [0.014,    0.01,     0.01,      0.01]   # laser vane is slightly thicker
    support_offset_y = [0.2179, -0.2179,  -0.2179,   0.2179]
    def __init__(self,  name='Gemini South Primary', undersized=False):
        outer = poppy.CircularAperture(radius=self.primary_diameter/2)
        outer.pupil_diam = 8.0*u.meter   # slightly oversized array

        # Secondary obscuration from pupil diagram provided by Gemini

        sr = self.obscuration_diameter/2
        if undersized:
            sr = 1.02375/2 # SM outer diameter (vs inner hole projected diameter)

        # FIXME could antialias using same technique as used for apodizer grids
        obscuration = poppy.AsymmetricSecondaryObscuration(
                            secondary_radius=sr,
                            support_angle=self.support_angles,
                            support_width=self.support_widths,
                            support_offset_y=self.support_offset_y)

        return super(GeminiPrimary,self).__init__(opticslist=[outer,obscuration], name=name)



    def display(self, npix=512, wavelength=1e-6, annotate=False, annotate_angle=45.0, crosshairs=False, **kwargs):
        ax = super(GeminiPrimary,self).display(npix=npix, wavelength=wavelength, crosshairs=crosshairs, **kwargs)

        if annotate:
            # if both transmission and phase were displayed, annotate both
            # else just one
            offset_len=15

            if not hasattr(ax,'__iter__'): ax = (ax,)
            for axis in ax:

                ang = np.deg2rad(annotate_angle)
                for diam in (self.primary_diameter, self.obscuration_diameter):
                    axis.annotate( xy = (diam/2*np.cos(ang),
                        diam/2*np.sin(ang)),
                        s="{:.3f} m".format(diam), color='red',
                        textcoords="offset points", xytext=(offset_len*np.cos(ang), offset_len*np.sin(ang)),
                        arrowprops=dict(arrowstyle="->", color='red'))

                for angle in self.support_angles:
                    d = self.primary_diameter
                    ang= np.deg2rad(angle)
                    #print(ang, offset_len*np.cos(ang), offset_len*np.sin(ang))
                    offset_len=40
                    axis.annotate( xy = (d/2*np.cos(ang), d/2*np.sin(ang)),
                        s="$\\theta=$ {:.1f} deg".format(angle), color='red',
                        textcoords="offset points", xytext=(offset_len*np.cos(ang), offset_len*np.sin(ang)),
                        arrowprops=dict(arrowstyle="->", color='red'))


class GPI_FPM(poppy.CircularOcculter):
    """ GPI Focal Plane Masks """
    # FPM, and **diameter** in milliarcsec from 
    # GPI Coronagraph Final Mask Documentation
    # (in GPIMAIN/Instrument Docs/coronagraph)
    FPM_table ={'Y': 156.2, 
                'J': 184.7,
                'H': 246.7,
                'K1': 306.3,
                'SCIENCE': 0.0}

    def __init__(self, name='H'):
        try:
            radius = self.FPM_table[name.upper()]/2
        except:
            raise ValueError("No GPI FPM named "+name)

        poppy.CircularOcculter.__init__(self, name='GPI FPM '+name, radius=radius*1e-3)
        self._default_display_size=2.7 # arcsec FOV for display
        # FIXME could antialias using the same skimage code as used for the NRM mask?


def GPI_Apodizer(name='H', *args, **kwargs):
    """ Factory function to return some suitable apodizer object 
    
    Parameters
    ------------
    name : string 
        Name of apodizer, e.g. 'H' or 'NRM' or 'CLEAR_GPI'
    satspots : bool
        Include satspot grid? (only applies to some apodizers)
    
    
    """

    if name=='NRM':
        return GPI_NRM(name=name, *args, **kwargs)
    elif name=='CLEAR_GPI':
        return GeminiPrimary(undersized=True, name='GPI Apodizer CLEAR_GPI', *args, **kwargs)
    else:
        return GPI_Coronagraphic_Apodizer(name=name, *args, **kwargs)


class GPI_Coronagraphic_Apodizer(poppy.AnalyticOpticalElement):
    """ GPI Apodizer for APLC 

    Parameters based on 'APLC Final Design REV0.3.pdf' by Remi Soummer, 2008-12-03
    GPI KT folderID=950

    """

    # Apodizers are designed to a 1% undersized pupil for alignment tolerances.
    # full pupil maps to 11.790 mm, 1% undersized to 11.672 mm
    magnification = GeminiPrimary.primary_diameter/0.011672
    # apodizer radial profile filename, grid line spacing, grid line width
    # grid parameters are in microns at the apodizer plane
    _apodizer_table = {
            'CLEAR': (None, 0, 0),
            'Y': ('GPI_Y_56_10', 7.5, 555),
            'J': ('GPI_J_56_10', 7.5, 555),
            'H': ('GPI_H_56_10', 7.5, 585),
            'HLIWA': ('GPI_H_69_10', 7.5, 585),
            'Hstar': ('GPI_H_69_10', 7.5, 585), # to be updated
            'K1': ('GPI_K1_56_10', 7.5, 585),
            'K2': ('GPI_K2_56_10', 7.5, 585),
            'NRM': ('GPI_K2_56_10', 7.5, 585),  # to be updated
            }

    def __init__(self, name='H', satspots=True, *args, **kwargs):
        poppy.AnalyticOpticalElement.__init__(self,planetype=poppy.poppy_core._PUPIL, name='GPI Apodizer '+name)
        self.pupil_diam=8.0*u.meter # default size for display
        import os
        self._apodname = name
        self._apod_params = self._apodizer_table[name]
        self.satspots=satspots
        if self._apodname != 'CLEAR' and self._apodname != 'CLEAR_GPI':
            import astropy.io.ascii as asc
            from scipy.interpolate import interp1d
            # Consider revising the following to use pkg_resources, but that may not matter
            # http://stackoverflow.com/questions/5897666/how-do-i-use-data-in-package-data-from-source-code
            my_path=os.path.abspath(os.path.dirname(__file__))
            self._apod_profile = asc.read(os.path.join(my_path,"data", self._apod_params[0]+".txt"))
            # apodizer profiles are given in 10 micron steps on the apodizer
            self._apod_radius =  np.arange(len(self._apod_profile))*0.00001*self.magnification
            self._apod_interpolator = interp1d(self._apod_radius, self._apod_profile['transmission'],
                    bounds_error=False, fill_value=0.0)


    def plot1d(self):
        import matplotlib.pyplot as plt
        plt.plot(self._apod_radius, self._apod_profile['transmission'], label=self.name)
        plt.xlabel('Radial separation at primary [m]')
        plt.ylabel('E field transmission')

    def get_transmission(self, wave):
        """ Compute the transmission inside/outside of the obscuration

        """
        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._PUPIL)

        self.transmission = np.ones(wave.shape)
        if self._apodname=='CLEAR': return self.transmission

        y, x = wave.coordinates()
        r = np.sqrt(x ** 2 + y ** 2) 

        # Draw the radial component from the lookup table
        self.transmission = self._apod_interpolator(r) 

        # Now draw the (thin!) sat spot lines
        if self.satspots:

            width = self._apod_params[1]*1e-6*self.magnification
            offset = self._apod_params[2]*1e-6*self.magnification # perpendicular to line grid, so we need to correct by sqrt(2)
            xoffset = offset * np.sqrt(2)

            nlines=11 # on each side of the origin

            for n in np.arange(nlines*2+1)-nlines:
                pxscl = wave.pixelscale.to(u.meter/u.pixel).value
                self._draw_subpixel_antialiased_45deg_line(self.transmission, y, x, 1, xoffset*n, width, pxscl)
                self._draw_subpixel_antialiased_45deg_line(self.transmission, y, x,-1, xoffset*n, width, pxscl)


        return self.transmission

    def _draw_subpixel_antialiased_45deg_line(self, array,  y,x, sign, offset, width, pixelscale):
        """
        Draw a subpixel aliased line, sort of

        Parameters
        ------------
        offset : float
            Y intercept of the line,  as in B for "y = M*x + B"
        width : float
            width of the line in 

        """

        # Find all pixels that are plausibly near to the line

        # need to use 2x the width here to catch fractionally-blocked pixels
        wclose  = np.where(np.abs(y - sign*(x+offset)) <= width*2)
        xclose = x[wclose]
        yclose = y[wclose]


        # Now, for each of those pixels we are going to divide it up into NxN subpixel
        # samples in order to better measure the line widths. 
        subpix_n=3


        # yes, for loops. This was confusing to code any other way. Consider 
        weights = np.zeros_like(xclose)

        # if we divide pixels up into N chunks, the coordinates of the
        # **centers of** the # extreme left/right chunks will 
        # be +- (N-1)/(2*N)
        # e.g. for 3 chunks, they are centered at -2/3 and +2/3
        outer = (subpix_n-1.)/(2.*subpix_n)
        for dx in np.linspace(-outer,outer, subpix_n)*pixelscale:
            for dy in np.linspace(-outer,outer, subpix_n)*pixelscale:
                weights[(((yclose+dy)*sign >= (xclose+dx+offset) ) &
                         ((yclose+dy)*sign <  (xclose+dx+offset+width)))] +=1
        weights /= subpix_n**2


        array[wclose] *= (1-weights)
        #array[wclose2] *= (1-weights)



class GPI_LyotMask(poppy.AnalyticOpticalElement):
    # From 2014 August email by Anand re Lyot mask geometry:
    #   For the sake of completeness, I quote from Remi's doc on apodizer designs:
    #  "The Gemini pupil at the Lyot plane is 9.825mm (without undersizing).'
    #   I presume Remi refers here to the GSOD after cardboard baffling of M2,
    #   viz. 7770.1mm, re-imaged down to the Lyot plane
    #   One final wrinkle is hat the spec is defined with the part at ambient temperature.

    magnification = 7.7701/.009825   # meters at primary/meters at Lyot

    # Locations of bad actuator tabs:
    #   tab locations in millimeters on the physical Lyot mask
    #   offset from center
    #   most masks have 4 tabs, only the 080m12_04_c has the 5th one
    #   3rd coordinate defines type of linear extension 1:radial, 2:to spider, 0: none
    #  REFERENCE:  CAD drawings and 'mask design notes 2012.pptx' in GPIMAIN/Instrument docs/coronagraph masks 
    bad_actuator_tab_locations = [ (3.977, 1.932,  1),
                                   (-2.614, 3.750, 1),
                                   (2.614, -3.977, 1),
                                   (-1.932,-1.250, 2),
                                   (1.023, -1.136, 0)]  # coupled to its neighbor, not fully broken
            #Name           R_out  R_in   spider_width  ntabs   
    lyot_table = {
            'Blank':       (0,     0,     0,     0, 0 ),
            'Open':        (5.1,   0,     0,     0, 0 ),
            '080m12_02':   (4.786, 1.388, 0.204, 4, 0.328),     # **not installed in GPI**
            '080m12_03':   (4.786, 1.388, 0.304, 4, 0.378),     # used with Y
            '080m12_03_06':(4.786, 1.388, 0.304, 4, 0.528),     # **not installed in GPI**
            '080_04':      (4.786, 1.388, 0.404, 3, 0.428),     # Older generation, does not have 4th tab
            '080m12_04':   (4.786, 1.388, 0.404, 4, 0.428),     # used with J and H
            '080m12_04_c': (4.786, 1.388, 0.604, 5, 0.378),     # extra 5th tab
            '080m12_06':   (4.786, 1.388, 0.604, 4, 0.528),     # 
            '080m12_06_03':(4.786, 1.388, 0.604, 4, 0.528),     # exact duplicate of above mask?
            '080m12_07':   (4.686, 1.488, 0.604, 4, 0.528),     #
            '080m12_10':   (4.377, 1.772, 1.004, 4, 0.728) }
    # Which ones are actually installed in the instrument currently?
    # REFERENCE: From TLC $tlc/config/CONFIG.IfsAssembly, retrieved on 2015-04-05
	#	LYOT_MASK_POS   Blank         0       0.0             0.0              0011
	#	LYOT_MASK_POS   080m12_03     36      -0.073249377   -0.047069377      0002
	#	LYOT_MASK_POS   080m12_04     72     -0.10638653     0.029958820     0004
	#	LYOT_MASK_POS   080_04        108     -0.11726604    -0.012915230     000A
	#	LYOT_MASK_POS   080m12_06     144     -0.052745055    -0.046359157      0016
	#	LYOT_MASK_POS   080m12_04_c   180     -0.20085273     -0.084516352     000E
	#	LYOT_MASK_POS   080m12_06_03  216     -0.11979799     -0.034771766     001C
	#	LYOT_MASK_POS   080m12_07     252     -0.14302189     -0.024171524      001A
	#	LYOT_MASK_POS   080m12_10     288     -0.14280973     -0.020752296     0014
	#	LYOT_MASK_POS   Open          324      0.0             0.0             0008

    # From Oct 25 2012 email from Jeff Chilcote to gpi_iandt:
    #   
    #   Here are the calculated Lyot & Mems rotations (in deg):
    #   
    #   Mems:92.42
    #   
    #   
    #   Name            Rotation             Remainder V1 V2 Average
    #   080_04          92.35473125 -0.065268750000001
    #   080m12_03   92.57181875 0.151818750000004
    #   080m12_04   92.01261875 -0.40738125
    #   080m12_06   91.6767375  -0.7432625
    #   080m12_04_c 92.35705            -0.062949999999986
    #   080m12_06_03    92.92439375 0.504393749999977
    #   080m12_07   93.1073625  0.687362499999992
    #   080m12_10   92.274825   -0.145175000000009
    #           
    #   Median Rotation 92.355890625    -0.064109374999994
    #   Mean Rotation   92.4099421875   -0.010057812500003
    #   
    def __init__(self, name='080m12_04', tabs=True):


        super(GPI_LyotMask,self).__init__(planetype=poppy.poppy_core._PUPIL, name='GPI Lyot '+name)
        self.pupil_diam=8.0 # default size for display
        # Read in some geometry info from the primary mirror class
        self.support_angles=GeminiPrimary.support_angles
        self.support_offset_y = GeminiPrimary.support_offset_y

        #'name': (primary radius,  secondary radius,  spiderwidth, ntabs, tabwidth)
        # All sizes are in mm on the physical parts; needs to be converted by
        # magnification to get the size reimaged onto the primary
        # case insensitive search
        for k in self.lyot_table.keys():
            if k.upper()==name.upper():
                cname=k
                break
        else:
            raise ValueError('Unknown Lyot stop name: '+name)
        params = self.lyot_table[cname]

        # convert params from mm to meters, then magnify onto the primary
        self.outer_radius = params[0]*1e-3*self.magnification
        self.inner_radius = params[1]*1e-3*self.magnification
        self.support_width = params[2]*1e-3*self.magnification
        self.ntabs= params[3]
        self.tabradius= params[4]*1e-3*self.magnification

        if tabs==False: self.ntabs=0  # allow user to turn tabs off if desired

        self.wavefront_display_hint = 'intensity' # preferred display for wavefronts at this plane

    def getPhasor(self, wave):
        """ Compute the transmission inside/outside of the obscuration

        Based on poppy.AsymmetricSecondaryObscuration but with the bad actuator tabs added.
        """
        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._PUPIL)

        self.transmission = np.ones(wave.shape)

        y, x = wave.coordinates()

        y *= -1 # Flip Y coordinate convention to match 
                # Lyot bad actuator tabs to AOWFS display bad actuators
        r = np.sqrt(x ** 2 + y ** 2)  #* wave.pixelscale

        self.transmission[r < self.inner_radius] = 0
        self.transmission[r > self.outer_radius] = 0

        for angle_deg, offset_y in zip(self.support_angles,
                self.support_offset_y):
            angle = np.deg2rad(angle_deg + 90)  # 90 deg offset is to start from the +Y direction
            width = self.support_width

            # calculate rotated x' and y' coordinates after rotation by that angle.
            # and application of offset
            xp =  np.cos(angle) * (x) + np.sin(angle) * (y-offset_y)
            yp = -np.sin(angle) * (x) + np.cos(angle) * (y-offset_y)

            self.transmission[(xp > 0) & (np.abs(yp) < width / 2)] = 0

            # TODO check here for if there are no pixels marked because the spider is too thin.
            # In that case use a grey scale approximation

        for itab in range(self.ntabs):

            offset_x = self.bad_actuator_tab_locations[itab][0]*1e-3*self.magnification
            offset_y = self.bad_actuator_tab_locations[itab][1]*1e-3*self.magnification
            #print (offset_x, offset_y)

            xo = x-offset_x
            yo = y-offset_y
            r = np.sqrt( xo**2+yo**2)
            self.transmission[r < self.tabradius] = 0

            if self.bad_actuator_tab_locations[itab][2] == 1: 
                # Extend tab radially outwards
                angle = np.arctan2(offset_y, offset_x) #-(np.pi/2)
                xp =  np.cos(angle) * xo + np.sin(angle) * yo
                yp = -np.sin(angle) * xo + np.cos(angle) * yo
                self.transmission[(xp > 0) & (np.abs(yp) < self.tabradius )] = 0
            elif self.bad_actuator_tab_locations[itab][2] == 2:
                # special case the one that hangs off the secondary support

                angle=np.deg2rad(self.support_angles[1]+90)
                angle += np.pi/2
                xp =  np.cos(angle) * xo + np.sin(angle) * yo
                yp = -np.sin(angle) * xo + np.cos(angle) * yo
                self.transmission[(xp > 0) & (xp < 0.5) & (np.abs(yp) < self.tabradius )] = 0

        return self.transmission

    def fancy_display(self, npix=512, wavelength=1e-6, **kwargs):
        """ Display the Lyot OPD and also the primary mirror pupil itself on the same scale

        """
        phasor, pixelscale = self.sample(wavelength=wavelength, npix=npix, what='complex',
                                         return_scale=True)

        pri_phasor, pri_pixelscale = GeminiPrimary().sample(wavelength=wavelength, 
                    npix=npix, what='complex',return_scale=True)


        # temporarily set attributes appropriately as if this were a regular OpticalElement
        self.amplitude = np.abs(phasor+pri_phasor)
        phase = np.angle(phasor) / (2 * np.pi)
        self.opd = phase * wavelength
        self.pixelscale = pixelscale

        asdasdasd
        #then call parent class display
        #super(GPI_LyotMask,self) pop.display(**kwargs)
        poppy.OpticalElement.display(self, **kwargs)

        # now un-set all the temporary attributes back, since this is analytic and
        # these are unneeded
        #self.pixelscale = None
        #self.opd = None
        #self.amplitude = None


    # From Bruce email 2012 Sep 26:
    #
    #    I have a piece of code (attached) that generates gemini pupils and
    #    lyot masks. it's a bit clunky and slow (it generates an oversampled
    #    pupil and then downsizes it to get grey pixels write; it could be done
    #    much more elegantly with modern IDL.) Note that for historical reasons
    #    it uses units and coordinates of "milimeters in the lyot plane" even
    #    though some parameters have names like dtel_m. It doesn't
    #    automatically generate the MEMS tabs - that requires some manual work.
    #
    #    here's an example to generate our lyot 4 mask:
    #    npixdtel=43*8. ; lisa geometry for full GS pupil
    #    dtel_m=9.826 ; this corresponds to the full GS pupil, projected into
    #    mm on
    #                  ; the lyot plane
    #    id_m=2.775 ; oversized as per remi's table of lyot masks
    #    od_m=9.571;undersized as per Remi's table of lyot masks
    #    spider_m=0.04*od_m ; Lyot_4_0 spiders
    #    spoff_m=0.2804 ; offset from center
    #    pixsize = 7.770 / npixdtel ; pixel size in primary meters per pixel
    #    npixout = npixdtel
    #
    #    lyot4
    #    =
    #    make_lyot(npixout,npixdtel,dtel_m,id_m,od_m,spider_m,spoff_m,greysize=4)
    #
    #    and here's an example to generate the raw Gemini pupil
    #
    #      npixdtel=43*8. ; lisa geometry for full GS pupil
    #    dtel_m=9.826 ; this corresponds to the full GS pupil, projected into
    #    mm on
    #                  ; the lyot plane
    #    id_m=1.29464
    #    od_m=9.826 ; no undersizing
    #    spider_m=0.018 ; raw gemini spiders
    #    spoff_m=0.2804 ; offset from center
    #    pixsize = 7.770 / npixdtel ; pixel size in primary meters per pixel
    #    npixout = npixdtel
    #
    #    gempupil
    #    =
    #    make_lyot(npixout,npixdtel,dtel_m,id_m,od_m,spider_m,spoff_m,greysize=6)





# Names of masks
# TLC Name -> Gemini header name
#   Blank -> LYOT_BLANK_G6232
#   Open -> LYOT_OPEN_G6231</to>
#   080m12_03 -> 080m12_03
#   080m12_04 -> 080m12_04
#   080_04 -> 080_04
#   080m12_06 -> 080m12_06
#   080m12_04_c -> 080m12_04_c
#   080m12_06_03 -> 080m12_06_03
#   080m12_07 -> 080m12_07
#   080m12_10 -> 080m12_10


class GPI_NRM(poppy.AnalyticOpticalElement):

    # Mask geometry information copied from NRM_mask_definitions.py
    # provided by Alex Greenbaum
    hole_diam = 0.595806967827971
    centers = [[-0.68712086181030,  2.514046357726287],
              [-0.251922728788131,  3.362423670611769],
              [1.8223921820303386,  0.157370753458909],
              [2.3417804300716787, -0.644378188031338],
              [-2.861816294382397, -1.733021136856143],
              [-0.869101033505584, -3.287947799633272],
              [-3.025663210535089,  1.567878988164694],
              [2.6921408318053010,  1.854772995499249],
              [3.2970144274045676, -0.595806967827971],
              [1.0355384147357893, -3.192100591765294]]
    rotdeg = 115.0

    def __init__(self, *args, **kwargs):
        #super(GPI_NRM,self).__init__(planetype=poppy.poppy_core._PUPIL, name='GPI NRM mask')
        poppy.AnalyticOpticalElement.__init__(self,planetype=poppy.poppy_core._PUPIL, name='GPI NRM mask')
        self.pupil_diam=8.0 # default size for display


        try:
            import skimage
            self._antialias=True
        except ImportError:
            self._antialias=False


    def _drawCircle(self, array, xind, yind, xcen, ycen, radius, value=1):
        """ Draw a circle in an array, of a given radius centered at (ycen,xcen)
        where the coordinates of the array indixes are (yind, xind)
        """

        if self._antialias:
            # use anti-aliased circle
            #  see http://scikit-image.org/docs/dev/auto_examples/plot_shapes.html#example-plot-shapes-py
            from skimage.draw import circle_perimeter_aa
            # assuming xind and yind are monotonic linear series,
            # convert to pixel coords
            pxscl = xind[0,1]-xind[0,0]
            xpix = int(round((xcen-xind[0,0])/pxscl))
            ypix = int(round((ycen-yind[0,0])/pxscl))
            rad =  int(round(radius/pxscl))
            rr,cc,val = circle_perimeter_aa(ypix,xpix,rad)
            array[rr,cc] = val
            # FIXME the above is imperfect because the method in skimage only works
            # on integer coordinates. TBD to see if we care enough to do anything better than
            # that.

        r = np.sqrt((xind-xcen) ** 2 + (yind-ycen) ** 2)
        array[r<radius]=value




    def getPhasor(self, wave):
        """ Compute the transmission inside/outside of the obscuration

        """
        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._PUPIL)

        self.transmission = np.zeros(wave.shape)
        y, x = wave.coordinates()

        c, s = np.cos(np.deg2rad(self.rotdeg)), np.sin(np.deg2rad(self.rotdeg))

        for coords in self.centers:
            xrot = c*coords[0] - s*coords[1]
            yrot = s*coords[0] + c*coords[1]
            self._drawCircle(self.transmission, x, y, xrot, yrot, self.hole_diam/2)

        return self.transmission


# Here is a copy of the POPPY routine, added just so I can slightly change how the
# npix is defined.  Once that functionality is merged into poppy this will be unnecessary
def _patched_inputWavefront(self, wavelength=2e-6):
    """Create a Wavefront object suitable for sending through a given optical system, based on
    the size of the first optical plane, assumed to be a pupil.

    If the first optical element is an Analytic pupil (i.e. has no pixel scale) then
    an array of 1024x1024 will be created (not including oversampling).

    Uses self.source_offset to assign an off-axis tilt, if requested.

    Parameters
    ----------
    wavelength : float
        Wavelength in meters

    Returns
    -------
    wavefront : poppy.Wavefront instance
        A wavefront appropriate for passing through this optical system.

    """

    if hasattr(self, 'npix'):
        npix = int(self.npix)
    else:
        npix = self.planes[0].shape[0] if self.planes[0].shape is not None else 1024
    diam = self.planes[0].pupil_diam if hasattr(self.planes[0], 'pupil_diam') else 8

    inwave = poppy.Wavefront(wavelength=wavelength,
            npix = npix,
            diam = diam,
            oversample=self.oversample)
    poppy.poppy_core._log.debug("Creating input wavefront with wavelength=%f, npix=%d, pixel scale=%f meters/pixel" % (wavelength, npix, diam/npix))

    if np.abs(self.source_offset_r) > 0:
        offset_x = self.source_offset_r *-np.sin(self.source_offset_theta*np.pi/180)  # convert to offset X,Y in arcsec
        offset_y = self.source_offset_r * np.cos(self.source_offset_theta*np.pi/180)  # using the usual astronomical angle convention
        inwave.tilt(Xangle=offset_x, Yangle=offset_y)
        poppy.poppy_core._log.debug("Tilted input wavefront by theta_X=%f, theta_Y=%f arcsec" % (offset_x, offset_y))
    return inwave
        # Monkey patch the class to allow defining a rescalable pupil based on this class's npix attribute
        # FIXME this is a hack and should be removed once this functionality is merged into poppy master


# Monkey patch the class to allow defining a rescalable pupil based on this class's npix attribute
# FIXME this is a hack and should be removed once this functionality is merged into poppy master

poppy.OpticalSystem.inputWavefront = _patched_inputWavefront


