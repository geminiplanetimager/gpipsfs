from .main import *
from .dms import GPITweeter, GPIWoofer

from poppy import display_psf

__version__ = '0.2.0'

# Imperfect display hack - set poppy display default
# for the intermediate planes to the GPI FOV size
import poppy
poppy.conf.default_image_display_fov=1.4 # (half size)
