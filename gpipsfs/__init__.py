from .main import *
from .dms import GPITweeter, GPIWoofer
from .utils import plot_contrast

from poppy import display_psf

__version__ = '0.3.0dev'

# Imperfect display hack - set poppy display default
# for the intermediate planes to the GPI FOV size
import poppy
poppy.conf.default_image_display_fov=2.8  # GPI IFS FOV size
