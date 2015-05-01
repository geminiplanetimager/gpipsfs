import poppy
import numpy as np

class MEMS_DM(poppy.AnalyticOpticalElement):
    def __init__(self, size=(48,48), printthrough=True):
        super(MEMS_DM,self).__init__(name='MEMS DM')
        self.planetype=poppy.poppy_core._PUPIL
        self.size = size
        self.commands = np.zeros(size)

        self.actuator_spacing=0.18 # meters/actuator projected onto the primary

    def set_shape(self, new_commands):
        raise NotImplementedError('Not implemented yet!')

    def getPhasor(self,wave):
        """ DM surface """

        # GPI tweeter actuators are reimaged to 18 cm subapertures

        # Boston DM print through info in:
        #  http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf
        #  ao4elt.lesia.obspm.fr/sites/ao4elt/IMG/ppt/Bifano.ppt

        # in horizontal direction, the print through is about 35/190 pixels = 18% of the width
        # in the vertical direction, closer to 31%, but it's more like 2 narrow bars each 10% wide
        # and there's a 10% wide dot in the middle of it too
        #printthrough properties:
        pt_col_width = 0.18
        pt_col_value = 15e-9  # 15 nm, estimated from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.513.4649&rep=rep1&type=pdf
        pt_row_width = 0.10
        pt_row_value = -15e-9

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

        #raise NotImplementedError('Not implemented yet!')

