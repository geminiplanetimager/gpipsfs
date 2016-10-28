import gpipsfs

def test_coron():
    gpi = gpipsfs.GPI()
    gpi.obsmode='H_coron'
    psf = gpi.calc_psf(monochromatic=1.6e-6)

    assert psf[0].data.sum() < 5e-4

def test_direct():
    gpi = gpipsfs.GPI()
    gpi.obsmode='H_direct'
    psf = gpi.calc_psf(monochromatic=1.6e-6)

    assert psf[0].data.sum() > 0.99

def test_unblocked():
    gpi = gpipsfs.GPI()
    gpi.obsmode='H_unblocked'
    psf = gpi.calc_psf(monochromatic=1.6e-6)
    assert psf[0].data.sum() > 0.35
    assert psf[0].data.sum() < 0.40

def test_obsmode():

    def check_modes(gpi, apod, occ, lyot, filt):
        assert gpi.apodizer == apod, 'Got unexpected apodizer value. Was {}, expected {}'.format(gpi.apodizer, apod)
        assert gpi.occulter == occ, 'Got unexpected occulter value. Was {}, expected {}'.format(gpi.occulter, occ)
        assert gpi.lyotmask == lyot, 'Got unexpected lyotmask value. Was {}, expected {}'.format(gpi.lyotmask, lyot)
        assert gpi.filter == filt, 'Got unexpected filter value. Was {}, expected {}'.format(gpi.filter, filt)




    gpi = gpipsfs.GPI()
    gpi.obsmode='H_direct'
    check_modes(gpi, 'CLEAR','SCIENCE','Open', 'H')

    gpi.obsmode='H_coron'
    check_modes(gpi, 'H','H','080m12_04', 'H')


    gpi.obsmode='K1_unblocked'
    check_modes(gpi, 'K1','SCIENCE','080m12_06_03', 'K1')

    gpi.obsmode='NRM_J'
    check_modes(gpi, 'NRM','SCIENCE','Open', 'J')


