# -------------------------------- ** -------------------------------- #
# Actual FITS parsing: retrieves information based on the instrument.  #
# Same structure, sometimes cards have different names thus different  #
# functions.                                                           #
# Works provided an instrument is given, otherwhise it can't set the   #
# correct cards.                                                       #
# -------------------------------- ** -------------------------------- #

import numpy as np


def parseWFCCD(hdul):
    """
    Parses information from a given HDU, for data produced at WFCCD
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CD1_1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(hdul[0].data[0], (1, hdul[0].data[0].shape[1]))
    error = np.reshape(hdul[0].data[1], (1, hdul[0].data[1].shape[1]))

    return (wave, flux, error)


def parseLDSS3(hdul):
    """
    Parses information from a given HDU, for data produced at WFCCD
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CDELT1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * 0.1

    return (wave, flux, error)


def parseIMACS(hdul):
    """
    Parses information from a given HDU, for data produced at IMACS
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CDELT1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * 0.1

    return (wave, flux, error)


def parseFIRE(hdul):
    """
    Parses information from a given HDU, for data produced at FIRE
    """
    data = hdul[5].data

    wave = data.field("WAVE")
    flux = data.field("FLUX")
    error = data.field("SIG")

    return (wave, flux, error)


def parseEFOSC(hdul):
    """
    Parses information from a given HDU, for data produced at EFOSC
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CDELT1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * 0.1

    return (wave, flux, error)


def parseLRS(hdul):
    """
    Parses information from a given HDU, for data produced at TNG LRS
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CDELT1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    r_wav = np.argwhere((wave >= 3700) & (wave <= 8000))  # reduced_wave,
    # TNG spectra are very noisy at the extremes of the wavelength range

    wave = np.reshape(wave[r_wav][:, 0], (1, wave[r_wav][:, 0].shape[0]))
    flux = np.reshape(
        hdul[0].data[r_wav][:, 0], (1, hdul[0].data[r_wav][:, 0].shape[0])
    )
    error = flux * 0.1

    return (wave, flux, error)


def parseGeneric(hdul):
    """
    Parses information from a generic HDU. Will fail most of the time,
    for every fail I will try to improve the function. This handles
    calibrated Gaia spectra at the minimum.
    """
    wave = hdul[1].data["wave"].reshape(1, -1)
    flux = hdul[1].data["flux"].reshape(1, -1)
    err = hdul[1].data["err"].reshape(1, -1)

    return wave, flux, err


def parseSDSS(hdul):
    """
    Parses information from SDSS spectra.
    """

    def revIVar(x, m):
        if x == 0:
            return m
        return np.sqrt(1 / x)

    vectRevIVar = np.vectorize(revIVar)

    data = np.array([np.array(i) for i in hdul[1].data])

    flux = data[:, 0].reshape(1, -1)
    wave = (10 ** data[:, 1]).reshape(1, -1)
    error = vectRevIVar(data[:, 2], max(flux)).reshape(1, -1)

    return (wave, flux, error)


def parseLAMOST(hdul):
    """
    Parses information from LAMOST spectra.
    """

    def revIVar(x, m):
        if x == 0 or x < 0:
            return m
        return np.sqrt(1 / x)

    vectRevIVar = np.vectorize(revIVar)

    data = np.array([np.array(i) for i in hdul[0].data])

    flux = data[0, :].reshape(1, -1)
    wave = data[2, :].reshape(1, -1)
    error = vectRevIVar(data[1, :], max(flux)).reshape(1, -1)

    return (wave, flux, error)


def parseAstrocook(hdul):
    """
    Parses information from spectra produced by Astrocook
    """
    data = hdul[1].data

    wave = data.x * 10  # Marz wants A, not nm
    flux = data.y
    error = flux * 0.01

    nanwave = np.isfinite(wave)
    flux[np.isnan(flux)] = np.nanmedian(flux)
    error[np.isnan(error)] = np.inf

    return (
        wave[nanwave].reshape(1, -1),
        flux[nanwave].reshape(1, -1),
        error[nanwave].reshape(1, -1),
    )


def parse6DF(hdul):
    """
    Parses information from spectra downloaded by 6dfGS
    """
    data = hdul[7].data

    wave = data[3]
    flux = data[0]
    error = flux * 0.1
    return (wave.reshape(1, -1), flux.reshape(1, -1), error.reshape(1, -1))


def parse2DF(hdul):
    """
    Parses information from spectra downloaded by 6dfGS
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CD1_1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    # Transform flux, should not matter for redshift identification but might
    #  as well consider it
    BZERO = hdul[0].header["BZERO"]
    BSCALE = hdul[0].header["BSCALE"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    flux = BSCALE * hdul[0].data + BZERO

    error = hdul[2].data

    return (wave.reshape(1, -1), flux.reshape(1, -1), error.reshape(1, -1))


def parseGoodman(hdul):
    """
    Parses information from a given HDU, for data produced at Goodman (Cerro Tololo)
    """
    start = hdul[0].header["CRVAL1"]
    step = hdul[0].header["CDELT1"]
    total = hdul[0].header["NAXIS1"]
    corr = hdul[0].header["CRPIX1"]

    wave = (np.arange(1, total + 1) - corr) * step + start
    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * 0.1

    return (wave, flux, error)
