#!/usr/bin/python3
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.io import fits
from fg_utils import format as fmt_spec
from tqdm import tqdm

import classes as c
import clogger
import defaults as d
import misc

# Custom logging
sys.excepthook = clogger.handle_exception


home = str(Path.home())

# TODO: Add handling of external errors
# TODO: Check if fallback for no INSTRUME cards is enough

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

# Attempts to gather appropriate data from QubricsDB.
# If mariadb module is not installed, falls back on mock data.

# -------------------------------- ** -------------------------------- #

if d.HAS_MDB:
    import import_mdb as mdb
else:
    import import_no_mdb as mdb

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def MarzConverter(**kwargs):
    """
    MarzConverter! Calls the appropriate method depending on arguments.
    Can accept a single fits, or a fits list in the format of an input file.
    The output will always be a single file, containing all required informations.
    Can accept up to 3 arguments:
     - first argument is mandatory (input file); can be a fits file or a file list (.txt, .dat, .mzc)
       Absolute or relative path should both work.
     - second and third argument are optional:
        if a single fits is given:
        - second arg. can either be a wavelength range (`wr=[w1,w2]`), or the outfile.
          Extension is automatically added and is not required
        - third arg. can only be the outfile, and is used only if second arg is a `wr`
        if a file list is give:
        - no third argument is used, second argument is outfile

    If no `wr` or outfile are given the outfile will be called `infile_Marz.fits`;
    moreover, the initial `wr` will be used.
    """
    argd = parse_ext_arguments(kwargs["sysargs"], mdb.user)

    if (
        (argd["infile"].endswith(".txt"))
        or argd["infile"].endswith(".dat")
        or argd["infile"].endswith(".mzc")
    ):
        if argd["pooling"]:
            multi_fits_to_file_pooled(argd)
        else:
            multi_fits_to_file(argd)
    else:
        fits_to_file(argd)


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def parse_ext_arguments(args, user):
    """Parses optional arguments for better handling and less issues with filenames and extensions."""
    outd = {"infile": args[1]}
    outd.update(d.DEFAULT_ARGS_DICT)  # Update with defaults to keep things tidy

    allowed_keys = list(d.DEFAULT_ARGS_DICT.keys())

    if len(args) > 2:
        for arg in args[2:]:
            k, v = arg.replace(" ", "").split("=")
            if k not in allowed_keys:
                raise ValueError(
                    "Unknown optional argument. Valid options: wr [None], outfile [None], pooling [T/F], npool [NCPU].\nCall as, e.g., 'MarzConverter infile npool=6'."
                )
            else:
                outd[k] = v

    if outd["qaccess"]:
        user.set_credentials()

    return outd


def fits_to_file(argd):
    """
    Reads a FITS file and calls the appropriate function.
    Writes the resulting FITS to file, ready for MARZ.
    If error is not found the original FITS, error is assumed .1
    of the original flux.
    """
    fits_in = Path("/home/francio/test.fits")  # "argd["infile"])
    path = fits_in.parent / fits_in.stem

    # Sets the correct outfile:
    name = (
        Path(argd["outfile"])
        if argd["outfile"] is not None
        else path.parent / f"{path.name}_Marz.fits"
    )

    spec_db_data = mdb.get_observation_data([path.name])
    fibre_hdu = misc.generate_fibres_data(spec_db_data)

    with fits.open(fits_in) as hdul:
        spec = parse_data(hdul, name)

    wave_range = misc.parse_wave_range(argd["wr"]) if argd["wr"] is not None else None
    if wave_range is not None:
        spec.region_extract(*wave_range, in_place=True)

    write_fits(spec, fibre=fibre_hdu, name=name)


# -------------------------------- ** -------------------------------- #


def multi_fits_to_file(argd):
    """
    Reads a file list and calls the appropriate function for each fits.
    """
    fits_in = Path(argd["infile"])
    path = fits_in.parent / fits_in.stem

    name = (
        argd["outfile"] + ".fits"
        if argd["outfile"] is not None
        else path.parent / f"{path.name}_Marz.fits"
    )

    spec_list = c.spec_collection()
    spec_files = misc.read_spec_list(fits_in)

    for spec_file in tqdm(spec_files):
        spec_file_name = spec_file[0].stem
        spec = fits_to_array(spec_file_name[0], wave_range=spec_file_name[1])
        spec_list.append(spec)

    spec_db_data = mdb.get_observation_data(name_list)
    fibre_hdu = misc.generate_fibres_data(spec_db_data)

    max_shape = max([s.x.shape[1] for s in spec_list])
    spec_list.pad_array(max_length=max_shape)
    spec_list.complete_wave()

    write_fits(spec_list, fibre=fibre_hdu, name=name)


# -------------------------------- ** -------------------------------- #


def _addToList(spec):
    specFileName = p.splitext(spec[0])[0].split("/")[-1]
    wave, flux, error = fits2array(spec[0], waveRange=spec[1])
    return specFileName, wave, flux, error


def multiFits2FilePooled(argd):
    """
    Reads a file list and calls the appropriate function for each fits.
    """

    path, _ = p.splitext(argd["infile"])
    name = (
        argd["outfile"]
        if argd["outfile"] is not None
        else path.split("/")[-1] + "_Marz.fits"
    )

    # Fix the double extension issue
    if not name.endswith(".fits"):
        name = name + ".fits"

    waveList, fluxList, errorList, nameList = [], [], [], []
    specFiles = readSpecList(argd["infile"])

    pool = Pool(int(argd["npool"]))
    r = list(tqdm(pool.imap(_addToList, specFiles)))

    for _r in r:
        nameList.append(_r[0])
        waveList.append(_r[1])
        fluxList.append(_r[2])
        errorList.append(_r[3])

    specDBData = getObservationData(nameList)
    fibreHDU = generateFibresData(specDBData)

    maxShape = max([s.shape[1] for s in waveList])
    waveList, fluxList, errorList = padArray(
        waveList, fluxList, errorList, maxLength=maxShape
    )

    completeWave(waveList)

    writeFits(fluxList, errorList, waveList, fibre=fibreHDU, name=name)


# -------------------------------- ** -------------------------------- #


def fits2array(fitsIn, waveRange=None):
    """
    Reads a FITS file and calls the appropriate function.
    Returns `wave`, `flux` and `error`.
    If error is not found the original FITS, error is assumed .1
    of the original flux.
    """
    with fits.open(fitsIn) as hduList:
        wave, flux, error = parseData(hduList)

    if waveRange is not None:
        return cutWavelength(wave, flux, error, waveRange)

    return wave, flux, error


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def parse_data(hdul, name):
    """
    Parses a FITS file given the instrument name.
    Very likely to fail if the fits is non standard, please report the log files
    if errors are encountered, possibly with the problematic header/spectrum.
    """
    header = hdul[0].header
    try:
        inst = header["INSTRUME"]
    except KeyError:
        inst = None

    if inst is None:
        try:
            inst = header["TELESCOP"]
        except KeyError:
            inst = None

    try:
        origin = header["ORIGIN"]
    except KeyError:
        origin = None

    if origin is not None and origin == "Astrocook":
        return fmt_spec.parse_astrocook_fits(hdul, name)

    if inst is None:
        return fmt_spec.parse_generic(hdul, name)
    elif inst in d.COMMON_FITS_HAS_ERR:
        return fmt_spec.parse_common_fits(hdul, name, has_err=True)
    elif re.search("LDSS3-.*", inst) is not None or inst in d.COMMON_FITS_NO_ERR:
        return fmt_spec.parse_common_fits(hdul, name)
    elif inst == "FIRE":
        return fmt_spec.parse_fire(hdul, name)
    elif inst == "LRS":  # TNG
        return fmt_spec.parse_lrs(hdul, name)
    elif inst == "SDSS 2.5-M":
        return fmt_spec.parse_sdss(hdul, name)
    elif inst == "LAMOST":
        return fmt_spec.parse_lamost(hdul, name)
    elif inst == "SuperCOSMOS I" or inst == "2dF" or inst == "6dF":
        return fmt_spec.parse_2df_6df(hdul, name)
    else:
        print("Can't parse fits file, please report log file (/tmp/MarzConverter.log)")
        return None


# -------------------------------- ** -------------------------------- #


def padArray(*args, maxLength=0):
    """
    Allows array of different dimensions to be stacked together.
    Fills missing data with zero, centres the shortest array.
    """
    returnList = []
    for arg in args:
        innerReturnList = np.array([padSingleArray(a, maxLength) for a in arg])
        returnList.append(np.vstack(innerReturnList))
    return np.array(returnList)


# -------------------------------- ** -------------------------------- #


def padSingleArray(array, ml=0):
    """
    Performs the padding operation.
    """
    paddedArray = np.zeros((1, ml))
    paddedArray[: array.shape[0], : array.shape[1]] = array
    return paddedArray


# -------------------------------- ** -------------------------------- #


def completeWave(array):
    for wave in array:
        maxWave = np.max(wave[np.nonzero(wave)])
        iszero = np.argwhere(wave == 0).flatten()

        if len(iszero) == 0:
            continue

        topWave = np.linspace(1.01 * maxWave, 1.1 * maxWave, len(iszero))
        wave[iszero] = topWave


# -------------------------------- ** -------------------------------- #


def writeFits(flux, error, wave, fibre=None, name="MarzConverterOutput.fits"):
    """
    Writes the process fits file, ready for Marz. Asks for overwrite permission!
    """
    primaryHDU = fits.PrimaryHDU(flux)
    varianceHDU = fits.ImageHDU(error, name="variance")
    waveHDU = fits.ImageHDU(wave, name="wavelength")

    if fibre is None:
        hduListOut = fits.HDUList([primaryHDU, varianceHDU, waveHDU])
    else:
        hduListOut = fits.HDUList([primaryHDU, varianceHDU, waveHDU, fibre])

    try:
        hduListOut.writeto(name)
        hduListOut.close()
    except OSError:
        overwrite = input("File already exists, overwrite (Y/n)? ")
        if overwrite.lower() == "y" or overwrite == "":
            hduListOut.writeto(name, overwrite=True)
            hduListOut.close()
        else:
            hduListOut.close()
    print("Saved to: " + name)
    return 0


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

if __name__ == "__main__":
    MarzConverter(sysargs=sys.argv)
