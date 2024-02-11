import os
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
from astropy.io import fits
from fg_utils.utilities import is_number

import parse


def generate_fibres_data(db_data):
    """
    Given data from QDB, produces fibre data.
    If no data are found on the DB (or the DB can't be accessed) neutral, mock
    data are generated on the fly (e.g. ra/dec = 0/0, `z_spec = -`)

    Parameters:
        db_data (list): List containing data from QubricsDB.

    Returns:
        fits.BinTableHDU: Fibres extension header data unit.
    """
    name, t, ra, dec, comm = [], [], [], [], []
    t = ["P"] * len(db_data)

    for data in db_data:
        z_name = float(data[4]) if is_number(data[4]) else -1
        name.append(str(data[0]) + " - " + str(round(z_name, 2)))
        ra.append(str(float(data[1]) * np.pi / 180))
        dec.append(str(float(data[2]) * np.pi / 180))
        comm.append(generate_comment(data))

    name_col = fits.Column(name="NAME", format="80A", array=name)
    type_col = fits.Column(name="TYPE", format="1A", array=t)
    ra_col = fits.Column(name="RA", format="1D", array=ra)
    dec_col = fits.Column(name="DEC", format="1D", array=dec)
    comm_col = fits.Column(name="COMMENT", format="80A", array=comm)

    out_cols = fits.ColDefs([name_col, type_col, ra_col, dec_col, comm_col])
    return fits.BinTableHDU().from_columns(out_cols, name="fibres")


def generate_comment(db_data):
    """
    Generates the comment string in the Fibres Extension.

    Parameters:
        db_data (list): List containing data from QubricsDB.

    Returns:
        str: Generated comment string.
    """
    t = db_data[3]
    z = db_data[4]
    tf = "P" if db_data[5] == "" else db_data[5]
    qf = db_data[6]
    n = db_data[7]
    return str(t) + " " + str(z) + " " + tf + qf + " - " + n


def read_spec_list(path):
    """
    Reads a list of spectra to process.

    Parameters:
        path (str): Path to the file containing the list of spectra.

    Returns:
        list: List of spectra to convert, each item containing a file path and optional wavelength range.
    """
    with open(path) as f:
        read_data = [line.strip("\n") for line in f.readlines()]

    spec_to_convert = []
    for data in read_data:
        split_args = data.split("wr")
        spec = Path(split_args[0].strip(' "'))
        # Hardcodes .fits/.fit as the only acceptable format
        if spec.suffix not in [".fits", ".fit"]:
            continue

        wr = parse.parse_wave_range(split_args[1]) if len(split_args) > 1 else None
        spec_to_convert.append([spec, wr])

    return spec_to_convert


def gen_name(argd, verbose=False):
    """
    Generate input and output file names based on the provided arguments.

    Parameters:
        argd (dict): Parsed command-line arguments.
        verbose (bool): If True, print output file path.

    Returns:
        tuple: Tuple containing input FITS file path, output FITS file path, and file path stem.
    """
    fits_in = Path(argd["infile"])
    path = fits_in.parent / fits_in.stem

    # Sets the correct outfile:
    name = (
        Path(argd["outfile"])
        if argd["outfile"] is not None
        else path.parent / f"{path.name}_Marz.fits"
    )

    if verbose:
        print(f"Saving to: {name}")

    # Check that we are indeed saving a FITS file
    if name.suffix != ".fits":
        name = name.parent / f"{name.name}.fits"

    return fits_in, name, path


def fits_to_gspec(fits_in, name, wave_range=None):
    """
    Reads a FITS file and calls the appropriate function.
    Returns `wave`, `flux`, and `error`.
    If error is not found, the original FITS, error is assumed .1
    of the original flux.

    Parameters:
        fits_in (str): Path to the input FITS file.
        name (str): Name of the instrument.
        wave_range (list): Wavelength range for region extraction.

    Returns:
        parse.gspec: An instance of the generic_spectrum class.
    """
    with fits.open(fits_in) as hduList:
        spec = parse.parse_data(hduList, name)

    if wave_range is not None:
        return spec.region_extract(*wave_range, in_place=True)

    return spec


def _add_to_list(spec):
    """
    Helper function to convert a spectrum to generic_spectrum.

    Parameters:
        spec (list): List containing a file path and optional wavelength range.

    Returns:
        parse.gspec: An instance of the generic_spectrum class.
    """
    return fits_to_gspec(spec[0], spec[0].stem, wave_range=spec[1])


@contextmanager
def tempinput(data, suffix=".mzc"):
    """
    Context manager for creating a temporary input file.

    Parameters:
        data (str): Data to be written to the temporary file.
        suffix (str): Suffix for the temporary file.

    Yields:
        str: Temporary file name.
    """
    temp = NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    temp.write(data)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)
