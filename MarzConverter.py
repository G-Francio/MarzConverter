#!/usr/bin/python3
import sys
from multiprocessing import Pool

from astropy.io import fits
from tqdm import tqdm

import classes as c
import defaults as d
import misc
import parse

# -------------------------------- ** -------------------------------- #
# #################################################################### #
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
    Convert data for Marz.

    Calls the appropriate method depending on arguments.
    Accepts a single FITS file or a FITS list in the format of an input file or as a sequence passed
    to the script. The output is always a single file containing all required information.
    The first argument is mandatory and is the input file. It can be either a FITS file, a file list
    compiled in a text file (with extension `.txt`, `.dat` or `.mzc`) or multiple `.fits` files passed
    as input one after the other (see examples below).
    Both absolute and relative paths should work. Additional optional arguments are allowed.

    Args:
        **kwargs: Additional keyword arguments. These can be:
        - "wr" (optional): The wavelength range. Should be formatted as `wr = [start, end]` with no spaces.
        - "outfile" (optional): Path to the output file. If not provided, the default is
          set to the user's home directory with the name "Converted_Marz.fits".
        - "pooling" (optional): A boolean indicating whether pooling is enabled or not. Default is False.
        - "npool" (optional): Number of pooling processes. If not provided, the default is half the available CPUs.
        - "qaccess" (optional): A boolean indicating whether access to the QUBRICS is granted. Default is False.

    Returns:
        Spec or SpecCollection: Processed data for Marz. Note: SpecCollection is simply a collection of Spec.

    Examples:
        MarzConverter input.fits
        MarzConverter input_list.txt
        MarzConverter input1.fits input2.fits input3.fits --optional_arg value
    """
    argd = parse.parse_ext_arguments(kwargs["sysargs"], mdb.user)

    if (
        (argd["infile"].endswith(".txt"))
        or argd["infile"].endswith(".dat")
        or argd["infile"].endswith(".mzc")
    ):
        return multi_fits_to_file(argd)
    elif "\n" in argd["infile"]:
        with misc.tempinput(argd["infile"]) as tempfilename:
            argd["infile"] = tempfilename
            if not argd["outfile"]:
                argd["outfile"] = d.HOME / "Converted_Marz.fits"
            return multi_fits_to_file(argd)
    else:
        return fits_to_file(argd)


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def fits_to_file(argd):
    """
    Process a single FITS file.

    Reads a FITS file and calls the appropriate function.
    Writes the resulting FITS to file, ready for MARZ.
    If an error is not found, the original FITS is assumed. The error is assumed to be 0.1 times the original flux.

    Args:
        argd (dict): Parsed command-line arguments.

    Returns:
        c.Spec: Processed data for Marz.
    """
    fits_in, name, path = misc.gen_name(argd, verbose=True)

    spec_db_data = mdb.get_observation_data([path.name])
    fibre_hdu = misc.generate_fibres_data(spec_db_data)

    with fits.open(fits_in) as hdul:
        spec = parse.parse_data(hdul, name)

    wave_range = parse.parse_wave_range(argd["wr"]) if argd["wr"] is not None else None
    if wave_range is not None:
        spec.region_extract(*wave_range, in_place=True)

    write_fits(spec, name, fibre_hdu)
    return spec


# -------------------------------- ** -------------------------------- #


def multi_fits_to_file(argd):
    """
    Process a list of FITS files.

    Reads a file list and calls the appropriate function for each FITS file.

    Args:
        argd (dict): Parsed command-line arguments.

    Returns:
        c.SpecCollection: Processed data for Marz.
    """
    fits_in, name, _ = misc.gen_name(argd, verbose=True)

    spec_list = c.SpecCollection()
    spec_files = misc.read_spec_list(fits_in)

    if argd["pooling"]:
        pool = Pool(int(argd["npool"]))
        spec_list.import_collection(
            list(tqdm(pool.imap(misc._add_to_list, spec_files)))
        )
    else:
        spec_list.import_collection(
            [misc._add_to_list(spec_file) for spec_file in tqdm(spec_files)]
        )

    spec_db_data = mdb.get_observation_data(spec_list.get_names())
    fibre_hdu = misc.generate_fibres_data(spec_db_data)

    max_shape = max([s.x.shape[0] for s in spec_list])
    spec_list.pad_array(max_shape)
    spec_list.complete_wave()
    spec_list.set_arrays()

    write_fits(spec_list, name, fibre_hdu)
    return spec_list


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def write_fits(spec, name, fibre=None):
    """
    Write processed FITS file for Marz.

    Writes the processed FITS file, ready for Marz. Asks for overwrite permission if the file already exists.

    Args:
        spec (c.Spec or c.SpecCollection): Processed data for Marz.
        name (str): Name of the output FITS file.
        fibre (astropy.io.fits.ImageHDU, optional): Fibre data to be included in the FITS file.

    Returns:
        int: 0 if successful.
    """
    primaryHDU = fits.PrimaryHDU(spec.y)
    varianceHDU = fits.ImageHDU(spec.dy, name="variance")
    waveHDU = fits.ImageHDU(spec.x, name="wavelength")

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
    print(f"Saved to: {name}")
    return 0


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

if __name__ == "__main__":
    MarzConverter(sysargs=sys.argv)
