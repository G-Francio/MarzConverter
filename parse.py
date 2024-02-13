import re
from pathlib import Path

import fg_utils.format as fmt_spec
from fg_utils.spec import generic_spectrum as gspec

import clogger as cl
import defaults as d

# from argparse import args


def parse_ext_arguments(args, user):
    """
    Parses optional arguments for better handling and less issues with filenames and extensions.

    Parameters:
        args (list): List of command-line arguments.
        user: An object representing the user.

    Returns:
        dict: A dictionary containing parsed optional arguments.
    """
    if len(args) == 1:
        raise cl.NoInputError("Please provide a `fits` to work on.")

    outd = {"infile": args[1]}
    outd.update(d.DEFAULT_ARGS_DICT)  # Update with defaults to keep things tidy

    allowed_keys = list(d.DEFAULT_ARGS_DICT.keys())

    if len(args) > 2:
        for arg in args[2:]:
            if Path(arg).is_file() and (arg.endswith("fits") or arg.endswith("fit")):
                outd["infile"] += f"\n{arg}"
            else:
                k, v = arg.replace(" ", "").split("=")
                if k not in allowed_keys:
                    raise ValueError(
                        "Unknown optional argument. Valid options: wr [None], outfile [None], pooling [T/F], npool [NCPU].\nCall as, e.g., 'MarzConverter infile npool=6'."
                    )
                else:
                    outd[k] = v

    if outd["qaccess"]:
        user.set_db_access()
        user.set_credentials()

    return outd


def parse_data(hdul, name) -> gspec:
    """
    Parses a FITS file given the instrument name.
    Very likely to fail if the fits is non-standard, please report the log files
    if errors are encountered, possibly with the problematic header/spectrum.

    Parameters:
        hdul: Header Data Unit List representing the FITS file.
        name (str): The name of the instrument.

    Returns:
        gspec: An instance of the generic_spectrum class.
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
        raise cl.FitsParseError(
            "Can't parse fits file, please report log file (/tmp/MarzConverter.log)"
        )


def parse_wave_range(str):
    """
    Parses the wavelength range.

    Parameters:
        str (str): The input string containing the wavelength range information.

    Returns:
        list: A list containing the parsed wavelength range values as floats.
    """
    wr = str.strip('wr=[ ]"')
    return [float(i.strip(" ")) for i in wr.split(",")]
