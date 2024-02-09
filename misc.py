from pathlib import Path

from astropy.io import fits
from fg_utils.utilities import is_number
from numpy import pi


def generate_fibres_data(db_data):
    """
    Given data from QDB, produces fibre data.
    If no data are found on the DB (or the DB can't be accessed) neutral, mock
    data are generated on the fly (e.g. ra/dec = 0/0, `z_spec = -`)
    """
    name, t, ra, dec, comm = [], [], [], [], []
    t = ["P"] * len(db_data)

    for data in db_data:
        z_name = float(data[4]) if is_number(data[4]) else -1
        name.append(str(data[0]) + " - " + str(round(z_name, 2)))
        ra.append(str(float(data[1]) * pi / 180))
        dec.append(str(float(data[2]) * pi / 180))
        comm.append(generate_comment(data))

    name_col = fits.Column(name="NAME", format="80A", array=name)
    type_col = fits.Column(name="TYPE", format="1A", array=t)
    ra_col = fits.Column(name="RA", format="1D", array=ra)
    dec_col = fits.Column(name="DEC", format="1D", array=dec)
    comm_col = fits.Column(name="COMMENT", format="80A", array=comm)

    out_cols = fits.ColDefs([name_col, type_col, ra_col, dec_col, comm_col])
    return fits.BinTableHDU().from_columns(out_cols, name="fibres")


# -------------------------------- ** -------------------------------- #


def generate_comment(db_data):
    """
    Generates the comment string in the Fibres Extension.
    """
    t = db_data[3]
    z = db_data[4]
    tf = "P" if db_data[5] == "" else db_data[5]
    qf = db_data[6]
    n = db_data[7]
    return str(t) + " " + str(z) + " " + tf + qf + " - " + n


# -------------------------------- ** -------------------------------- #


def parse_wave_range(str):
    """
    Parses the wavelength range.
    """
    wr = str.strip('wr=[ ]"')
    return [float(i.strip(" ")) for i in wr.split(",")]


# -------------------------------- ** -------------------------------- #


def read_spec_list(path):
    """
    Reads a list of spectra to process.
    """
    with open(path) as f:
        read_data = [line.strip("\n") for line in f.readlines()]

    spec_to_convert = []
    for data in read_data:
        split_args = data.split("wr")
        spec = Path(split_args[0].strip(' "'))
        # Hardcodes .fits/.fit as only acceptable format
        if spec.suffix not in [".fits", ".fit"]:
            continue

        wr = parse_wave_range(split_args[1]) if len(split_args) > 1 else None
        spec_to_convert.append([spec, wr])

    return spec_to_convert
