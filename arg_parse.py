import argparse
import os
from pathlib import Path

home = Path.home()

workers = os.cpu_count() if os.cpu_count() is not None else 1
if "sched_getaffinity" in dir(os):
    workers = len(os.sched_getaffinity(0))

parser = argparse.ArgumentParser(
    prog="MarzConverter.py",
    description="Convert spectra in a Marz (https://samreay.github.io/Marz/) compatible format.",
)

parser.add_argument(
    dest="infile",
    metavar="input_files",
    type=str,
    nargs="+",
    help="Input fits(s) file(s).",
)

parser.add_argument(
    "--wr",
    dest="wr",
    metavar="WAVELENGTH RANGE",
    type=str,
    nargs="*",
    help="Wavelength range for the spectum. To be given in Angstrom, for example: --wr 4000 8000.",
)

parser.add_argument(
    "--outfile",
    dest="outfile",
    type=str,
    default=home / "Converted_Marz.fits",
    help="Sets the output file.",
)

parser.add_argument(
    "--pooling",
    dest="pooling",
    action="store_true",
    help="Allow or disallows multiprocess though pooling.",
)

parser.add_argument(
    "--npool",
    dest="npool",
    metavar="npool",
    action="store",
    default=workers // 2,
    help="Number of CPUs used for pooling. If not set but used with `--pooling`, then --npool = ncpu // 2",
)

parser.add_argument(
    "--qaccess",
    dest="qaccess",
    action="store_true",
    help="Indicates whether the user has access to the QUBRICS database.",
)

args = parser.parse_args()

# Additional optional arguments are allowed and can be specified in the kwargs dictionary:
#     - "wr" (optional): A placeholder; currently not implemented.
#     - "outfile" (optional): Path to the output file. If not provided, the default is
#       set to the user's home directory with the name "Converted_Marz.fits".
#     - "pooling" (optional): A boolean indicating whether pooling is enabled or not. Default is False.
#     - "npool" (optional): Number of pooling processes. If not provided, the default is None.
#     - "qaccess" (optional): A boolean indicating whether quick access is enabled or not. Default is False.
