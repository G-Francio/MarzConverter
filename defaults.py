import os
from importlib import util
from pathlib import Path

# Check if mariadb is installed
HAS_MDB = True if util.find_spec("mariadb") is not None else False

# For multiprocessing purposes
NCPU = os.cpu_count()

# Access to the DB
HAS_QUBRICS_ACCESS = False

# Types in the QUBIRCS DB
DEFAULT_TYPE_DICT = {
    1: "QSO",
    2: "Type2",
    3: "BLLac",
    4: "Galaxy",
    5: "Star",
    6: "EmLines",
    7: "Uncertain",
    8: "Extended?",
    9: "closeneighb",
}

DEFAULT_ARGS_DICT = {
    "wr": None,
    "outfile": None,
    "pooling": False,
    "npool": None,
    "qaccess": False,
}

COMMON_FITS_HAS_ERR = ["WFCCD/WF4K-1"]
COMMON_FITS_NO_ERR = ["IMACS Short-Camera", "MagE", "EFOSC", "Goodman Spectr"]

HOME = Path.home()


# Functions
def typedict(key):
    """
    Retrieves the value associated with the given key from the DEFAULT_TYPE_DICT dictionary.

    Parameters:
    - key (str): The key for which the corresponding value will be retrieved.

    Returns:
    str: The value associated with the provided key in the DEFAULT_TYPE_DICT dictionary. If the key is not found, an empty string is returned.
    """
    try:
        return DEFAULT_TYPE_DICT[key]
    except KeyError:
        return ""


def get_default_data(nameList):
    """
    Generates default data if getting data fails.

    Parameters:
    - nameList (list): A list of names for which mock data will be generated.

    Returns:
    numpy.ndarray: An array containing mock observation data for each name in the input list.
    """
    observation_data_fallback = []
    for name in nameList:
        mockData = get_default_data_single(name)
        observation_data_fallback.append(mockData)
    return observation_data_fallback


def get_default_data_single(name):
    """
    Generates mock data for a single name if getting data fails.

    Parameters:
    - name (str): The name for which mock data will be generated.

    Returns:
    list: A list containing mock observation data for the given name.
    """
    return [name, "0", "0", "-", "-", "-", "-", "-"]
