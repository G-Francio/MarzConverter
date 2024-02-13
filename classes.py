import getpass
import warnings

import numpy as np


class QUser:
    """
    A class representing a user with optional database access credentials.

    Attributes:
        uname (str): The username for database access.
        pwd (str): The password for database access.
        has_db_access (bool): A flag indicating whether the user has database access.

    Methods:
        set_credentials(): Set the database access credentials.
        get_credentials(): Get the database access credentials.
        set_db_access(): Enable database access for the user.
        get_db_access(): Check if the user has database access.
    """

    def __init__(self, uname=None, pwd=None, has_db_access=False):
        """
        Initialize a QUser object.

        Parameters:
            uname (str): The username for database access.
            pwd (str): The password for database access.
            has_db_access (bool): A flag indicating whether the user has database access.
        """
        self.has_db_access = has_db_access
        self.uname = None if not has_db_access else uname
        self.pwd = None if not has_db_access else pwd

    def set_credentials(self):
        """
        Set the database access credentials.
        """
        if self.has_db_access:
            self.uname = input("Username: ")
            self.pwd = getpass.getpass("Password: ")
        else:
            warnings.warn(
                """User does not appear to have access to the QUBIRCS database.\n
            If you do have access, please start the script with the optional argument:
            `qaccess=True` to login and pull data from the DB."""
            )

    def get_credentials(self):
        """
        Get the database access credentials.

        Returns:
            tuple: A tuple containing the username and password.
        """
        return self.uname, self.pwd

    def set_db_access(self):
        """
        Enable database access for the user.
        """
        self.has_db_access = True

    def get_db_access(self):
        """
        Check if the user has database access.

        Returns:
            bool: True if the user has database access, False otherwise.
        """
        return self.has_db_access


class SpecCollection:
    """
    A class representing a collection of spectra.

    Attributes:
        collection (list): A list containing individual spectra.
        x (numpy.ndarray): Array of x-axis values for all spectra.
        y (numpy.ndarray): Array of y-axis values for all spectra.
        dy (numpy.ndarray): Array of error values for all spectra.

    Methods:
        __iter__(): Iterate over the spectra in the collection.
        append(spec): Add a spectrum to the collection.
        get_names(): Get the names of all spectra in the collection.
        import_collection(lst): Import a list of spectra into the collection.
        pad_array(max_length=0): Pad arrays in each spectrum to a specified maximum length.
        pad_single_attr(spec, attr, max_length): Pad a single attribute array in a spectrum.
        complete_wave(): Complete the wavelength array in each spectrum with non-zero values.
        set_arrays(): Set x, y, and dy arrays for the entire collection.
    """

    def __init__(self):
        """
        Initialize a SpecCollection object.
        """
        self.collection = []
        self.x = None
        self.y = None
        self.dy = None

    def __iter__(self):
        """
        Iterate over the spectra in the collection.

        Yields:
            Spectrum: A spectrum from the collection.
        """
        for spec in self.collection:
            yield spec

    def append(self, spec):
        """
        Add a spectrum to the collection.

        Parameters:
            spec: The spectrum to be added to the collection.
        """
        self.collection.append(spec)

    def get_names(self):
        """
        Get the names of all spectra in the collection.

        Returns:
            list: A list containing the names of all spectra in the collection.
        """
        return [spec.name for spec in self.collection]

    def import_collection(self, lst):
        """
        Import a list of spectra into the collection.

        Parameters:
            lst (list): A list of spectra to be imported into the collection.
        """
        self.collection = lst

    def pad_array(self, max_length=0):
        """
        Pad arrays in each spectrum to a specified maximum length.

        Parameters:
            max_length (int): The maximum length to pad the arrays to.
        """
        _attrs = ["x", "y", "dy"]
        for spec in self.collection:
            for _attr in _attrs:
                self.pad_single_attr(spec, _attr, max_length)

    def pad_single_attr(self, spec, attr, max_length):
        """
        Pad a single attribute array in a spectrum to a specified maximum length.

        Parameters:
            spec: The spectrum containing the attribute array to be padded.
            attr (str): The attribute name (e.g., 'x', 'y', 'dy').
            max_length (int): The maximum length to pad the attribute array to.
        """
        current = getattr(spec, attr)
        padded_array = np.zeros(max_length) * current.unit
        padded_array[: current.shape[0]] = current
        setattr(spec, attr, padded_array)

    def complete_wave(self):
        """
        Complete the wavelength array in each spectrum with non-zero values.
        """
        for spec in self.collection:
            is_zero = np.where(spec.x.value == 0)[0]
            if len(is_zero) == 0:
                continue

            max_wave = np.max(spec.x.value[np.nonzero(spec.x.value)])

            top_wave = np.linspace(1.01 * max_wave, 1.1 * max_wave, len(is_zero))
            spec.x[is_zero] = top_wave * spec.x.unit

    def set_arrays(self):
        """
        Set x, y, and dy arrays for the entire collection.
        """
        self.x = [spec.x for spec in self.collection]
        self.y = [spec.y for spec in self.collection]
        self.dy = [spec.dy for spec in self.collection]
