#!/usr/bin/python3

# Logging, see here: https://stackoverflow.com/questions/3383865/how-to-log-error-to-file-and-not-fail-on-exception
# and here: https://stackoverflow.com/questions/8050775/using-pythons-logging-module-to-log-all-exceptions-and-errors
# and here: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python/16993115#16993115

import logging, sys, warnings

logging.basicConfig(
    filename="/tmp/MarzConverter.log",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(handler)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.critical(
        "Uncaught exception, consider reporting the log (/tmp/MarzConverter.log).",
        exc_info=(exc_type, exc_value, exc_traceback),
    )


sys.excepthook = handle_exception

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

from multiprocessing import Pool
from astropy.io import fits
from os import cpu_count
from pathlib import Path

home = str(Path.home())

import numpy as np
import os.path as p, re

from tqdm import tqdm

NCPU = (
    cpu_count()
)

USER = None
PWD  = None

# TODO: Add handling of external errors
# TODO: Check if fallback for no INSTRUME cards is enough

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

# Attempts to gather appropriate data from QubricsDB.
# If mariadb module is not installed, falls back on mock data.


def getFallbackData(nameList):
    """
    Generates mock data if everything fails.
    """
    observationDataFallback = []
    for name in nameList:
        mockData = getFallbackDataSingle(name)
        observationDataFallback.append(mockData)
    return np.array(observationDataFallback)


# -------------------------------- ** -------------------------------- #


def getFallbackDataSingle(name):
    return [name, "0", "0", "-", "-", "-", "-", "-"]


# -------------------------------- ** -------------------------------- #

try:
    import mariadb as mdb

    def getUserCredential(user = USER, pwd = PWD):
        """
        Read user credentials for QubricsDB.
        """
        if user is None or pwd is None:
            user = input("Username: ")
            pwd  = input("Password: ")
        return user, pwd

    USER, PWD = getUserCredential()

    # -------------------------------- ** -------------------------------- #

    def DBConnect(user, pw):
        """
        Connects to QubricsDB, returns a cursor and the connection object.
        """
        conn = mdb.connect(
            user=user, password=pw, host="127.0.0.1", port=3306, database="Qubrics"
        )
        cur = conn.cursor()
        return cur, conn

    # -------------------------------- ** -------------------------------- #

    def typedict(key):
        d = {
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

        try:
            return d[key]
        except KeyError:
            return ""

    # -------------------------------- ** -------------------------------- #

    def getObservationData(namelist):
        """
        Queries the DB for complementary information available per target.
        Will either find data based on the QID or, if observed by us, from the named spectrum.
        """
        _nameListQIDs = [e for e in namelist if e.isdecimal()]
        _nameListFiles = [e for e in namelist if not e.isdecimal()]

        data = _getDataFromQID(_nameListQIDs)
        if len(data) == 0:
            return _getObservationData(_nameListFiles)

        return data

    # -------------------------------- ** -------------------------------- #        

    def _getDataFromQID(namelist):
        """
        Queries the DB is the name is numerical (and thus a qid).
        This allows to download complementary informations even if the target is not from
        our observations but from literature.
        """
        if len(namelist) == 0:
            return []

        _cred = getUserCredential()
        cur, conn = DBConnect(_cred[0], _cred[1])

        qidList = []

        # Can't query all info at once, in the unlikely case that a numerical name is not in the DB.
        # Also, this will return wrong info in that case. Can't do anything about it, need to check
        #  beforehand!

        for _qid in namelist:
            cur.execute("SELECT qid, RAd, DECd, otypeid, z_spec FROM Qubrics.All_info WHERE qid = ?", (_qid,))
            queryRes = [*cur]
            # Check if I get exactly one result out of the query
            if len(queryRes) == 1:
                qidList.append(np.array([*queryRes[0]] + ["", "A", ""]))
            elif len(queryRes) == 0:
                qidList.append(getFallbackDataSingle(_qid))
            elif len(queryRes) > 1:
                warnings.warn("Too many results from query, this should never happen!\nUsing fallback data.")

        conn.close()

        return np.array(qidList)

    # -------------------------------- ** -------------------------------- #

    def _getObservationData(nameList):
        """
        Gathers data for observations that needs conversion.
        `nameList` is a list of fileNames taken from Drive/QSOCandidates.
        """
        _cred = getUserCredential()

        # index = []
        observationData = []
        fileNameList = []

        try:
            cur, conn = DBConnect(_cred[0], _cred[1])
            cur.execute("SELECT file FROM Qubrics.Observations")
            for c in cur:
                fileNameList.append(c)
            fileNameList = np.array(fileNameList)[:, 0]

            for name in nameList:
                foundIndices = np.where(np.char.find(fileNameList, "/" + name) > 0)
                if len(foundIndices[0]) > 1:
                    print("Multiple matches for {}, using fallback data.".format(name))
                    observationData.append(getFallbackDataSingle(name))
                elif len(foundIndices[0]) == 0:
                    print(
                        "Can't find spec data on DB for {}, fallback on mock data.".format(
                            name
                        )
                    )
                    observationData.append(getFallbackDataSingle(name))
                else:
                    fileName = fileNameList[foundIndices[0][0]]
                    cur.execute(
                        "SELECT target_qid, RAd, DECd, otypeid, z_spec, targetflag, qflag, notes FROM Qubrics.Observations WHERE file = '{}'".format(
                            fileName
                        )
                    )
                    for c in cur:
                        c_objtype = (
                            c[0],
                            c[1],
                            c[2],
                            typedict(c[3]),
                            c[4],
                            c[5],
                            c[6],
                            c[7],
                        )
                        observationData.append(c_objtype)

            conn.close()
            return np.array(observationData)
        except mdb.OperationalError:
            print("Connection to the DB failed, mock data will be used")
            print(
                "Connect to the DB before running the script for possible additional data."
            )
            return getFallbackData(nameList)           

except ModuleNotFoundError:
    print("MariaDB module not found, fallback on mock data.")

    def getObservationData(nameList):
        """
        Fallback method if mariadb is not found.
        Returns a list of neutral items in order to preserve the same
        script used with the DB itself.
        """
        return getFallbackData(nameList)


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

# This is called at the beginning of the script


def MarzConverter(**kwargs):
    """
    MarzConverter! Calls the appropriate method depending on arguments.
    Can accept a single fits, or a fits list. The output will always be a single file,
    containing all required informations.
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
    argd = parseExtArguments(kwargs["sysargs"])
    if (
        (argd["infile"].endswith(".txt"))
        or argd["infile"].endswith(".dat")
        or argd["infile"].endswith(".mzc")
    ):
        if argd["pooling"]:
            multiFits2FilePooled(argd)
        else:
            multiFits2File(argd)
    else:
        fits2File(argd)


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def parseExtArguments(args):
    """Parses optional arguments for better handling and less issues with filenames and extensions."""
    outd = {"infile": args[1], "wr": None, "outfile": None, 'pooling':False, "npool": NCPU}
    allowedKeys = list(outd.keys())
    if len(args) > 2:
        for arg in args[2:]:
            k, v = arg.replace(" ", "").split("=")
            if k not in allowedKeys:
                raise ValueError(
                    "Unknown optional argument. Valid options: wr [None], outfile [None], pooling [T/F], npool [NCPU].\nCall as, e.g., 'MarzConverter infile npool=6'."
                )
            outd[k] = v
    return outd


def fits2File(argd):
    """
    Reads a FITS file and calls the appropriate function.
    Writes the resulting FITS to file, ready for MARZ.
    If error is not found the original FITS, error is assumed .1
    of the original flux.
    """
    fitsIn = argd["infile"]

    path, _ = p.splitext(fitsIn)
    # Sets the correct outfile:
    name = (
        argd["outfile"]
        if argd["outfile"] is not None
        else path.split("/")[-1] + "_Marz.fits"
    )

    specDBData = getObservationData([path.split("/")[-1]])
    fibreHDU = generateFibresData(specDBData)

    with fits.open(fitsIn) as hduList:
        wave, flux, error = parseData(hduList)

    waveRange = parseWR(argd["wr"]) if argd["wr"] is not None else None
    if waveRange is not None:
        wave, flux, error = cutWavelength(wave, flux, error, waveRange)

    writeFits(flux, error, wave, fibre=fibreHDU, name=name)


# -------------------------------- ** -------------------------------- #


def multiFits2File(argd):
    """
    Reads a file list and calls the appropriate function for each fits.
    """
    path, _ = p.splitext(argd["infile"])
    name = (
        argd["outfile"] + ".fits"
        if argd["outfile"] is not None
        else path.split("/")[-1] + "_Marz.fits"
    )

    waveList, fluxList, errorList, nameList = [], [], [], []
    specFiles = readSpecList(argd["infile"])

    for spec in tqdm(specFiles):
        specFileName = p.splitext(spec[0])[0].split("/")[-1]
        wave, flux, error = fits2array(spec[0], waveRange=spec[1])
        nameList.append(specFileName)
        waveList.append(wave)
        fluxList.append(flux)
        errorList.append(error)

    specDBData = getObservationData(nameList)
    fibreHDU = generateFibresData(specDBData)

    maxShape = max([s.shape[1] for s in waveList])
    waveList, fluxList, errorList = padArray(
        waveList, fluxList, errorList, maxLength=maxShape
    )

    completeWave(waveList)

    writeFits(fluxList, errorList, waveList, fibre=fibreHDU, name=name)


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
        argd["outfile"] + ".fits"
        if argd["outfile"] is not None
        else path.split("/")[-1] + "_Marz.fits"
    )

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
    print(fitsIn)
    with fits.open(fitsIn) as hduList:
        wave, flux, error = parseData(hduList)

    if waveRange is not None:
        return cutWavelength(wave, flux, error, waveRange)

    return wave, flux, error


# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #


def parseData(hdul):
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
        return parseAstrocook(hdul)

    if inst is None:
        return parseGeneric(hdul)
    elif inst == "WFCCD/WF4K-1":
        return parseWFCCD(hdul)
    elif re.search("LDSS3-.*", inst) is not None or inst == "MagE":
        return parseLDSS3(hdul)
    elif inst == "IMACS Short-Camera":
        return parseIMACS(hdul)
    elif inst == "FIRE":
        return parseFIRE(hdul)
    elif inst == "EFOSC":
        return parseEFOSC(hdul)
    elif inst == "LRS":  # TNG
        return parseLRS(hdul)
    elif inst == "SDSS 2.5-M":
        return parseSDSS(hdul)
    elif inst == "LAMOST":
        return parseLAMOST(hdul)
    elif inst == "SuperCOSMOS I":
        return parse6DF(hdul)
    else:
        print("Can't parse fits file, please report log file (/tmp/MarzConverter.log)")


# -------------------------------- ** -------------------------------- #


def cutWavelength(wave, flux, error, waveRange):
    """
    Reduces a spectrum to a given wavelength range.
    """
    wr = waveRange
    reduced_wav = np.argwhere((wave >= wr[0]) & (wave <= wr[1]))
    wave = wave[:, reduced_wav[:, 1]]
    flux = flux[:, reduced_wav[:, 1]]
    error = error[:, reduced_wav[:, 1]]

    return wave, flux, error


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


def readSpecList(fileIn):
    """
    Reads a list of spectra to process.
    """
    with open(fileIn) as f:
        readData = [line.strip("\n") for line in f.readlines()]

    spec2Convert = []
    for data in readData:
        splitArgs = data.split("wr")
        spec = splitArgs[0].strip(' "')
        wr = parseWR(splitArgs[1]) if len(splitArgs) > 1 else None
        spec2Convert.append([spec, wr])

    return np.array(spec2Convert, dtype=object)


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


# -------------------------------- ** -------------------------------- #


def parseWR(str):
    """
    Parses the wavelength range.
    """
    wr = str.strip('wr=[ ]"')
    return [float(i.strip(" ")) for i in wr.split(",")]


# -------------------------------- ** -------------------------------- #

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def generateFibresData(DBData):
    """
    Given data from QDB, produces fibre data.
    If no data are found on the DB (or the DB can't be accessed) neutral, mock
    data are generated on the fly (e.g. ra/dec = 0/0, `z_spec = -`)
    """
    name, t, ra, dec, comm = [], [], [], [], []
    t = ["P"] * len(DBData)

    for data in DBData:
        z_name = float(data[4]) if is_number(data[4]) else -1
        name.append(str(data[0]) + ' - ' + str(round(z_name, 2)))
        ra.append(str(float(data[1]) * np.pi / 180))
        dec.append(str(float(data[2]) * np.pi / 180))
        comm.append(generateComment(data))

    nameCol = fits.Column(name="NAME", format="80A", array=name)
    typeCol = fits.Column(name="TYPE", format="1A", array=t)
    raCol = fits.Column(name="RA", format="1D", array=ra)
    decCol = fits.Column(name="DEC", format="1D", array=dec)
    commCol = fits.Column(name="COMMENT", format="80A", array=comm)

    outCols = fits.ColDefs([nameCol, typeCol, raCol, decCol, commCol])
    return fits.BinTableHDU().from_columns(outCols, name="fibres")


# -------------------------------- ** -------------------------------- #


def generateComment(DBData):
    """
    Generates the comment string in the Fibres Extension.
    """
    t = DBData[3]
    z = DBData[4]
    tf = "P" if DBData[5] == "" else DBData[5]
    qf = DBData[6]
    n = DBData[7]
    return str(t) + " " + str(z) + " " + tf + qf + " - " + n


# -------------------------------- ** -------------------------------- #
# Actual FITS parsing: retrieves information based on the instrument.  #
# Same structure, sometimes cards have different names thus different  #
# functions.                                                           #
# Works provided an instrument is given, otherwhise it can't set the   #
# correct cards.                                                       #
# -------------------------------- ** -------------------------------- #


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
    Parses information from a generic HDU. Very likely to fail, will hopefully be improved!
    """
    wave = []
    flux = []
    for w, f in hdul[1].data:
        wave.append(w)
        flux.append(f)

    wave = np.array(wave)
    flux = np.array(flux)

    wave = np.reshape(wave, (1, wave.shape[0]))
    flux = np.reshape(flux, (1, flux.shape[0]))
    error = flux * 0.1

    return (wave, flux, error)


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

    wave = data.x * 10
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
    error = flux * .1
    return (wave.reshape(1, -1), flux.reshape(1, -1), error.reshape(1, -1))

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

if __name__ == "__main__":
    MarzConverter(sysargs=sys.argv)
