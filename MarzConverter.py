#!/usr/bin/python3

# Logging, see here: https://stackoverflow.com/questions/3383865/how-to-log-error-to-file-and-not-fail-on-exception
# and here: https://stackoverflow.com/questions/8050775/using-pythons-logging-module-to-log-all-exceptions-and-errors
# and here: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python/16993115#16993115

import logging, sys

logging.basicConfig(filename = '/tmp/MarzConverter.log', level = logging.DEBUG, 
                    format = '%(asctime)s %(levelname)s %(name)s %(message)s')

logger  = logging.getLogger(__name__)
handler = logging.StreamHandler(stream = sys.stdout)
logger.addHandler(handler)

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.critical("Uncaught exception, consider reporting the log (/tmp/MarzConverter.log).", exc_info = (exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

from astropy.io import fits
import numpy as np
import os.path as p, re

# TODO: Add handling of external errors
# TODO: Check if fallback for no INSTRUME cards is enough
# TODO: Centre on screen if cut wave -> Theoretically correct, in practice does not work

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
    return [name, '0', '0', '-', '-', '-', '-', '-']

# -------------------------------- ** -------------------------------- #

try:
    import mariadb as mdb

    def getUserCredential():
        """
        Read user credentials for QubricsDB.
        Requires the cred.auth file in the same folder as MarzConverter,
        or manually setting a different path.
        """
        with open('./cred.auth') as f:
            cred = f.readlines()
        return cred[0].strip('\n'), cred[1].strip('\n')

    # -------------------------------- ** -------------------------------- #

    def DBConnect(user, pw):
        """
        Connects to QubricsDB, returns a cursor and the connection object.
        """
        conn = mdb.connect(user = user, password = pw, host = "127.0.0.1", port = 3306, database = "Qubrics")
        cur  = conn.cursor()
        return cur, conn

    # -------------------------------- ** -------------------------------- #

    def getObservationData(nameList):
        """
        Gathers data for observations that needs conversion.
        `nameList` is a list of fileNames taken from Drive/QSOCandidates.
        """
        _cred = getUserCredential()

        #index = []
        observationData = []
        fileNameList = []

        try:
            cur, conn = DBConnect(_cred[0], _cred[1])
            cur.execute('SELECT file FROM Qubrics.Observations')
            for c in cur:
                fileNameList.append(c)
            fileNameList = np.array(fileNameList)[:, 0]

            for name in nameList:
                foundIndices = np.where(np.char.find(fileNameList, '/' + name) > 0)
                if len(foundIndices[0]) > 1:
                    print("Multiple matches for {}, using fallback data.".format(name))
                    observationData.append(getFallbackDataSingle(name))
                elif len(foundIndices[0]) == 0:
                    print("Can't find spec data on DB for {}, fallback on mock data.".format(name))
                    observationData.append(getFallbackDataSingle(name))
                else:
                    fileName = fileNameList[foundIndices[0][0]]
                    cur.execute("SELECT target_qid, RAd, DECd, objtype, z_spec, targetflag, qflag, notes FROM Qubrics.Observations WHERE file = '{}'".format(fileName))
                    for c in cur:
                        observationData.append(c)

            conn.close()
            return np.array(observationData)
        except mdb.OperationalError:
            print("Connection to the DB failed, mock data will be used")
            print("Connect to the DB before running the script for possible additional data.")
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
    args = kwargs['sysargs']
    if (any('txt' in s for s in args) or
        any('dat' in s for s in args) or
        any('mzc' in s for s in args)):
        multiFits2File(args)
    else:
        fits2File(args)

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

def fits2File(args):
    """
    Reads a FITS file and calls the appropriate function.
    Writes the resulting FITS to file, ready for MARZ.
    If error is not found the original FITS, error is assumed .1
    of the original flux.
    """
    fitsIn    = args[1]
    waveRange = parseWR(args[2]) if any('wr' in s for s in args) else None

    path, _ = p.splitext(fitsIn)
    name = path.split('/')[-1] + '_Marz.fits'

    specDBData = getObservationData([path.split('/')[-1]])
    fibreHDU   = generateFibresData(specDBData)

    # Sets the correct outfile:
    if len(args) == 3 and not any('wr' in s for s in args):
        args[-1] = args[-1].split('.fits')[0]
        name = args[-1] + '.fits'
    elif len(args) == 4:
        args[-1] = args[-1].split('fits')[0]
        name = args[-1] + '.fits'

    hduList = fits.open(fitsIn)

    wave, flux, error = parseData(hduList, fitsIn)
    hduList.close()

    if waveRange is not None:
        wave, flux, error = cutWavelength(wave, flux, error, waveRange)

    writeFits(flux, error, wave, fibre = fibreHDU, name = name)

# -------------------------------- ** -------------------------------- #

def multiFits2File(args):
    """
    Reads a file list and calls the appropriate function for each fits.
    """
    path, _ = p.splitext(args[1])
    name = path.split('/')[-1] + '_Marz.fits' if len(args) < 3 else args[2] + '.fits'
    
    waveList, fluxList, errorList, nameList = [], [], [], []
    specFiles = readSpecList(args[1])
    
    for spec in specFiles:
        specFileName = p.splitext(spec[0])[0].split('/')[-1]
        wave, flux, error = fits2array(spec[0], waveRange = spec[1])
        nameList.append(specFileName)
        waveList.append(wave)
        fluxList.append(flux)
        errorList.append(error)

    specDBData = getObservationData(nameList)
    fibreHDU   = generateFibresData(specDBData)

    maxShape = max([s.shape[1] for s in waveList])
    waveList, fluxList, errorList = padArray(waveList, fluxList, errorList, maxLength = maxShape)
    writeFits(fluxList, errorList, waveList, fibre = fibreHDU, name = name)

# -------------------------------- ** -------------------------------- #

def fits2array(fitsIn, waveRange = None):
    """
    Reads a FITS file and calls the appropriate function.
    Returns `wave`, `flux` and `error`.
    If error is not found the original FITS, error is assumed .1
    of the original flux.
    """

    hduList = fits.open(fitsIn)

    wave, flux, error = parseData(hduList, fitsIn)
    hduList.close()

    if waveRange is not None:
        return cutWavelength(wave, flux, error, waveRange)

    return wave, flux, error

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

def parseData(hdul, fileName):
    """
    Parses a FITS file given the instrument name.
    Very likely to fail if the fits is non standard, please report the log files
    if errors are encountered, possibly with the problematic header/spectrum.
    """
    header = hdul[0].header
    try:
        inst = header['INSTRUME']
    except KeyError:
        inst = None
    
    if inst == 'WFCCD/WF4K-1':
        return parseWFCCD(hdul)
    elif re.search('LDSS3-.*', inst) is not None:
        return parseLDSS3(hdul)
    elif inst == 'IMACS Short-Camera':
        return parseIMACS(hdul)
    elif inst == 'FIRE':
        return parseFIRE(hdul)
    elif inst == 'EFOSC':
        return parseEFOSC(hdul)
    elif inst == 'LRS': # TNG
        return parseLRS(hdul)
    elif inst is None:
        return parseGeneric(hdul)
    else:
        print("Can't parse fits file, please report log file (/tmp/MarzConverter.log)")

# -------------------------------- ** -------------------------------- #

def cutWavelength(wave, flux, error, waveRange):
    """
    Reduces a spectrum to a given wavelength range.
    """
    wr          = waveRange
    reduced_wav = np.argwhere((wave >= wr[0]) & (wave <= wr[1]))
    wave        = wave[:, reduced_wav[:, 1]]
    flux        = flux[:, reduced_wav[:, 1]]
    error       = error[:, reduced_wav[:, 1]]

    return wave, flux, error

# -------------------------------- ** -------------------------------- #

def padArray(*args, maxLength = 0):
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

def padSingleArray(array, ml = 0):
    """
    Performs the padding operation.
    """
    paddedArray = np.zeros((1, ml))
    paddedArray[:array.shape[0], :array.shape[1]] = array
    return paddedArray

# -------------------------------- ** -------------------------------- #

def readSpecList(fileIn):
    """
    Reads a list of spectra to process.
    """
    with open(fileIn) as f:
        readData = [line.strip('\n') for line in f.readlines()]

    spec2Convert = []
    for data in readData:
        splitArgs = data.split('wr')
        spec = splitArgs[0].strip(' "')
        wr   = parseWR(splitArgs[1]) if len(splitArgs) > 1 else None
        spec2Convert.append([spec, wr])

    return np.array(spec2Convert, dtype = object)

# -------------------------------- ** -------------------------------- #

def writeFits(flux, error, wave, fibre = None, name = 'MarzConverterOutput.fits'):
    """
    Writes the process fits file, ready for Marz. Asks for overwrite permission!
    """
    primaryHDU  = fits.PrimaryHDU(flux)
    varianceHDU = fits.ImageHDU(error, name = 'variance')
    waveHDU     = fits.ImageHDU(wave,  name = 'wavelength')

    if fibre is None:
        hduListOut = fits.HDUList([primaryHDU, varianceHDU, waveHDU])
    else:
        hduListOut = fits.HDUList([primaryHDU, varianceHDU, waveHDU, fibre])

    try:
        hduListOut.writeto(name)
        hduListOut.close()
    except OSError:
        overwrite = input("File already exists, overwrite (Y/n)? ")
        if overwrite.lower() == "y":
            hduListOut.writeto(name, overwrite = True)
            hduListOut.close()
        else:
            hduListOut.close()

# -------------------------------- ** -------------------------------- #

def parseWR(str):
    """
    Parses the wavelength range.
    """
    wr = str.strip('wr=[ ]"')
    return [float(i.strip(' ')) for i in wr.split(',')]

# -------------------------------- ** -------------------------------- #

def generateFibresData(DBData):
    """
    Given data from QDB, produces fibre data.
    If no data are found on the DB (or the DB can't be accessed) neutral, mock
    data are generated on the fly (e.g. ra/dec = 0/0, `z_spec = -`)
    """
    name, t, ra, dec, comm = [], [], [], [], []
    for data in DBData:
        name.append(data[0])
        t.append('P')
        ra.append(str(float(data[1])*np.pi/180))
        dec.append(str(float(data[2])*np.pi/180))
        comm.append(generateComment(data))

    nameCol = fits.Column(name = 'NAME',     format = '80A', array = name)
    typeCol = fits.Column(name = 'TYPE',     format = '1A',  array = t)
    raCol   = fits.Column(name = 'RA',       format = '1D',  array = ra)
    decCol  = fits.Column(name = 'DEC',      format = '1D',  array = dec)
    commCol = fits.Column(name = 'COMMENT',  format = '80A', array = comm)

    outCols = fits.ColDefs([nameCol, typeCol, raCol, decCol, commCol])
    return fits.BinTableHDU().from_columns(outCols, name = 'fibres')

# -------------------------------- ** -------------------------------- #

def generateComment(DBData):
    """
    Generates the comment string in the Fibres Extension.
    """
    t  = DBData[3]
    z  = DBData[4]
    tf = 'P' if DBData[5] == '' else DBData[5]
    qf = DBData[6]
    n  = DBData[7]
    return t + ' ' + z + ' ' + tf + qf + ' - ' + n

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
    start = hdul[0].header['CRVAL1']
    step  = hdul[0].header['CD1_1']
    total = hdul[0].header['NAXIS1']
    corr  = (hdul[0].header['CRPIX1'] - 1) * step

    wave  = np.arange(start - corr, start + total*step - corr, step)
    wave  = np.reshape(wave, (1, wave.shape[0]))
    flux  = np.reshape(hdul[0].data[0], (1, hdul[0].data[0].shape[1]))
    error = np.reshape(hdul[0].data[1], (1, hdul[0].data[1].shape[1]))

    return (wave, flux, error)


def parseLDSS3(hdul):
    """
    Parses information from a given HDU, for data produced at WFCCD
    """
    start = hdul[0].header['CRVAL1']
    step  = hdul[0].header['CDELT1']
    total = hdul[0].header['NAXIS1']

    corr  = (hdul[0].header['CRPIX1'] - 1) * step

    wave  = np.arange(start - corr, start + total*step - corr, step)
    wave  = np.reshape(wave, (1, wave.shape[0]))
    flux  = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * .1

    return (wave, flux, error)


def parseIMACS(hdul):
    """
    Parses information from a given HDU, for data produced at IMACS
    """
    start = hdul[0].header['CRVAL1']
    step  = hdul[0].header['CDELT1']
    total = hdul[0].header['NAXIS1']
    corr  = (hdul[0].header['CRPIX1'] - 1) * step

    wave  = np.arange(start - corr, start + total*step - corr, step)
    wave  = np.reshape(wave, (1, wave.shape[0]))
    flux  = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * .1

    return (wave, flux, error)


def parseFIRE(hdul):
    """
    Parses information from a given HDU, for data produced at FIRE
    """
    data = hdul[5].data

    wave  = data.field('WAVE')
    flux  = data.field('FLUX')
    error = data.field('SIG')

    return (wave, flux, error)


def parseEFOSC(hdul):
    """
    Parses information from a given HDU, for data produced at EFOSC
    """
    start =  hdul[0].header['CRVAL1']
    step  =  hdul[0].header['CDELT1']
    total =  hdul[0].header['NAXIS1']
    corr  = (hdul[0].header['CRPIX1'] - 1) * step

    wave  = np.arange(start - corr, start + total*step - corr, step)
    wave  = np.reshape(wave, (1, wave.shape[0]))
    flux  = np.reshape(hdul[0].data, (1, hdul[0].data.shape[0]))
    error = flux * .1

    return (wave, flux, error)


def parseLRS(hdul):
    """
    Parses information from a given HDU, for data produced at TNG LRS
    """
    start =  hdul[0].header['CRVAL1']
    step  =  hdul[0].header['CDELT1']
    total =  hdul[0].header['NAXIS1']
    corr  = (hdul[0].header['CRPIX1'] - 1) * step

    wave  = np.arange(start - corr, start + total*step - corr, step)
    r_wav = np.argwhere((wave >= 3700) & (wave <= 8000)) #reduced_wave,
    # TNG spectra are very noisy at the extremes of the wavelength range

    wave  = np.reshape(wave[r_wav][:, 0], (1, wave[r_wav][:, 0].shape[0]))
    flux  = np.reshape(hdul[0].data[r_wav][:, 0], (1, hdul[0].data[r_wav][:, 0].shape[0]))
    error = flux * .1

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
    
    wave  = np.reshape(wave, (1, wave.shape[0]))
    flux  = np.reshape(flux, (1, flux.shape[0]))
    error = flux * .1

    return (wave, flux, error)

# -------------------------------- ** -------------------------------- #
# #################################################################### #
# -------------------------------- ** -------------------------------- #

if __name__ == '__main__':
    MarzConverter(sysargs = sys.argv)

