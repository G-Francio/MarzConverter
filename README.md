# MarzConverter
Simple script to convert spectra to a [Marz](https://skymapper.anu.edu.au/static/sm_asvo/marz/index.html#/overview) compliant format. It is designed to pull data from the QUBRICS internal DB, but can be used without access to it.

Code is probably terrible, not a software dev!

The script was written and tested on Linux Ubuntu with Python 3 and will likely not work with Windows or Python 2.x; it might work on Mac.

Requires:
* [Astropy](https://www.astropy.org/) (with its requirements)
* (optional) [Mariadb Python Connector](https://mariadb.com/resources/blog/how-to-connect-python-programs-to-mariadb/) (with its requirements)

## Usage
The script requires a mandatory input argument (the input file) and some optional arguments. It produces a single output file which includes all converted spectra (making use of the `FIBRE` FITS extension allowed by Marz).
The first argument is mandatory and is the input spectrum or list of spectra. Accepted formats are:
* `.fits` if a single spectra is given;
* `.txt`, `.dat`, `.mzc` text file if multiple spectra need conversion (see below or the example for format).
By default, converted spectra are saved in the same folder the original spectrum is in.

Optional arguments should be passed by keyword. Valid arguments are:
* `outfile`: allow to specify the name and path of the output file
* `wr`: allows to specify a wavelenght range (only useful for single .fits files)
* `pooling`: enables (`pooling=True`) or disables (`pooling=False`, default behaviour)
* `npool`: if pooling is enabled, spcifies the number of spawned sub-processes. Default is half the number of available CPUs (`os.cpu_count//2`)

If the input is a text file, it should contain a list of spectra (one per line), together with the (optional) appropriate wavelength range. Example of a possibile text file is reported below, or in `TextExample.mzc`:

```
/path/to/first/spec.fits
/path/to/second/spec.fits
/path/to/third/spec.fits wr = [4000, 9000]
```

The wavelength range can be set using the following: `"wr = [w_1, w_2]"`, or, without quotes as `wr=[w_1,w_2]`. The same syntax is used if a list of spectra is given by appending the same string to the spectrum path inside the text file list. Each spectrum can be given an appropriate wavelength range.

The wavelength range has to be given in appropriate units (same as the input spectrum).

The variable `QUBRICSACCESS` can be set to `False` to stop the script from repeatedly asking the user credentials when running.

### Example
`python MarzConverter.py /path/to/spec.fits wr=[w_1,w_2] outfile=/path/to/outSpec`  
`python MarzConverter.py /path/to/specList.txt outfile=/path/to/outSpec`
