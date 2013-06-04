"""
Runs sextractor-based photometry on HIPPIES data.
Run from the directory where the par*_W directories live. Also needs a
sextractor/ directory in the same place with photometry.conf file inside.

Example usage:

$> python 04_hippies_phot.py clean

Command-line options:

clean       Delete old catalogs, forcing sextractor re-run
help        Display this message
"""
from __future__ import division
import subprocess
import os
import glob
import sys
import pyfits
import SextractorTools as SexTools
from HippiesConfig import *
from StringIO import StringIO
from numpy import savetxt
from numpy.lib.recfunctions import append_fields


def parse_cl_args(argv):
    """
    Processes the command line arguments supplied to hippies_process_new,
    replacing defaults, checking for file existence, etc.
    @param argv: Command-line arguments sent to hippies_process_new
    @return: A dictionary of the command-line arguments
    """
    argv = [arg.strip('-') for arg in argv]
    _known_args = ['clean', 'help']
    args = dict([(arg, arg in argv) for arg in _known_args])
    return args


def sextractField(field):
    os.chdir(field)
    sciFiles = glob.glob(field+'*drz_sci.fits')
    # Find which is reference filter. If none found, skip this field
    refFile = ''
    for refFilt in PHOT_DETECT_FILTERS:
        for sciFile in sciFiles:
            if refFilt in sciFile:
                refFile = sciFile
    if refFile == '':
        message('Unable to find reference filter for {}'.format(field),
                msgtype='WARN')
        os.chdir('..')
        return

    # Run sextractor on each field
    for sciFile in sciFiles:
        sextractFile(sciFile, detectFile=refFile)

    # Combine individual sextractor catalogs
    message('Combining sextractor catalogs for field: {}'.format(field))
    combineSextractorCatalogs(field)

    os.chdir('..')


def sextractFile(sciFile, detectFile=None):
    rmsFiles = sciFile.replace('_sci', '_rms')
    sciFiles = sciFile
    weight_thresh = str(RMS_BAD_PX_VALUE)

    ## comma-separate detection and photometry images if dual-image mode
    if detectFile is not None:
        rmsFiles = detectFile.replace('_sci', '_rms') + ',' + rmsFiles
        sciFiles = detectFile + ',' + sciFile
        weight_thresh = weight_thresh + ',' + str(RMS_BAD_PX_VALUE)

    catFile = sciFile.replace('.fits', '.cat')
    if os.path.exists(catFile):
        message('Found catalog {}. Use "clean" to force re-run'.format(catFile))
        return catFile

    filt = pyfits.getval(sciFile, 'FILTER')
    try:
        gain = pyfits.getval(sciFile, 'CCDGAIN') * pyfits.getval(sciFile,
                                                                 'EXPTIME')
    except KeyError:
        gain = pyfits.getval(sciFile, 'ADCGAIN') * pyfits.getval(sciFile,
                                                                 'EXPTIME')
    magZP = FILTER_ZEROPOINTS[filt]

    cmd = [SEX_COMMAND, sciFiles, '-c', '../'+SEX_PHOTOMETRY_CONF,
           '-weight_image', rmsFiles, '-weight_thresh', weight_thresh,
           '-catalog_name', catFile, '-mag_zeropoint', str(magZP),
           '-gain', str(gain)]
    subprocess.call(cmd)

    ## Ensure the catalog file was created
    if not os.path.exists(catFile):
        raise IOError('Unable to find catalog. Suspected Sextractor failure.')

    return catFile


def combineSextractorCatalogs(field):
    catList = glob.glob(field+'*drz_sci.cat')
    filters = [cat.split('_')[2] for cat in catList]
    firstCat = SexTools.read_catalog(catList[0])
    filtCols = []
    for filt, catFile in zip(filters[1:], catList[1:]):
        newCat = SexTools.read_catalog(catFile)
        if (firstCat['NUMBER'] != newCat['NUMBER']).any():
            message('Catalogs for field {} come from '.format(field) +
                    'different sextractor runs. Delete catalogs and run again.',
                    msgtype='ERROR')
            return
        for col in newCat.dtype.names:
            if ((firstCat[col] != newCat[col]).any() and
                        col not in IGNORE_DUPLICATE_COLUMNS):
                filtCols += [col]
                firstCat = append_fields(firstCat, names=col+'_'+filt,
                                         data=newCat[col], asrecarray=True)
    filtCols = set(filtCols)
    message('Filter-specific columns found: {}'.format(', '.join(filtCols)))
    # Rename columns that were filter-specific for the original array
    firstFilt = filters[0]
    firstCat.dtype.names = [name+'_'+firstFilt if name in filtCols else name
                            for name in firstCat.dtype.names]
    masterFilename = field+'_master.cat'
    outData = StringIO('')
    savetxt(outData, firstCat, fmt='%s')
    outData.seek(0)
    with open(masterFilename, 'w') as f:
        f.write('#' + '\t'.join(firstCat.dtype.names) + '\n')
        f.write(outData.read())


if __name__ == '__main__':
    ## Process command-line arguments
    args = parse_cl_args(sys.argv)

    if args['help']:
        print(__doc__)
        exit(0)

    if args['clean']:
        cleanFiles = glob.glob(FIELD_PREFIX+'*W/*.cat')
        message('Cleaning old catalog files.')
        for cleanFile in cleanFiles:
            try:
                os.remove(cleanFile)
            except OSError:
                pass

    fields = glob.glob(FIELD_PREFIX+'*W')
    ## Run sextractor on everything
    for field in fields:
        message('Processing field: {}'.format(field))
        sextractField(field)

