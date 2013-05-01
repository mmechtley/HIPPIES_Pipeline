"""
Processes new observations for the HIPPIES HST program.
Please run the script from the hippies data root directory, which should
contain two directories: sextractor/ and unprocessed/

Example usage:
$> python 01_hippies_process_new.py overwrite

Available command-line options:
overwrite      Overwrite existing files? Default is no.
help           Display this message
"""

import sys                # for arguments
import pyfits            # for reading fits
import glob                # for finding all flt files
import os                # for creating dirs, getting working dir, etc
import shutil
from HippiesConfig import *


def parse_cl_args(argv):
    """
    Processes the command line arguments supplied to hippies_process_new,
    replacing defaults, checking for file existence, etc.
    @param argv: Command-line arguments sent to hippies_process_new
    @return: A dictionary of the command-line arguments
    """
    argv = [arg.strip('-') for arg in argv]
    _known_args = ['overwrite', 'help']
    args = dict([(arg, arg in argv) for arg in _known_args])
    return args


def is_processed(filtDir, baseName):
    """
    Determines if the field has been processed by checking for existing files
    @param filtDir: Path to the observation. e.g. ../parXXX_W/F105W
    @param baseName: Base name of the observation. e.g. ib8d0aatq
    @return: True if the observation has been processed, else False
    """
    return os.path.isfile('{}/do_not_touch/{}_flt.fits'.format(filtDir,
                                                               baseName))


def construct_field_name(instrumentName, aperRA, aperDec):
    """
    Returns the default field name string for a given file handle. Constructed
    as first two elements of sexigesimal format
    """
    ## Grab first two components of sexigesimal format
    raHr = int(aperRA / 15)
    decDeg = int(aperDec)
    decMin = int(abs(aperDec - decDeg) * 60)
    fRA = '{:02d}{:02d}'.format(raHr, int(abs(aperRA * 4 - raHr * 60)))
    fDec = '{}{:02d}{:02d}'.format(('+' if aperDec >= 0.0 else '-'),
                                   abs(decDeg), decMin)

    return 'par' + fRA + fDec + '_' + instrumentName[0]


def find_field_dir(aperRA, aperDec, allFieldsDict):
    """
    Given a field center RA and Dec, and a dictionary of fieldName:
    [(ra, dec), ...], find which one this RA/Dec belongs to
    @param aperRA: right ascension coordinate of aperture center (degrees)
    @param aperDec: declination coordinate of aperture center (degrees)
    @param allFieldsDict: dictionary of {fieldName : [(ra, dec), ...]}
    @return: directory name of field aperRA/aperDec belongs to
    """
    foundDir = ''
    for dirName, raDecPairs in allFieldsDict.iteritems():
        if fields_are_same([(aperRA, aperDec)], raDecPairs):
            foundDir = dirName
    # Something is horribly wrong if we can't find which field this file
    # belongs to
    assert foundDir != ''
    return foundDir


def combine_adjacent_fields(fieldDict):
    """
    Combines pointings that are spatially close into single fields
    @param fieldDict: dictionary with form {fieldName: [(ra1,dec1), ...]}
    @return: dictionary with same form, but adjacent fields combined
    """
    outDict = fieldDict.copy()
    for f1Name, f1Coords in fieldDict.iteritems():
        for f2Name, f2Coords in fieldDict.iteritems():
            if (f1Name != f2Name and fields_are_same(f1Coords, f2Coords) and
                    f1Name in outDict and f2Name in outDict):
                outDict[f1Name] = outDict[f1Name] + f2Coords
                del outDict[f2Name]
    return outDict


def fields_are_same(coordList1, coordList2):
    """
    Checks if two fields (with our naming scheme) are closely adjacent
    @param coordList1: list of 2-tuples containing field pointing centers
    @param coordList2: Same as coordList1, for second field
    """
    permutations = [(c1, c2) for c1 in coordList1 for c2 in coordList2]
    for c1, c2 in permutations:    # c1, c2 are 2-tuple
        sqrdist = (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2
        if sqrdist < FIELD_CONNECT_DISTANCE ** 2:
            return True
    return False


if __name__ == '__main__':
    # Put command-line arguments into a nice little dictionary
    args = parse_cl_args(sys.argv)

    if (os.path.exists('sextractor') and os.path.exists('pipeline') and
            not os.path.exists('unprocessed')):
        message("unprocessed/ directory does not exist. I'll attempt to" +
                "create it for you. Dump all your groovy *_flt.fits files " +
                "inside, then re-run this script.", msgtype='ERROR')
        try:
            os.mkdir('unprocessed')
        except OSError:
            pass
        exit(1)

    # First, make sure we're running this from the right place
    if 'unprocessed' not in os.listdir('.') or args['help']:
        print(__doc__)
        exit(1)

    os.chdir('unprocessed')

    # Find all _flt fits files
    filesToProcess = glob.glob('*_flt.fits')

    # Get a list of all files, to find adjacency (match par*)
    proc_pattern = '../{}*/F*/do_not_touch/*_flt.fits'.format(FIELD_PREFIX)
    allFiles = filesToProcess + glob.glob(proc_pattern)

    # A dictionary of [short field name] : [(ra,dec)]
    message('Collecting field location information')
    allFields = dict()
    for flt_file in allFiles:
        instrumentName = pyfits.getval(flt_file, 'INSTRUME', ext=0)
        aperRA = pyfits.getval(flt_file, 'RA_APER', ext=1)
        aperDec = pyfits.getval(flt_file, 'DEC_APER', ext=1)

        fieldName = construct_field_name(instrumentName, aperRA, aperDec)
        allFields[fieldName] = [(aperRA, aperDec)]

    message('Calculating field adjacency')
    oldFields = {}
    while set(allFields.keys()) != set(oldFields.keys()):
        oldFields = allFields.copy()
        allFields = combine_adjacent_fields(allFields)

    for flt_file in filesToProcess:
        primHeader = pyfits.getheader(flt_file, ext=0)

        baseName = primHeader['ROOTNAME']
        instName = primHeader['INSTRUME']
        detector = primHeader['DETECTOR']
        aperRA = pyfits.getval(flt_file, 'RA_APER', ext=1)
        aperDec = pyfits.getval(flt_file, 'DEC_APER', ext=1)

        # Skip files that aren't WFC3
        if instName != 'WFC3':
            message('File {} is not a WFC3 observation, skipping it.'.format(
                flt_file), msgtype='WARN')
            continue

        # Skip files that don't contain science data
        if primHeader['FILETYPE'] != 'SCI':
            message('File {} '.format(flt_file) +
                    'is not a science file. Perhaps it is named incorrectly?',
                    msgtype='WARN')
            continue

        # Find the filter used for the file
        filt = ''
        for key in ['FILTER', 'FILTER1', 'FILTER2']:
            if primHeader.get(key, default='').startswith('F'):
                filt = primHeader[key]
        # If file had no FILTER keyword whose value starts with F, something is
        # badly amiss
        if filt == '':
            message('No filter name found for file {}.'.format(flt_file) +
                    'Ensure header information is correct', msgtype='WARN')
            continue

        # Search the field location dictionary to find out which field dir to
        # place this file in
        fieldDir = find_field_dir(aperRA, aperDec, allFields)
        obsDir = os.path.join('../../', fieldDir)
        filtDir = os.path.join(obsDir, filt)

        # Check to make sure this observation has not been processed
        if is_processed(filtDir, baseName) and not args['overwrite']:
            message('File {} '.format(flt_file) +
                    'has already been processed, skipping it. To force ' +
                    'reprocessing, run again with overwrite', msgtype='WARN')
            continue
        else:
            message('Processing observation {} in {}'.format(flt_file, obsDir))

        # Create directory tree for this observation set. Ignore errors that
        # complain of already existing dirs
        try:
            os.makedirs(os.path.join(filtDir, 'do_not_touch'))
        except OSError:
            message('Directory {} already exists.'.format(filtDir))

        # Copy the unprocessed file to the /do_not_touch/ dir
        shutil.move(flt_file, os.path.join(filtDir, 'do_not_touch', flt_file))

    message('Processing complete. To run multidrizzle on all observations, ' +
            'run 02_hippies_drizzle_all.py')
