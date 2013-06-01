"""
Drizzles all observations for HST HIPPIES project. Run from the directory where
all the parXXXX+XXXX_W directories live.
Requires a /sextractor/ directory as a sibling to the parXXXX+XXXX_W
directories, which must include the sextractor parameter files for finding
compact objects (getcompact.conf)

Available command-line options:
clean      Clean existing files first? Default is no.
help       Display this message
"""

import os
import shutil
import re
import sys
import glob
import subprocess
import multidrizzle
import pyfits
from numpy import seterr
from HippiesConfig import *
import CatalogTools as CatTools
import SextractorTools as SexTools

# Additional Multidrizzle parameters for single-image drizzle step. Do not
# modify these
_driz_params_single = {'static': False, 'skysub': False, 'median': False,
                       'blot': False, 'driz_cr': False, 'driz_combine': False,
                       'clean': False}


def drizzle_field(fieldPath, sys_args=None):
    """
    Drizzle a field, beginning with WFC3 IR F105W or F098M and then drizzling
    other filters to the same grid.
    @param fieldPath: Base path to the field to be drizzled.
    @param sys_args: dictionary of command-line supplied arguments
    """
    ## Get list of filter dirs that are children of this dir
    filts = filter_list(fieldPath)
    ## Set which filter will be the reference
    refFilt = ref_filter(filts)
    if refFilt == '':
        message('Field {} '.format(fieldPath) +
                'has no suitable reference filter (F105W or F098M).',
                msgtype='WARN')
        return

    # Drizzle the other filters to the same grid
    for filt in filts:
        filtPath = os.path.join(fieldPath, filt)

        # Clean directory if needed
        if is_processed(filtPath) and sys_args['clean']:
            clean_dir(filtPath)

        # Drizzle if the dir is clean, otherwise just warn that we won't process
        if not is_processed(filtPath):
            message('Beginning drizzle process for {}'.format(filtPath))

            refFile = ''
            if filt != refFilt:
                refFile = '_'.join([fieldPath, refFilt, 'drz_sci.fits'])
                refFile = os.path.join(fieldPath, refFilt, refFile)
                # Symlink the reference file to this filter's dir
                srcsci = os.path.relpath(refFile, filtPath)
                destsci = os.path.join(filtPath,
                                       'ref_' + os.path.basename(refFile))
                srcweight = srcsci.replace('sci', 'weight')
                destweight = destsci.replace('sci', 'weight')
                os.symlink(srcsci, destsci)
                os.symlink(srcweight, destweight)
                refFile = os.path.basename(destsci)

            drizzle_filt(filtPath, ref_file=refFile, sys_args=sys_args)
        else:
            message('Filter {} has already '.format(filtPath) +
                    'been drizzled. To re-drizzle, run with "clean."')

    ## Trim off last line or column to deal with Multidrizzle bug
    trim_files(fieldPath)
    return


def drizzle_filt(filt_path, ref_file='', sys_args=None):
    """
    Runs multidrizzle on the selected filter using HIPPIES default settings, or
    settings from supplied file
    @param filt_path: Path to the filter directory
    @param ref_file: Optional file to use as the reference grid
    @return: Filename of the final drizzled science image
    """
    old_dir = os.getcwd()
    out_file = base_filename(filt_path)

    ## Change to working dir
    os.chdir(filt_path)

    out_file = os.path.join(os.getcwd(), out_file)

    flt_pattern = '*_flt.fits'
    flt_files = glob.glob('./do_not_touch/'+flt_pattern)
    instInfo = instrument_info(flt_files[0])

    ref_file = os.path.basename(ref_file)

    # Skip anything that's not WFC3
    if instInfo['instrument'] == 'WFC3':
        # Copy files to working directory
        for src in flt_files:
            shutil.copy2(src, './')
        flt_files = glob.glob(flt_pattern)

        # Run first drizzle pass, doing only driz_separate to get single_sci
        # files
        driz_params_single = dict(DRIZ_PARAMS_WFC3.items() +
                                  _driz_params_single.items())
        ## For slave filters (reference file is a drizzled file)
        if os.path.isfile(ref_file):
            driz_params_single['refimage'] = ref_file

        md = multidrizzle.Multidrizzle(input=flt_pattern,
                                       output=out_file,
                                       **driz_params_single)
        md.build()
        md.run()

        # If no ref image is supplied (reference filter) use first single_sci
        # image
        singleScis = glob.glob('*_single_sci.fits')
        if ref_file == '':
            ref_file = 'refimage_sci.fits'
            shutil.copy2(singleScis[0], ref_file)
            shutil.copy2(singleScis[0].replace('sci', 'wht'),
                         ref_file.replace('sci', 'weight'))

        ## Calculate delta shifts for all the single sci files
        ## TODO: Build tree of image + refimage, by distance. Calculate shift
        # from one to next, set absolute shift as sum of relative shifts of
        # chain
        shifts = [calc_shift(ref_file, singleSci) for singleSci in singleScis]

        ## Write out shift file and delete all other files
        write_shift_file(ref_file=ref_file, file_list=flt_files,
                         shift_list=shifts)
        keep_files = glob.glob('*.cat') + ['shifts.txt', ref_file]
        clean_dir('.', exclude=keep_files)
        ## Re-copy FLTs
        for src in glob.glob('./do_not_touch/' + flt_pattern):
            shutil.copy2(src, './')

        ## Now redrizzle with shifts
        md = multidrizzle.Multidrizzle(input=flt_pattern,
                                       output=out_file,
                                       refimage=ref_file,
                                       shiftfile='shifts.txt',
                                       **DRIZ_PARAMS_WFC3)
        md.build()
        md.run()

    sci_file = out_file + '_drz_sci.fits'

    ## Return to original working dir
    os.chdir(old_dir)

    return sci_file


def calc_shift(refFile, otherFile):
    """
    Runs sextractor to find compact objects, then finds the shift between the
    compact object catalogs
    """
    refCatFile = refFile.replace('.fits', '.cat')
    otherCatFile = otherFile.replace('.fits', '.cat')
    sex_conf = os.path.join('..', '..', SEX_COMPACT_CONF)

    ## Find weight files
    refWeight = refFile.replace('sci', 'weight')
    otherWeight = otherFile.replace('sci', 'weight')
    if not os.path.exists(refWeight):
        refWeight = refFile.replace('sci', 'wht')
    if not os.path.exists(otherWeight):
        otherWeight = otherFile.replace('sci', 'wht')

    message('Calculating shift for file {} using reference {}'.format(
        otherFile, refFile))
    for f in (refFile, otherFile, sex_conf):
        if not os.path.exists(f):
            raise OSError('Could not find file {}. '.format(f) +
                          'While working in {}'.format(os.getcwd()) +
                          'Unable to calculate shift.')

    ## Run sextractor for reference file
    if not os.path.exists(refCatFile):
        cmd = [SEX_COMMAND, refFile, '-c', sex_conf, '-weight_image',
               refWeight, '-catalog_name', refCatFile]
        subprocess.call(cmd)

    ## Run sextractor for file to shift
    if not os.path.exists(otherCatFile):
        cmd = [SEX_COMMAND, otherFile, '-c', sex_conf, '-weight_image',
               otherWeight, '-catalog_name', otherCatFile]
        subprocess.call(cmd)

    refCat = SexTools.read_catalog(refCatFile,
                                   keepCols=['X_IMAGE', 'Y_IMAGE'])
    refCat = refCat.view(float).reshape(-1, 2)
    otherCat = SexTools.read_catalog(otherCatFile,
                                     keepCols=['X_IMAGE', 'Y_IMAGE'])
    otherCat = otherCat.view(float).reshape(-1, 2)

    ## Returns a 2-tuple of dx, dy
    center = (pyfits.getval(otherFile, 'CRPIX1'),
              pyfits.getval(otherFile, 'CRPIX2'))
    offset, offerr = CatTools.findOffsetMCMC(refCat, otherCat,
                                             rotOrigin=center,
                                             **SHIFT_MCMC_PARAMS)
    return offset


def write_shift_file(ref_file, file_list, shift_list, file_name='shifts.txt'):
    """
    Writes the multidrizzle shift file given a reference file name, list of
    files, and list of shifts. Shift file is written as delta shift file in
    pixel coordinates. Rotation is assumed 0 and scale assumed 1
    """
    ## Write out shift file for multidrizzle
    with open(file_name, 'w') as sf:
        sf.writelines(['# frame: output\n',
                       '# refimage: ' + ref_file + '\n',
                       '# form: delta\n',
                       '# units: pixels\n'])
        lines = []
        for fname, shift in zip(file_list, shift_list):
            lines += ['{}   {:0.2f}   {:0.2f}   {:0.2f}   1.0\n'.format(fname,
                                                                        *shift)]
        sf.writelines(lines)


def trim_files(field):
    """
    Trims files to a common size. Multidrizzle contains a bug:
    When drizzling with a reference image, the resultant image size is sometimes
    1 pixel too large. This messes up sextractor dual-image photometry.
    This function cuts off the offending line or column.
    @param field: directory to act upon
    """
    sciFiles = glob.glob('{0}/F*/{0}*_drz_sci.fits'.format(field))
    outX = min([pyfits.getval(sciFile, 'NAXIS1') for sciFile in sciFiles])
    outY = min([pyfits.getval(sciFile, 'NAXIS2') for sciFile in sciFiles])
    trimSlice = (slice(0, outY), slice(0, outX))

    for sciFile in sciFiles:
        wtFile = sciFile.replace('sci', 'weight')
        for filename in [sciFile, wtFile]:
            with pyfits.open(filename, mode='update') as f:
                f[0].data = f[0].data[trimSlice]
                f[0].update_header()
                f.flush()


def filter_list(fieldPath):
    """
    Returns the list of filter directories for a given field dir.
    """
    return [filt for filt in os.listdir(fieldPath)
            if re.match('^F\d', filt) and
            os.path.isdir(os.path.join(fieldPath, filt))]


def ref_filter(filtList):
    """
    Returns the name of the reference filter, given a list of filters
    (e.g. from filter_list, etc).
    """
    for filt in DRIZZLE_REF_FILTERS:
        if filt in filtList:
            return filt
    return ''


def parse_cl_args(argv):
    """
    Processes the command line arguments supplied to hippies_drizzle_all,
    replacing defaults, checking for file existence, etc.
    @param argv: Command-line arguments sent to hippies_drizzle_all
    @return: A dictionary of the command-line arguments
    """
    argv = [arg.strip('-') for arg in argv]
    _known_args = ['clean', 'help']
    args = dict([(arg, arg in argv) for arg in _known_args])
    return args


def is_processed(filt_path):
    """
    Checks to see if the supplied directory has already been drizzled
    """
    return len(glob.glob(os.path.join(filt_path, '*_drz_sci.fits'))) > 0


def clean_dir(filtDir, exclude=None):
    """
    Clean the supplied directory of multidrizzle products. Used for
    "clean" to first delete any previously-drizzled files
    """
    # selecting everything that is not a directory handles removing regular
    # files & symlinks
    oldFiles = [os.path.join(filtDir, filename)
                for filename in os.listdir(filtDir)
                if not os.path.isdir(os.path.join(filtDir, filename))]

    if exclude is not None:
        keepers = [os.path.join(filtDir, filename) for filename in exclude]
        oldFiles = set(oldFiles) - set(keepers)

    try:
        for filename in oldFiles:
            os.remove(filename)
    except OSError:
        pass
    message('Cleaned directory {}, keeping files {}'.format(filtDir, exclude))


def base_filename(filtPath):
    """
    Construct output base filename. Looks like parXXXX_W_FILT
    """
    pathElems = os.path.abspath(filtPath).split(os.path.sep)
    return pathElems[-2] + '_' + pathElems[-1] + '_drz.fits'


def instrument_info(fileName):
    """
    Gets the HST instrument info (Instrument, Detector) for a given file
    @return: Dictionary of instrument-specific values pulled from fits header
    """
    with pyfits.open(fileName) as fh:
        info = {'instrument': fh[0].header['INSTRUME'],
                'detector': fh[0].header['DETECTOR']}
    return info


if __name__ == '__main__':
    seterr(all='warn')
    ## Be helpful and go up a level if we run this from /unprocessed/
    if os.path.split(os.getcwd())[-1] == 'unprocessed':
        os.chdir('..')

    ## First, make sure there are some directories to run on.
    fields = [field for field in os.listdir(os.getcwd()) if
              field.startswith(FIELD_PREFIX)]

    args = parse_cl_args(sys.argv)

    sextractorConf = os.path.join(os.getcwd(), SEX_COMPACT_CONF)
    if len(fields) < 1 or args['help'] or not os.path.exists(sextractorConf):
        print(__doc__)
        exit(0)

    ## Ensure that reference file environment vars are set
    if 'iref' not in os.environ:
        message('iref environment variable not set. Please set it to enable ' +
                'WFC3IR drizzling.', msgtype='ERROR')
        exit(1)

    undrizzled = []
    for field in fields:
        drizzle_field(field, args)
        ## Sanity check: Make sure all filters were actually drizzled
        for filt in filter_list(field):
            filtPath = os.path.join(field, filt)
            sciName = base_filename(filtPath).replace('_drz', '_drz_sci')
            drizName = os.path.join(filtPath, sciName)
            if not os.path.exists(drizName):
                undrizzled += [filtPath]

    ## Output sanity check
    if len(undrizzled) > 0:
        message('Sanity check: The following filters were not drizzled! ' +
                'See Multidrizzle output for errors.', msgtype='WARN')
        message(', '.join(undrizzled), msgtype='WARN')

    message('Drizzling complete. To create RMS maps, run ' +
            '03_hippies_make_rms.py')
