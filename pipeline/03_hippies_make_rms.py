"""
Creates RMS maps from Multidrizzle weight maps by measuring the autocorrelation
in areas of blank sky.

Available command-line options:
clean            Clean existing files first? Default is no.
rmsboxsize NN    Size (in pixels) of the boxes used to measure autocorrelation
                 Default is 256
numrmsboxes NN   Number of boxes to average when calculating autocorrelation
                 Default is 10
"""

from __future__ import division
import sys
import os
import numpy as np
import astroRMS
import pyfits
from scipy.ndimage import uniform_filter, minimum_filter, maximum_filter

_field_prefix = 'par'
## Name of the region file for calculating RMS within each field dir
_region_name = 'rms_calc.reg'
_rms_bad_px_val = 10000


def parse_cl_args(argv):
    """
    Processes the command line arguments supplied to hippies_phot,
    replacing defaults, checking for file existence, etc.

    @param argv: Command-line arguments sent to hippies_phot
    @return: A dictionary of the command-line arguments
    """
    argv = [arg.strip('-') for arg in argv]
    _known_args = ['help', 'overwrite']
    _known_val_args = ['rmsboxsize', 'numrmsboxes']
    args = dict([(arg, arg in argv) for arg in _known_args])
    for arg in _known_val_args:
        if arg in argv:
            keypos = argv.index(arg)
            args[arg] = argv[keypos+1]
    ## Set defaults and do additional per-argument processing
    args.setdefault('rmsboxsize', 256)
    args.setdefault('numrmsboxes', 10)
    args['rmsboxsize'] = int(args['rmsboxsize'])
    args['numrmsboxes'] = int(args['numrmsboxes'])
    return args


def selectRMSBoxSlices(sciFile, weightFiles, boxSize=256, numBoxes=10):
    ## Ensure weightFiles is a list
    if isinstance(weightFiles, basestring):
        weightFiles = [weightFiles]

    with pyfits.open(sciFile) as f:
        imgBoxcar = uniform_filter(f[0].data, size=boxSize)

    ## Allocate our mask
    mask = np.zeros(imgBoxcar.shape, dtype=bool)
    ## 'Or' combine masks of zero-weight pixels from all weight files
    for weightFile in weightFiles:
        with pyfits.open(weightFile) as w:
            mask |= (w[0].data <= 0)

    # Smooth over mask image with min filter, to ignore small areas of bad
    # pixels
    smoothSize = boxSize // 10
    mask = minimum_filter(mask, size=smoothSize)
    ## Expand zone of avoidance of bad pixels, so we don't pick boxes that
    ## contain them. mode=constant, cval=True means treat all borders as
    ## if they were masked-out pixels
    mask = maximum_filter(mask, size=(boxSize + smoothSize), mode='constant',
                          cval=True)

    imgBoxcar = np.ma.array(imgBoxcar, mask=mask)

    boxlocs = []
    for box in xrange(numBoxes):
        ## Find the location of the minimum value of the boxcar image,
        ## excluding masked areas. This will be a pixel with few nearby sources
        ## within one box width
        minloc = imgBoxcar.argmin()
        minloc = np.unravel_index(minloc, imgBoxcar.shape)
        lowerLeft = tuple(int(x - boxSize / 2) for x in minloc)
        ## Negative values of lowerLeft mean argmin ran out of unmasked pixels
        if lowerLeft[0] < 0 or lowerLeft[1] < 0:
            message('Ran out of good pixels when placing RMS calculation ' +
                    'regions for file {}.'.format(sciFile), msgtype='WARN')
            lowerLeft = tuple(slc.start + smoothSize for slc in boxlocs[-1])
        minSlice = tuple(slice(x, x + boxSize) for x in lowerLeft)
        boxlocs += [minSlice]

        ## Zone of avoidance is twice as big, since we are picking box centers
        ## Use clip to ensure avoidance slice stays within array bounds
        lowerLeft = tuple(int(x - boxSize) for x in minloc)
        avoidSlice = tuple(slice(np.clip(x, 0, extent),
                                 np.clip(x + 2 * boxSize, 0, extent))
                           for x, extent in zip(lowerLeft, imgBoxcar.shape))

        ## Add this box to the mask
        imgBoxcar[avoidSlice] = np.ma.masked

    return boxlocs


def createRMSMap(field, filt, calcSlices):
    ## Open drizzled image and weight image
    fNames = getFileNamesDict(field, filt)
    sciData = pyfits.getdata(os.path.join(field, fNames['sci']))
    weightImage = pyfits.open(os.path.join(field, fNames['weight']))
    weightData = weightImage[0].data

    weightToVarScale = 0.0
    ## Calculate weightmap to rmsmap scaling for each slice, average
    message('Autocorrelating regions in {}/{}:'.format(field, filt))
    for sliceY, sliceX in calcSlices:
        autoCorrDict = astroRMS.calc_RMS(sciData[sliceY, sliceX],
                                         weightData[sliceY, sliceX])
        weightToVarScale += autoCorrDict['weightScale'] / len(calcSlices)
        sliceStr = '[{:d}:{:d},{:d}:{:d}]'.format(sliceX.start, sliceX.stop,
                                                  sliceY.start, sliceY.stop)
        message('RMS scale from region {}: {:.5f}'.format(
            sliceStr, autoCorrDict['weightScale']))

    ## Bad pixels are those with 0 or negative weight.
    bpMask = weightImage[0].data <= 0
    # Now transform the weight map into an RMS map using weightToVarScale.
    # Set bad pixels to very large RMS
    rmsData = np.divide(1.0, np.sqrt(weightImage[0].data * weightToVarScale))
    rmsData[bpMask] = _rms_bad_px_val

    sciHeader = pyfits.getheader(os.path.join(field, fNames['sci']))

    ## Save files
    pyfits.writeto(os.path.join(field, fNames['rms']), header=sciHeader,
                   data=rmsData, clobber=True)
    weightImage.close()

    statSection = rmsData[~bpMask]
    rmsStats = {
        'image': fNames['rms'],
        'mean': np.mean(statSection),
        'stddev': np.std(statSection),
        'min': np.min(statSection),
        'max': np.max(statSection),
        'numPix': statSection.size
    }

    return rmsStats


def writeSlicesToRegionFile(filename, slices):
    regFile = open(filename, 'w')
    regFile.write('image\n')
    for sliceY, sliceX in slices:
        regFile.write('box({:.2f},{:.2f},{:d},{:d},0)\n'.format(
            0.5 * (sliceX.stop + sliceX.start - 1),
            0.5 * (sliceY.stop + sliceY.start - 1),
            sliceX.stop - sliceX.start,
            sliceY.stop - sliceY.start))
    regFile.close()


def readSlicesFromRegionFile(fileName):
    retSlices = []
    regFile = open(fileName)
    for line in regFile:
        ## Get just the arguments passed to box()
        boxArgs = line.strip('box()\n ').replace(',', ' ').split()
        try:
            ## Convert arguments to numeric types
            boxArgs = [float(arg) for arg in boxArgs[0:2]] + [int(arg) for arg
                                                              in boxArgs[2:]]
        except ValueError:
            ## Just skip over malformed lines: comments, circles, etc.
            continue
        ## Calculate start and stop as in writeSlicesToRegionFile()
        llCorner = [(int(boxArgs[0] * 2) - boxArgs[2] + 1) // 2,
                    (int(boxArgs[1] * 2) - boxArgs[3] + 1) // 2]
        urCorner = [llCorner[0] + boxArgs[2], llCorner[1] + boxArgs[3]]
        ## Convert to slices for numpy, recalling numpy arrays are Y,X
        retSlices += [
            (slice(llCorner[1], urCorner[1]), slice(llCorner[0], urCorner[0]))]
    regFile.close()
    return retSlices


def linkDrizzledFilters(field):
    os.chdir(field)
    filters = [path for path in os.listdir('.') if path.startswith('F')]
    linkedFilters = list(filters)  # Copy list
    for filt in filters:
        fNames = getFileNamesDict(field, filt)
        try:
            if os.path.exists(os.path.join(filt, fNames['sci'])):
                os.symlink(os.path.join(filt, fNames['sci']), fNames['sci'])
            if os.path.exists(os.path.join(filt, fNames['weight'])):
                os.symlink(os.path.join(filt, fNames['weight']),
                           fNames['weight'])
        ## Simply ignore error for the symlink already existing
        except OSError:
            pass
        if not os.path.exists(fNames['sci']) or not os.path.exists(
                fNames['weight']):
            linkedFilters.remove(filt)
    os.chdir('..')
    return linkedFilters


def getFileNamesDict(field, filt):
    fileTypes = ('sci', 'weight', 'rms', 'mskrms')
    return dict(
        [(ftype, '_'.join([field, filt, 'drz', '{}.fits'.format(ftype)]))
         for ftype in fileTypes])


def message(msg, msgtype='INFO'):
    print('[{}] {}'.format(msgtype, msg))


if __name__ == '__main__':
    np.seterr(divide='warn')
    ## Be helpful and go up a level if we run this from /unprocessed/
    if os.path.split(os.getcwd())[-1] == 'unprocessed':
        os.chdir('..')

    ## Get command-line arguments
    args = parse_cl_args(sys.argv[1:])

    ## We'll find all the fields first
    fields = [field for field in os.listdir(os.getcwd()) if
              field.startswith(_field_prefix)]

    for field in fields:
        # First task is to symlink all drizzled sci and weight files into the
        # par_XXX_W dirs
        filters = linkDrizzledFilters(field)

        ## Find optimal regions to calculate RMS scaling from
        if not os.path.exists(os.path.join(field, _region_name)):
            allNames = [getFileNamesDict(field, filt) for filt in filters]
            weightPaths = [os.path.join(field, fNames['weight']) for fNames in
                           allNames]
            rmsCalcSlices = selectRMSBoxSlices(
                os.path.join(field, allNames[0]['sci']),
                weightPaths,
                args['rmsboxsize'],
                args['numrmsboxes'])

            ## Save out the RMS calculations region file
            writeSlicesToRegionFile(os.path.join(field, _region_name),
                                    rmsCalcSlices)

        else:
            rmsCalcSlices = readSlicesFromRegionFile(
                os.path.join(field, _region_name))

        ## Now we need to derive RMS maps for each filter in each field.
        for filt in filters:
            fNames = getFileNamesDict(field, filt)
            if (os.path.exists(os.path.join(field, fNames['rms']))
                and not args['overwrite']):
                message('RMS mask exists: {}. '.format(fNames['rms']) +
                        'To re-calculate, run with overwrite.',
                        msgtype='WARN')
                continue

            rmsStats = createRMSMap(field, filt, rmsCalcSlices)

    message('RMS map creation complete. To run photometry on all fields, ' +
            'run 04_hippies_phot.py')