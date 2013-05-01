"""
HIPPIES pipeline configuration variables. All user-configurable parameters can
be changed in this file.
"""

# Prefix for field directories. Field RA/Dec are appended to this.
FIELD_PREFIX = 'par'

# In degrees, tolerance for connecting fields. Right now it's 90 arcsec, about
# half the WFC3IR field of view. Don't make this too big, or the drizzling
# process will have trouble matching catalogs.
FIELD_CONNECT_DISTANCE = 90. / 3600.

# Multidrizzle parameters for final image combination
DRIZ_PARAMS_WFC3 = {'driz_sep_scale': 0.1, 'driz_sep_rot': 0.0,
                    'driz_final_scale': 0.1, 'driz_final_rot': 0.0,
                    'driz_final_pixfrac': 0.8, 'driz_final_wht_type': 'EXP',
                    'driz_final_kernel': 'square', 'build': False}

# List of reference filters, in order of preference. Every field should have
# one of these filters, and deepest IR filters should be first. Filters are
# checked in order, so if F105W is not found, F098M will be tried, and so on.
DRIZZLE_REF_FILTERS = ['F105W', 'F098M']

# Sextractor shell command. Can replace with /usr/local/bin/sex etc if you want
# to use a specific sextractor binary
SEX_COMMAND = 'sex'

# Sextractor configuration file used for finding compact objects, during the
# drizzle exposure alignment step.
SEX_COMPACT_CONF = 'sextractor/getcompact.conf'

# Sextractor configuration file used for photometry
SEX_PHOTOMETRY_CONF = 'sextractor/photometry.conf'

# Name of the region file for calculating RMS within each field dir
RMS_REGION_NAME = 'rms_calc.reg'

# Value to assign to bad (zero-weight) pixels when creating RMS maps
RMS_BAD_PX_VALUE = 10000

# Photometry detection filters, i.e. detection image for dual-image mode
# sextractor. Filters are checked in order, so if F160W is not found, F125W will
# be tried, and so on.
PHOT_DETECT_FILTERS = ['F160W', 'F125W']

# Filter Photometric zeropoints, in AB magnitudes. From:
# http://www.stsci.edu/hst/wfc3/phot_zp_lbn
FILTER_ZEROPOINTS = {'F300X': 24.9565,
                     'F475X': 26.1579,
                     'F600LP': 25.8746,
                     'F850LP': 23.8338,
                     'F606W': 26.0691,
                     'F098M': 25.6674,
                     'F105W': 26.2687,
                     'F110W': 26.8223,
                     'F125W': 26.2303,
                     'F160W': 25.9463}

# Sometimes columns are erroneously duplicated in the master catalogs due to
# floating point roundoff errors. They can be manually specified here to force
# them not to be duplicated in the final catalogs
IGNORE_DUPLICATE_COLUMNS = ['ALPHA_J2000', 'DELTA_J2000']


# This function is used by all scripts to print messages to the terminal.
def message(msg, msgtype='INFO'):
    print('\n[{}] {}\n'.format(msgtype, msg))
