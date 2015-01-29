__version__ = "unknown"

from wndcharm import *

try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass

try:
    from _svn_version import __svn_version__
    __version__ = __version__+'-'+__svn_version__
except ImportError:
    # We're running in a tree that doesn't have a _svn_version.py, so we don't know what our version is.
    pass


# The numbers *must* be consistent with what's defined in wndchrm C-codebase.
feature_vector_major_version = 2
# Feature vector lengths in current version
# #define NUM_LC_FEATURES  4059
# #define NUM_L_FEATURES   2919
# #define NUM_C_FEATURES   2199
# #define NUM_DEF_FEATURES 1059
# These are definitions for Version 2 features.
feature_vector_minor_version_from_num_features = {
    1059:1,
    2919:2,
    2199:3,
    4059:4
}
# // original lengths prior to Version 2:
# // no Gini coefficient, no inverse otsu features
# // #define NUM_LC_FEATURES  4008
# // #define NUM_L_FEATURES   2873
# // #define NUM_C_FEATURES   2160
# // #define NUM_DEF_FEATURES 1025
feature_vector_minor_version_from_num_features_v1 = {
    1025:1,
    2873:2,
    2160:3,
    4008:4
}
feature_vector_minor_version_from_vector_type = {
    'short':1,
    'long':2,
    'short_color':3,
    'long_color':4
}
feature_vector_num_features_from_vector_type = {
    'short':1059,
    'long':2919,
    'short_color':2199,
    'long_color':4059
}
