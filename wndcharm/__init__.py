#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (C) 2015 National Institutes of Health 
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Written by:  Christopher Coletta (github.com/colettace)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__version__ = "unknown"

from wndcharm import *

try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass

try:
    from _git_hash import __git_hash__
    __version__ = __version__+ '+' + __git_hash__
except:
    pass

class _diagnostics( object ):
    """Report the versions of various Python packages WND-CHARM
    depends on/is often used with"""

    def __init__( self ):
        self.module_list = ['wndcharm', 'numpy', 'scipy', 'matplotlib', 'sklearn', \
                'skimage', 'IPython', 'tifffile', 'PIL', 'pandas']

    def get_package_versions( self ):
        """Runs through self.module_list, tries to import,
        then gets .__version__ or .VERSION"""

        ret = []
        import sys
        ret.append( ('python', sys.version ) )
        
        for name in self.module_list:
            m = None
            ver = None
            try: # 1. can we import it?
                m = __import__( name )
                try: #2. does it have a __version__?
                    ver = m.__version__
                except AttributeError:
                    try: # 3. Is it PIL which has a .VERSION instead?
                        ver = m.VERSION
                    except AttributeError:
                        ver = 'version not available'
            except ImportError:
                pass

            ret.append( ( name, ver ) )
        return ret

    def __call__( self ):
        return self.get_package_versions()

    def __str__( self ):
        outstr = "WND-CHARM Python API Diagnostics\n"
        outstr += "================================\n"

        from sys import executable
        outstr += "Executable:" + '\n\t' + str( executable ) + '\n'

        from os import getenv
        outstr += 'PYTHONPATH environment variable:\n\t' + \
                getenv( 'PYTHONPATH', '<unset>') + '\n'

        import wndcharm

        outstr += 'WND-CHARM library path:\n\t' + wndcharm.__file__ + '\n'
        
        outstr += 'WND-CHARM features major version: ' + str(feature_vector_major_version) + '\n'
        for type,foo in sorted(feature_vector_minor_version_from_vector_type.items(), key=lambda x: x[1]):
        	outstr += '\t'+type+': '
        	outstr += str(feature_vector_major_version) + '.' + str(feature_vector_minor_version_from_vector_type[type])
        	outstr += ' (N = '+str(feature_vector_num_features_from_vector_type[type]) + ')\n'
        

        outstr += 'Package versions:\n'
        retval = self.get_package_versions()

        for name, ver in retval:
            outstr += '\t' + str( name ).ljust(10) + '\t' + str( ver ).replace( '\n', ' ') + '\n'
        return outstr

diagnostics = _diagnostics()

# The numbers *must* be consistent with what's defined in wndchrm C-codebase.
feature_vector_major_version = wndcharm.StdFeatureComputationPlans.feature_vector_major_version
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
feature_vector_num_features_from_vector_version = {
    '1.1':1025,
    '1.2':2873,
    '1.3':2160,
    '1.4':4008,
    '2.1':1059,
    '2.2':2919,
    '2.3':2199,
    '2.4':4059,
    '3.1':1059,
    '3.2':2919,
    '3.3':2199,
    '3.4':4059,
    '4.1':1059,
    '4.2':2919,
    '4.3':2199,
    '4.4':4059,
}
