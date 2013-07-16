__version__ = "unknown"

__all__ = [ 'pychrm', 'FeatureSet' ]

from pychrm import *

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



