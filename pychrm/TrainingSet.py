#!/usr/bin/env python

try:
	from . import pychrm
except:
	import pychrm

# FeatureRegistration.py is where the SWIG wrapped objects get put into a dict
# for use in signature calculation
try:
	from . import FeatureRegistration 
except:
	import FeatureRegistration

# FeatureNameMap.py contains mapping from old style names to new style
# and the function TranslateFeatureNames()
try:
	from . import FeatureNameMap
except:
	import FeatureNameMap

import numpy as np
from StringIO import StringIO
try:
	import cPickle as pickle
except:
	import pickle

# os module has function os.walk( ... )
import os
import os.path 
import re
import itertools
import copy


# Initialize module level globals
Algorithms = None
Transforms = None
small_featureset_featuregroup_strings = None
large_featureset_featuregroup_strings = None
small_featureset_featuregroup_list = None
large_featureset_featuregroup_list = None

error_banner = "\n*************************************************************************\n"

def initialize_module(): 
	# If you're going to calculate any signatures, you need this stuff
	# FIXME: Maybe rig something up using a load on demand?
	global Algorithms
	global Transforms
	global small_featureset_featuregroup_strings
	global large_featureset_featuregroup_strings
	global small_featureset_featuregroup_list
	global large_featureset_featuregroup_list

	Algorithms = FeatureRegistration.LoadFeatureAlgorithms()
	Transforms = FeatureRegistration.LoadFeatureTransforms()

	feature_lists = FeatureRegistration.LoadSmallAndLargeFeatureSetStringLists()

	small_featureset_featuregroup_strings = feature_lists[0]
	full_list = "\n"
	large_featureset_featuregroup_strings = full_list.join( feature_lists )

	small_featureset_featuregroup_list = []
	for fg_str in small_featureset_featuregroup_strings.splitlines():
		fg = ParseFeatureGroupString( fg_str )
		if fg:
			small_featureset_featuregroup_list.append( fg )

	large_featureset_featuregroup_list = []
	for fg_str in large_featureset_featuregroup_strings.splitlines():
		fg = ParseFeatureGroupString( fg_str )
		if fg:
			large_featureset_featuregroup_list.append( fg )

#############################################################################
# class definition of FeatureVector
#############################################################################
# FeatureVector inherits from class object and is thus a Python "new-style" class
class FeatureVector(object):
	"""
	FIXME: Would be nice to define a [] operator that returns
	       a tuple = (self.names[n], self.values[n])
	"""
	names = None
	values = None
	#================================================================
	def __init__( self, data_dict = None):
		"""@brief: constructor
		"""
		self.names = []
		self.values = []

		if data_dict:
			if "values" in data_dict:
				self.values = data_dict[ "values" ]
			if "names" in data_dict:
				self.names = data_dict[ "names" ]
	#================================================================
	def is_valid( self ):
		"""
		@brief: an instance should know all the criteria for being a valid FeatureVector
		"""
		if len( self.values ) != len( self.names ):
			raise RuntimeError( "Instance of {0} is invalid: ".format( self.__class__ ) + \
			  "different number of values ({0}) and names ({1}).".format( \
			  len( self.feature_values ), len( self.feature_names ) ) )
		return True

#############################################################################
# class definition of FeatureWeights
#############################################################################
class FeatureWeights( FeatureVector ):
	"""
	"""
	def __init__( self, data_dict = None ):
		# call parent constructor
		super( FeatureWeights, self ).__init__( data_dict )

	#================================================================
	@classmethod
	def NewFromFile( cls, weights_filepath ):
		"""@brief written to read in files created by wndchrm -vw/path/to/weightsfile.txt"""

		weights = cls()
		with open( weights_filepath, 'r' ) as weights_file:
			for line in weights_file:
				# split line "number <space> name"
				feature_line = line.strip().split( " ", 1 )
				weights.values.append( float( feature_line[0] ) )
				weights.names.append( feature_line[1] )
		return weights

	#================================================================
	def EliminateZeros( self ):
		"""@breif Good for Fisher scores, N/A for Pearson scores - FIXME: subclass!"""

		new_weights = FeatureWeights()
		scores = zip( self.names, self.values )
		nonzero_scores = [ (name, weight) for name, weight in scores if weight != 0 ]
		new_weights.names, new_weights.values = zip( *nonzero_scores )
		return new_weights

	#================================================================
	def Threshold( self, num_features_to_be_used  ):
		"""@breif Good for Fisher scores, N/A for Pearson scores - FIXME: subclass!"""

		new_weights = FeatureWeights()
		raw_featureweights = zip( self.names, self.values )
		# raw_featureweights is now a list of tuples := [ (name1, value1), (name2, value2), ... ]

		# sort from max to min
		# sort by the second item in the tuple, i.e., index 1
		sort_func = lambda feat_a, feat_b: cmp( feat_a[1], feat_b[1] ) 

		sorted_featureweights = sorted( raw_featureweights, sort_func, reverse = True )
		# take top N features
		use_these_feature_weights = \
				list( itertools.islice( sorted_featureweights, num_features_to_be_used ) )
		
		# we want lists, not tuples!
		new_weights.names, new_weights.values =\
		  [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]
		return new_weights

	#================================================================
	def PickleMe( self, pathname = None ):

		outfile_pathname = ""
		if pathname != None:
			outfile_pathname = pathname
		else:
			outfile_pathname = "feature_weights_len_{0}.weights.pickled".format( len( self.names ) )

		if os.path.exists( outfile_pathname ):
			print "Overwriting {0}".format( outfile_pathname )
		else:
			print "Writing {0}".format( outfile_pathname )
		with open( outfile_pathname, 'wb') as outfile:
			pickle.dump( self.__dict__, outfile, pickle.HIGHEST_PROTOCOL )

  #=================================================================================
	@classmethod
	def NewFromPickleFile( cls, pathname ):
		"""
		The pickle is in the form of a dict
		"""
		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".weights.pickled" ):
			raise ValueError( "File isn't named with .weights.pickled extension: {0}".format( pathname ) )

		print "Loading feature weights from pickled file {0}".format( pathname )
		weights = None
		with open( pathname, "rb" ) as pkled_in:
			weights = cls( pickle.load( pkled_in ) )

		return weights

	

#############################################################################
# class definition of Signatures
#############################################################################
class Signatures( FeatureVector ):
	"""
	"""

	path_to_source_image = None
	options = ""

	#================================================================
	def __init__( self ):
		"""@brief: constructor
		"""
		# call parent constructor
		super( Signatures, self ).__init__()

	#================================================================
	@classmethod
	def SmallFeatureSet( cls, imagepath, options = None ):
		"""@brief Equivalent of invoking wndchrm train in c-chrm
		@argument path - path to a tiff file
		@return An instance of the class Signatures for image with sigs calculated."""

		print "====================================================================="
		print "Calculating small feature set for file:"
		global small_featureset_featuregroup_list
		return cls.NewFromFeatureGroupList( imagepath, small_featureset_featuregroup_list, options )

	#================================================================
	@classmethod
	def LargeFeatureSet( cls, imagepath, options = None ):
		"""@brief Equivalent of invoking wndchrm train -l in c-chrm
		@argument path - path to a tiff file
		@return An instance of the class Signatures for image with sigs calculated."""

		print "====================================================================="
		print "Calculating large feature set for file:"

		the_sigs = cls()

		# Try to find a file that corresponds to this image with signatures in it.
		# derive the name of the sig/pysig file that would exist if there were one
		# given the options requested

		sigpath = None
		root, extension = os.path.splitext( imagepath )
		options_str = options if options else ""
		if not "-l" in options_str:
			options_str += "-l"

		if os.path.exists( root + options_str + ".pysig" ):
			sigpath = root + options_str + ".pysig"
			the_sigs = cls.FromSigFile( imagepath, sigpath, options_str )
		elif os.path.exists( root + options_str + ".sig" ):
			sigpath = root + options_str + ".sig" 
			the_sigs = cls.FromSigFile( imagepath, sigpath, options_str )
		else:
			sigpath = root + options_str + ".pysig"


		# Check to see what, if anything was loaded from file. The file could be corrupted
		# incomplete, or calculated with different options, e.g., -S1441
		# FIXME: Here's where you'd calculate a small subset of features
		# and see if they match what was loaded from file
		if len( the_sigs.names ) <= 0:
			# All hope is lost. Calculate sigs
			global large_featureset_featuregroup_list
			if options == None:
				the_sigs = cls.NewFromFeatureGroupList( imagepath, large_featureset_featuregroup_list, "-l" )
			else:
				# FIXME: dummyproofing: does options already contain '-l'?
				retval = cls.NewFromFeatureGroupList( imagepath, large_featureset_featuregroup_list, \
						options + "-l" )
		# There are historically 2885 signatures in the large feature set 
		elif len( the_sigs.names ) < 2873:
			# FIXME: Calculate the difference of features between what was loaded
			# from file, and what is supposed to be present in large feature set
			# and calculate only those. They might also need to be placed in an order
			# that makes sense.
			err_msg = "There are only {0} signatures in file '{1}'".format( \
			  len( the_sigs.names ), sigpath )
			err_msg += "where there are supposed to be at least 2873."
			raise ValueError( err_msg )
		return the_sigs


	#================================================================
	@classmethod
	def NewFromFeatureGroupList( cls, path_to_image, feature_groups, options = None ):
		"""@brief calculates signatures

		"""

		print path_to_image
		original = pychrm.ImageMatrix()
		if 1 != original.OpenImage( path_to_image, 0, None, 0, 0 ):
			raise ValueError( 'Could not build an ImageMatrix from {0}, check the path.'.\
			    format( path_to_image ) )

		im_cache = {}
		im_cache[ '' ] = original

		# instantiate an empty Signatures object
		signatures = cls()
		signatures.path_to_source_image = path_to_image
		signatures.options = options

		for fg in feature_groups:
			print "Group {0}".format( fg.name )
			returned_feature_vals = fg.CalculateFeatures( im_cache )
			count = 0
			for value in returned_feature_vals:
				signatures.names.append( fg.name + " [{0}]".format( count ) )
				signatures.values.append( value )	
				count += 1

		return signatures

	#================================================================
	@classmethod
	def NewFromFeatureNameList( cls, path_to_image, feature_names, options = None ):
		"""@brief calculates signatures

		"""

		work_order, num_features = GenerateWorkOrderFromListOfFeatureStrings( feature_names )
		sig = cls.NewFromFeatureGroupList( path_to_image, work_order, options )
		return sig.FeatureReduce( feature_names ) 

	#================================================================
	@classmethod
	def FromPickle( cls, path ):
		"""
		FIXME: Implement!
		"""
		raise NotImplementedError()

	#================================================================
	def PickleMe( self ):
		"""
		FIXME: Implement!
		"""
		raise NotImplementedError()

	#================================================================
	@classmethod
	def FromSigFile( cls, image_path, sigfile_path, options = None ):
		"""@argument sigfile_path must be a .sig or a .pysig file
		
		@return  - An instantiated signature class with feature names translated from
		           the old naming convention, if applicable.
		@remarks - old style sig files don't know about their options other than 
		           the information contained in the name. In the future pysig files
		           may keep that info within the file. Thus, for now, the argument
		           options is something set externally rather than gleaned from
		           reading the file."""

		print "Loading features from sigfile {0}".format( sigfile_path )

		signatures = cls()
		signatures.path_to_source_image = image_path
		signatures.options = options
 
		with open( sigfile_path ) as infile:
			linenum = 0
			for line in infile:
				if linenum == 0:
					# The class id here may be trash
					signatures.class_id = line.strip()
				elif linenum == 1:
					# We've already assigned self.path_to_source_image
					# the path in the sig file may be trash anyway
					#signatures.path_to_source_image = line.strip()
					pass
				else:
					value, name = line.strip().split( ' ', 1 )
					signatures.values.append( float( value ) )	
					signatures.names.append( name )
				linenum += 1
			print "Loaded {0} signatures.".format( len( signatures.values ) )

		# Check if the feature name follows the old convention
		signatures.names = FeatureNameMap.TranslateToNewStyle( signatures.names ) 
		
		return signatures
	
	#================================================================
	def WriteFeaturesToASCIISigFile( self, filepath = None ):
		"""Write a sig file.
		
		If filepath is specified, you get to name it whatever you want and put it
		wherever you want. Otherwise, it's named according to convention and placed 
		next to the image file in its directory."""

		self.is_valid()

		outfile_path = ""
		if not filepath or filepath == "":
			if not self.path_to_source_image or self.path_to_source_image == "":
				raise ValueError( "Can't write sig file. No filepath specified in function call, and no path associated with this instance of Signatures." )
			outfile_path = self.path_to_source_image

			path, filename = os.path.split( outfile_path )
			if not os.path.exists( path ):
				raise ValueError( 'Invalid path {0}'.format( path ) )

			filename_parts = filename.rsplit( '.', 1 )
			if self.options and self.options is not "":
				outfile_path = "{0}{1}.pysig".format( filename_parts[0],\
																					self.options if self.options else "" )
			else:
				outfile_path = "{0}.pysig".format( filename_parts[0] )
			outfile_path = os.path.join( path, outfile_path )

		if os.path.exists( outfile_path ):
			print "Overwriting {0}".format( outfile_path )
		else:
			print 'Writing signature file "{0}"'.format( outfile_path )

		with open( outfile_path, "w" ) as out_file:
			# FIXME: line 2 contains class membership, just hardcode a number for now
			out_file.write( "0\n" )
			out_file.write( "{0}\n".format( self.path_to_source_image ) )
			for i in range( 0, len( self.names ) ):
				out_file.write( "{val:0.6f} {name}\n".format( val=self.values[i], name=self.names[i] ) )

	#================================================================
	def FeatureReduce( self, requested_features ):

		if self.names == requested_features:
			return self

		selfs_features = set( self.names )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = error_banner + "Feature Reduction error:\n"
			err_str += "The signatures set for image file '{0}' ".format( self.path_to_source_image )
			err_str += "is missing {0}".format( len( missing_features_from_req ) )
			err_str += "/{1} features that were requested in the feature reduction list.".format(\
			             len( requested_features ) )
			raise ValueError( err_str )

		print "Performing feature reduce/reorder on signatures..."
		# The featurenames sets are either out of order or contain extraneous features.
		dictionary = dict( zip( self.names, self.values ) )
		reduced_sigs = Signatures()
		reduced_sigs.path_to_source_image = self.path_to_source_image
		reduced_sigs.options = self.options
			
		for new_name in requested_features:
			reduced_sigs.names.append( new_name )
			reduced_sigs.values.append( dictionary[ new_name ] )
		
		return reduced_sigs

	def Normalize( self, training_set ):
		"""
		FIXME: should there be some flag that gets set if this sig has 
		already been normalized??
		
		@return: a newly instantiated Signatures
		"""

		if training_set.featurenames_list != self.names:
			raise ValueError("Can't normalize signature for {0} against training_set {1}: Features don't match."\
		  .format( self.path_to_source_image, training_set.source_path ) )
	
		new_sig = copy.deepcopy( self )
		for i in range( len( self.values ) ):
			new_sig.values[i] -= training_set.feature_minima[i]
			new_sig.values[i] /= (training_set.feature_maxima[i] -training_set.feature_minima[i])
			new_sig.values[i] *= 100

		return new_sig
		

# end definition class Signatures

#############################################################################
# class definition of FeatureGroup
#############################################################################
class FeatureGroup:
	"""
	Attributes Name, Alg and Tforms are references to the SWIG objects
	"""

	name = ""
	algorithm = None
	transform_list = []
	def __init__( self, name_str = "", algorithm = None, tform_list = [] ):
		#print "Creating new FeatureGroup for string {0}:".format(name_str)
		#print "\talgorithm: {0}, transform list: {1}".format( algorithm, tform_list )
		self.name = name_str 
		self.algorithm = algorithm
		self.transform_list = tform_list
	def CalculateFeatures( self, cached_pixel_planes ):
		"""Returns a tuple with the features"""
		pixel_plane = None
		try:
			#print "transforms: {0}".format( self.Tforms )
			pixel_plane = RetrievePixelPlane( cached_pixel_planes, self.transform_list )
		except:
			raise
		return self.algorithm.calculate( pixel_plane )


#############################################################################
# global functions
#############################################################################
def RetrievePixelPlane( image_matrix_cache, tform_list ):
	"""
	Returns the image matrix prescribed in tform_list
	If it already exists in cache, just return.
	If it doesn't exist calculates it
	Recurses through the compound transform chain in tform_list
	"""
	#print "passed in: {0}".format( tform_list )
	requested_transform = " ".join( [ tform.name for tform in tform_list ] )
	#print "requesting pixel plane: '{0}'".format( requested_transform )
	if requested_transform in image_matrix_cache:
		return image_matrix_cache[ requested_transform ]
	
	# Required transform isn't in the cache, gotta make it
	# Pop transforms off the end sequentially and check to see if
	# lower-level transforms have already been calculated and stored in cache

	# Can't begin if there isn't at least the raw (untransformed) pixel plane
	# already stored in the cache
	if image_matrix_cache is None or len(image_matrix_cache) == 0:
		raise ValueError( "Can't calculate features: couldn't find the original pixel plane" +\
		                  "to calculate features {0}.".format( self.Name ) )

	sublist = tform_list[:]
	sublist.reverse()
	top_level_transform = sublist.pop()
	intermediate_pixel_plane = RetrievePixelPlane( image_matrix_cache, sublist )

	tformed_pp = top_level_transform.transform( intermediate_pixel_plane )
	#assert( intermediate_pixel_plane ), "Pixel Plane returned from transform() was NULL"
	image_matrix_cache[ requested_transform ] = tformed_pp
	return tformed_pp

#================================================================
def ParseFeatureGroupString( name ):
	"""Takes a string input, parses, and returns an instance of a FeatureGroup class"""
	#TBD: make a member function of the FeatureGroup
	# get the algorithm

	# The ability to comment out lines with a hashmark
	if name.startswith( '#' ):
		return None
	global Algorithms
	global Transforms
	string_rep = name.rstrip( ")" )
	parsed = string_rep.split( ' (' )
	
	alg = parsed[0]
	if alg not in Algorithms:
		raise KeyError( "Don't know about a feature algorithm with the name {0}".format(alg) )
	
	tform_list = parsed[1:]
	try:
		tform_list.remove( "" )
	except ValueError:
		pass
	if len(tform_list) != 0:
		for tform in tform_list:
			if tform not in Transforms:
				raise KeyError( "Don't know about a transform named {0}".format( tform ) )

	tform_swig_obj_list = [ Transforms[ tform ] for tform in tform_list ]

	return FeatureGroup( name, Algorithms[ alg ], tform_swig_obj_list )

#================================================================
def GenerateWorkOrderFromListOfFeatureStrings( feature_list ):
	"""
	Takes list of feature strings and chops off bin number at the first space on right, e.g.,
	"feature alg (transform()) [bin]" ... Returns a list of FeatureGroups.

	WARNING, RETURNS A TUPLE OF TWO THINGS!

	@return work_order - list of FeatureGroup objects
	@return output_features_count - total number of individual features contained in work_order
	"""

	feature_group_strings = set()
	output_features_count = 0

	for feature in feature_list:
		split_line = feature.rsplit( " ", 1 )
		# add to set to ensure uniqueness
		feature_group_strings.add( split_line[0] )

	# iterate over set and construct feature groups
	work_order = []
	for feature_group in feature_group_strings:
		fg = ParseFeatureGroupString( feature_group )
		output_features_count += fg.algorithm.n_features
		work_order.append( fg )

	return work_order, output_features_count

#############################################################################
# class definition of TrainingSet
#############################################################################
class TrainingSet:
	"""
  """

	# source_path - could be the name of a .fit, or pickle file from which this
	# instance was generated, could be a directory
	#  source_path is essentially a name
	# might want to make separate name member in the future
	source_path = ""
	num_classes = -1
	num_features = -1
	num_images = -1

	# For C classes, each with Ni images and M features:
	# If the dataset is contiguous, C = 1

	# A list of numpy matrices, length C (one Ni x M matrix for each class)
	# The design is such because it's useful to be able to quickly collect feature statistics
	# across an image class excluding the other classes
	data_list = None

	#: A list of strings, length C
	classnames_list = None

	#: A list of strings length M
	featurenames_list = None

	#: a list of lists, length C, where each list is length Ni, contining pathnames of tiles/imgs
	imagenames_list = None

	# The following class members are optional:
	# normalized_against is a string that keeps track of whether or not self has been
	# normalized. For test sets, value will be the source_path of the training_set.
	# For training sets, value will be "itself"
	normalized_against = None

	# Stored feature maxima and minima go in here
	# only exist, if self has been normalized against itself
	feature_maxima = None
	feature_minima = None

	# A list of floats against which marg probs can be multiplied
	# to obtain an interpolated value
	interpolation_coefficients = None

	# keep track of all the options (-l -S###, etc)
	# FIXME: expand to have all options kept track of individually
	feature_options = None

	def __init__( self, data_dict = None):
		"""
		TrainingSet constructor
		"""

		self.data_list = []
		self.classnames_list = []
		self.featurenames_list = []
		self.imagenames_list = []

		if data_dict != None:
			if "source_path" in data_dict:
				self.source_path = data_dict[ 'source_path' ]
			if "num_classes" in data_dict:
				self.num_classes = data_dict[ 'num_classes' ]
			if "num_features" in data_dict:
				self.num_features = data_dict[ 'num_features' ]
			if "num_images" in data_dict:
				self.num_images = data_dict[ 'num_images' ]
			if "data_list" in data_dict:
				self.data_list = data_dict[ 'data_list' ]
			if "classnames_list" in data_dict:
				self.classnames_list = data_dict[ 'classnames_list' ]
			if "featurenames_list" in data_dict:
				self.featurenames_list = data_dict[ 'featurenames_list' ]
			if "imagenames_list" in data_dict:
				self.imagenames_list = data_dict[ 'imagenames_list' ]
			if "feature_maxima" in data_dict:
				self.feature_maxima = data_dict[ 'feature_maxima' ]
			if "feature_minima" in data_dict:
				self.feature_minima = data_dict[ 'feature_minima' ]
			if "interpolation_coefficients" in data_dict:
				self.interpolation_coefficients = data_dict[ 'interpolation_coefficients' ]

  #=================================================================================
	@classmethod
	def NewFromPickleFile( cls, pathname ):
		"""
		The pickle is in the form of a dict
		FIXME: Shouldn't call Normalize if feature_maxima/minima are in the data_dict
		"""
		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".fit.pickled" ):
			raise ValueError( 'Not a pickled TrainingSet file: {0}'.format( pathname ) )

		print "Loading Training Set from pickled file {0}".format( pathname )
		the_training_set = None
		with open( pathname, "rb" ) as pkled_in:
			the_training_set = cls( pickle.load( pkled_in ) )

		# it might already be normalized!
		# FIXME: check for that
		# the_training_set.Normalize()

		return the_training_set

  #=================================================================================
	@classmethod
	def NewFromFitFile( cls, pathname ):
		"""
		Helper function which reads in a c-chrm fit file, builds a dict with the info
		Then calls the constructor and passes the dict as an argument
		"""
		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".fit" ):
			raise ValueError( 'Not a .fit file: {0}'.format( pathname ) )

		pickled_pathname = pathname + ".pychrm"

		print "Creating Training Set from legacy WND-CHARM text file file {0}".format( pathname )
		with open( pathname ) as fitfile:
			data_dict = {}
			data_dict[ 'source_path' ] = pathname
			data_dict[ 'data_list' ] = []
			data_dict[ 'imagenames_list' ] = []
			data_dict[ 'featurenames_list' ] = []
			data_dict[ 'classnames_list' ] = []
			data_dict[ 'imagenames_list' ] = []
			data_dict[ 'data_list' ] = []
			tmp_string_data_list = []

			name_line = False
			line_num = 0
			feature_count = 0
			image_pathname = ""
			num_classes = 0
			num_features = 0

			for line in fitfile:
				if line_num is 0:
					num_classes = int( line )
					data_dict[ 'num_classes' ] = num_classes
					# initialize list for string data
					for i in range( num_classes ):
						tmp_string_data_list.append( [] )
						data_dict[ 'imagenames_list' ].append( [] )
				elif line_num is 1:
					num_features = int( line )
					data_dict[ 'num_features' ] = num_features
				elif line_num is 2:
					data_dict[ 'num_images' ] = int( line )
				elif line_num <= ( num_features + 2 ):
					data_dict[ 'featurenames_list' ].append( line.strip() )
					feature_count += 1
				elif line_num == ( num_features + 3 ):
					pass # skip a line
				elif line_num <= ( num_features + 3 + num_classes ):
					data_dict[ 'classnames_list' ].append( line.strip() )
				else:
					# Read in features
					# Comes in alternating lines of data, then tile name
					if not name_line:
						# strip off the class identity value, which is the last in the array
						split_line = line.strip().rsplit( " ", 1)
						#print "class {0}".format( split_line[1] )
						zero_indexed_class_id = int( split_line[1] ) - 1
						tmp_string_data_list[ zero_indexed_class_id ].append( split_line[0] )
					else:
						image_pathname = line.strip()
						data_dict[ 'imagenames_list' ][ zero_indexed_class_id ].append( image_pathname )
					name_line = not name_line
				line_num += 1

		string_data = "\n"
		
		for i in range( num_classes ):
			print "generating matrix for class {0}".format( i )
			#print "{0}".format( tmp_string_data_list[i] )
			npmatr = np.genfromtxt( StringIO( string_data.join( tmp_string_data_list[i] ) ) )
			data_dict[ 'data_list' ].append( npmatr )

		# Can the class names be interpolated?
		tmp_vals = []
		for class_index in range( num_classes ):
			m = re.search( r'(\d*\.?\d+)', data_dict[ 'classnames_list' ][class_index] )
			if m:
				tmp_vals.append( float( m.group(1) ) )
			else:
				tmp_vals = None
				break
		if tmp_vals:
			data_dict[ 'interpolation_coefficients' ] = tmp_vals

		# Instantiate the class
		the_training_set = cls( data_dict )

		# normalize the features
		#the_training_set.Normalize()
		# no wait, don't normalize until we feature reduce!
		
		return the_training_set

  #=================================================================================
	@classmethod
	def NewFromSignature( cls, signature, ts_name = "TestSet", ):
		"""@brief Creates a new TrainingSet from a single signature
		Was written with performing a real-time classification in mind.
		"""

		try:
			signature.is_valid()
		except:
			raise

		new_ts = cls()
		new_ts.source_path = ts_name
		new_ts.num_classes = 1
		new_ts.num_features = len( signature.feature_values )
		new_ts.num_images = 1
		new_ts.classnames_list.append( "UNKNOWN" )
		new_ts.featurenames_list = signature.names
		new_ts.imagenames_list.append( [ inputimage_filepath ] )
		numpy_matrix = np.array( signature.values )
		new_ts.data_list.append( numpy_matrix )

		return new_ts

  #=================================================================================
	@classmethod
	def NewFromDirectory( cls, top_level_dir_path, feature_set = "large", write_sig_files_todisk = True ):
		"""
		@brief A quick and dirty implementation of the wndchrm train command
		Build up the self.imagenames_list, then pass it off to a sig classifier function
		"""
		print "Creating Training Set from directories of images {0}".format( top_level_dir_path )
		if not( os.path.exists( top_level_dir_path ) ):
			raise ValueError( 'Path "{0}" doesn\'t exist'.format( top_level_dir_path ) )
		if not( os.path.isdir( top_level_dir_path ) ):
			raise ValueError( 'Path "{0}" is not a directory'.format( top_level_dir_path ) )

		num_images = 0
		num_classes = 0
		classnames_list = []
		imagenames_list = []

		for root, dirs, files in os.walk( top_level_dir_path ):
			if root == top_level_dir_path:
				if len( dirs ) <= 0:
					# no class structure
					file_list = []
					for file in files:
						if '.tif' in file:
							file_list.append( os.path.join( root, file ) )
					if len( file_list ) <= 0:
						# nothing here to process!
						raise ValueError( 'No tiff files in directory {0}'.format( root ) )
					classnames_list.append( root )
					num_classes = 1
					num_images = len( file_list )
					imagenames_list.append( file_list )
					break
			else:
				file_list = []
				for file in files:
					if '.tif' in file:
						file_list.append( os.path.join( root, file ) )
				if len( file_list ) <= 0:
					# nothing here to process!
					continue
				# this class's name will be "subdir" in /path/to/topleveldir/subdir
				root, dirname = os.path.split( root )
				classnames_list.append( dirname )
				num_images += len( file_list )
				num_classes += 1
				imagenames_list.append( file_list )

		if num_classes <= 0:
			raise ValueError( 'No valid images or directories of images in this directory' )

		# instantiate a new training set
		new_ts = cls()
		new_ts.num_images = num_images
		new_ts.num_classes = num_classes
		new_ts.classnames_list = classnames_list
		new_ts.imagenames_list = imagenames_list
		new_ts.source_path = top_level_dir_path
		new_ts._ProcessSigCalculationSerially( feature_set, write_sig_files_todisk )
		if feature_set == "large":
			# FIXME: add other options
			new_ts.feature_options = "-l"
		return new_ts

  #=================================================================================
	@classmethod
	def NewFromFileOfFiles( cls, fof_path, feature_set = "large", write_sig_files_todisk = True ):
		"""FIXME: Implement!"""
		raise NotImplementedError()

  #=================================================================================
	@classmethod
	def NewFromSQLiteFile(cls, path):
		"""FIXME: Implement!"""
		raise NotImplementedError()

  #=================================================================================
	def _ProcessSigCalculationSerially( self, feature_set = "large", write_sig_files_to_disk = True, options = None ):
		"""
		Work off the self.imagenames_list
		"""
		# FIXME: check to see if any .sig, or .pysig files exist that match our
		#        Signature calculation criteria, and if so read them in and incorporate them

		sig = None
		class_id = 0
		for class_filelist in self.imagenames_list:
			for sourcefile in class_filelist:
				if feature_set == "large":
					sig = Signatures.LargeFeatureSet( sourcefile, options )
				elif feature_set == "small":
					sig = Signatures.SmallFeatureSet( sourcefile, options )
				else:
					raise ValueError( "sig calculation other than small and large feature set hasn't been implemented yet." )
				# FIXME: add all the other options
				# check validity
				if not sig:
					raise ValueError( "Couldn't create a valid signature from file {0} with options {1}".format( sourcefile, options ) )
				sig.is_valid()
				if write_sig_files_to_disk:
					sig.WriteFeaturesToASCIISigFile()
				self.AddSignature( sig, class_id )
			class_id += 1


  #=================================================================================
	def _ProcessSigCalculationParallelly( self, feature_set = "large", write_sig_files_todisk = True ):
		"""
		FIXME: When we figure out concurrency
		"""
		pass

  #=================================================================================
	def Normalize( self, training_set = None ):
		"""
		By convention, the range of values are normalized on an interval [0,100]
		FIXME: edge cases, clipping, etc
		"""

		if not( self.normalized_against ):

			full_stack = np.vstack( self.data_list )
			total_num_imgs, num_features = full_stack.shape
			self.feature_maxima = [None] * num_features
			self.feature_minima = [None] * num_features

			for i in range( num_features ):
				feature_max = np.max( full_stack[:,i] )
				self.feature_maxima[ i ] = feature_max
				feature_min = np.min( full_stack[:,i] )
				self.feature_minima[ i ] = feature_min
				for class_matrix in self.data_list:
					class_matrix[:,i] -= feature_min
					class_matrix[:,i] /= (feature_max - feature_min)
					class_matrix[:,i] *= 100
			self.normalized_against = "itself"

		if training_set:

			# sanity checks
			if self.normalized_against:
				raise ValueError( "Test set {0} has already been normalized against {1}."\
						.format( self.source_path, self.normalized_against ) )
			if training_set.featurenames_list != self.featurenames_list:
				raise ValueError("Can't normalize test_set {0} against training_set {1}: Features don't match."\
						.format( self.source_path, training_set.source_path ) )

			for i in range( self.num_features ):
				for class_matrix in self.data_list:
					class_matrix[:,i] -= training_set.feature_minima[i]
					class_matrix[:,i] /= (training_set.feature_maxima[i] -training_set.feature_minima[i])
					class_matrix[:,i] *= 100

			self.normalized_against = training_set.source_path
			

  #=================================================================================
	def FeatureReduce( self, requested_features ):
		"""
		Returns a new TrainingSet that contains a subset of the features
		arg requested_features is a tuple of features
		the returned TrainingSet will have features in the same order as they appear in
		     requested_features
		"""

		# Check that self's faturelist contains all the features in requested_features

		selfs_features = set( self.featurenames_list )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = error_banner + "Feature Reduction error:\n"
			err_str += "The training set '{0}' is missing ".format( self.source_path )
			err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
					len( missing_features_from_req ), len( requested_features ) )
			err_str += "\nDid you forget to convert the feature names into their modern counterparts?"

			raise ValueError( err_str )

		# copy everything but the signature data
		reduced_ts = TrainingSet()
		reduced_ts.source_path = self.source_path + "(feature reduced)"
		reduced_ts.num_classes = self.num_classes
		assert reduced_ts.num_classes == len( self.data_list )
		new_num_features = len( requested_features )
		reduced_ts.num_features = new_num_features
		reduced_ts.num_images = self.num_images
		reduced_ts.imagenames_list = self.imagenames_list[:] # [:] = deepcopy
		reduced_ts.classnames_list = self.classnames_list[:]
		reduced_ts.featurenames_list = requested_features[:]
		reduced_ts.interpolation_coefficients = self.interpolation_coefficients[:]
		reduced_ts.feature_maxima = [None] * new_num_features
		reduced_ts.feature_minima = [None] * new_num_features

		# copy feature minima/maxima
		if self.feature_maxima and self.feature_minima:
			new_index = 0
			for featurename in requested_features:
				old_index = self.featurenames_list.index( featurename )
				reduced_ts.feature_maxima[ new_index ] = self.feature_maxima[ old_index ]
				reduced_ts.feature_minima[ new_index ] = self.feature_minima[ old_index ]
				new_index += 1

		# feature reduce
		for fat_matrix in self.data_list:
			num_imgs_in_class, num_old_features = fat_matrix.shape
			# NB: double parentheses required when calling numpy.zeros(), i guess it's a tuple thing
			new_matrix = np.zeros( ( num_imgs_in_class, new_num_features ) )
			new_column_index = 0
			for featurename in requested_features:
				fat_column_index = self.featurenames_list.index( featurename )
				new_matrix[:,new_column_index] = fat_matrix[:,fat_column_index]
				new_column_index += 1
			reduced_ts.data_list.append( new_matrix )

		return reduced_ts

  #=================================================================================
	def AddSignature( self, signature, class_id_index = None ):
		"""
		@argument signature is a valid signature
		@argument class_id_index identifies the class to which the signature belongs
		"""
		
		if (self.data_list == None) or ( len( self.data_list ) == 0 ) :
			# If no class_id_index is specified, sig goes in first matrix in the list by default
			# make sure there's something there when you try to dereference that index
			self.data_list = []
			self.data_list.append( None )

			self.featurenames_list = signature.names
			self.num_features = len( signature.names )
		else:
			# Make sure features are in order.
			# Feature Reduce will throw an exception if they can't be placed in order
			signature = signature.FeatureReduce( self.featurenames_list )

		# signatures may be coming in out of class order
		while (len( self.data_list ) - 1) < class_id_index:
			self.data_list.append( None )

		if self.data_list[ class_id_index ] == None:
			self.data_list[ class_id_index ] = np.array( signature.values )
		else:
			# vstack takes only one argument, a tuple, thus the extra set of parens
			self.data_list[ class_id_index ] = np.vstack( ( self.data_list[ class_id_index ] ,\
					np.array( signature.values ) ) )


  #=================================================================================
	def CalculateFisherScores( self ):
		"""
		FIXME: implement!
		"""
		raise NotImplementedError()

  #=================================================================================
	def PickleMe( self, pathname = None ):
		"""
		FIXME: pathname needs to end with suffix '.fit.pickled'
		       or TrainingSet.FromPickleFile() won't read it.
		"""

		outfile_pathname = ""
		if pathname != None:
			outfile_pathname = pathname
		else:
			# try to generate a path based on member source_path
			if self.source_path == None or self.source_path == "":
				raise ValueError( "Can't pickle this training set: its 'source_path' member"\
						"is not defined, and you did not specify a file path for the pickle file." )
			if os.path.isdir( self.source_path ):
				# this trainingset was generated from a directory
				# naming convention is /path/to/topleveldir/topleveldir-options.fit.pickled
				root, top_level_dir = os.path.split( self.source_path )
				if self.feature_options != None and self.feature_options != "":
					outfile_pathname = os.path.join( self.source_path, \
							                  top_level_dir + self.feature_options + ".fit.pickled" )
				else:
					outfile_pathname = os.path.join( self.source_path, \
					                      top_level_dir + ".fit.pickled" )
			else:
				# was genearated from a file, could have already been a pickled file
				if self.source_path.endswith( "fit.pickled" ):
					outfile_pathname = self.source_path
				elif self.source_path.endswith( ".fit" ):
					outfile_pathname = self.source_path + ".pickled"
				else:
					outfile_pathname = self.source_path + ".fit.pickled"	

		if os.path.exists( outfile_pathname ):
			print "Overwriting {0}".format( outfile_pathname )
		else:
			print "Writing {0}".format( outfile_pathname )
		with open( outfile_pathname, 'wb') as outfile:
			pickle.dump( self.__dict__, outfile, pickle.HIGHEST_PROTOCOL )


	def DumpNumpyArrays():
		raise NotImplementedError( 'sorry, not implemented yet!' )

# END TrainingSet class definition

#=================================================================================
class ImageClassificationResult:
	path_to_source_image = None

	normalization_factor = None
	marginal_probabilities = []
	predicted_class_name = None
	ground_truth_class_name = None
	interpolated_value = None

	#: For the future:
	kernel_location = None
	#==============================================================
	def PrintToSTDOUT( self, line_item = False ):
		"""
		"""
		if line_item:
			# img name:
			output_str = self.path_to_source_image
			# normalization factor:
			output_str += "\t{val:0.3g}\t".format( val=self.normalization_factor )
			# marginal probabilities:
			output_str += "\t".join(\
					[ "{val:0.3f}".format( val=prob ) for prob in self.marginal_probabilities ] )
			output_str += "\t"
			# actual class:
			if self.ground_truth_class_name:
				output_str += self.ground_truth_class_name + "\t"
			else:
				output_str += "*\t"
			# predicted class:
			output_str += self.predicted_class_name + "\t"
			# interpolated value, if applicable
			if self.interpolated_value:
				output_str += "{val:0.3f}".format( val=self.interpolated_value )
			print output_str
		else:
			print "Image:             \t{0}".format( self.path_to_source_image )
			print "Normalization Factor:\t{0}".format( self.normalization_factor )
			print "Marg. Probabilities:\t" + "\t".join(\
					[ "{val:0.3f}".format( val=prob ) for prob in self.marginal_probabilities ] )
			print "Ground truth class:\t {0}".format( self.ground_truth_class_name ) 
			print "Predicted class:\t {0}".format( self.ground_truth_class_name ) 
			if self.interpolated_value:
				print "Interpolated Value:\t{0}".format( self.interpolated_value )
				#print "Interpolated Value:\t{val=0.3f}".format( val=self.interpolated_value )

class TestSetClassificationResult:
	training_set = None
	test_set = None
	individual_results = []
	classification_accuracy = None

	num_classifications = 0
	num_correct_classifications = 0

	def __init__( self, training_set, test_set ):
		self.training_set = training_set
		self.test_set = test_set

	def GenerateStats( self ):
		for indiv_result in self.individual_results:
			self.num_classifications += 1
			if indiv_result.ground_truth_class_name == indiv_result.predicted_class_name:
				self.num_correct_classifications += 1

		self.classification_accuracy = float( self.num_correct_classifications) / float( self.num_classifications )

	def PrintToSTDOUT( self ):
		if self.classification_accuracy == None:
			self.GenerateStats()

		print "==========================================="
		print "Classification accuracy for this split: {0}".format( self.classification_accuracy )


######################################################################################
# GLOBAL FUNCTIONS
######################################################################################

def _ClassifyOneImageWND5( trainingset, testimg, feature_weights ):
	"""
	Don't call this function directly, use the wrapper functions ClassifyTestSetWND5() or
	ClassifyOneSignatureWND5(), both of which have dummyproofing.

	If you're using this function, your training set data is not continuous
	for N images and M features:
	  trainingset is list of length L of N x M numpy matrices
	  testtile is a 1 x M list of feature values
	NOTE: the trainingset and test image must have the same number of features!!!
	AND: the features must be in the same order!!
	Returns a tuple with norm factor and list of length L of marginal probabilities
	FIXME: what about tiling??
	"""

	#print "classifying..."
	epsilon = np.finfo( np.float ).eps

	num_features_in_testimg = len( testimg ) 
	weights_squared = np.square( feature_weights )

	# initialize
	class_similarities = [0] * trainingset.num_classes

	for class_index in range( trainingset.num_classes ):
		#print "Calculating distances to class {0}".format( class_index )
		num_tiles, num_features = trainingset.data_list[ class_index ].shape
		assert num_features_in_testimg == num_features,\
		"num features {0}, num features in test img {1}".format( num_features, num_test_img_features )

		# create a view
		sig_matrix = trainingset.data_list[ class_index ]
		wnd_sum = 0
		num_collisions = 0

		#print "num tiles: {0}, num_test_img_features {1}".format( num_tiles, num_test_img_features )
		for tile_index in range( num_tiles ):
			#print "{0} ".format( tile_index )
			# epsilon checking for each feature is too expensive
			# do this quick and dirty check until we can figure something else out
			dists = np.absolute( sig_matrix[ tile_index ] - testimg )
			w_dist = np.sum( dists )
			if w_dist < epsilon:
				num_collisions += 1
				continue
			dists = np.multiply( weights_squared, np.square( dists ) )
			w_dist = np.sum( dists )
			# The exponent -5 is the "5" in "WND5"
			class_similarities[ class_index ] += w_dist ** -5
		#print "\n"

		class_similarities[ class_index ] /= ( num_tiles - num_collisions )

	result = ImageClassificationResult()
	norm_factor = sum( class_similarities )
	result.normalization_factor = norm_factor 
	result.marginal_probabilities = [ x / norm_factor for x in class_similarities ]

	return result
#=================================================================================
def ClassifyOneSignatureWND5( training_set, test_sig, feature_weights ):
	"""
	@brief: A wrapper function for _ClassifyOneImageWND5 that does dummyproofing
	@return: classification result
	@returntype: ImageClassificationResult()
	"""

	# check to see if sig is valid
	test_sig.is_valid()

	train_set_len = len( training_set.featurenames_list )
	test_set_len = len( test_sig.names )
	feature_weights_len = len( feature_weights.names )

	if train_set_len != test_set_len or \
	   train_set_len != feature_weights_len or \
	   test_set_len  != feature_weights_len:
		raise ValueError( "Can't classify: one or more of the inputs has a different number of" \
		 "features than the others: training set={0}, test set={1}, feature weights={2}".format( \
				train_set_len, test_set_len, feature_weights_len ) + ". Perform a feature reduce." )

	print "Classifying image '{0}' ({1} features) against test set '{2}' ({3} features)".\
	   format( test_sig.path_to_source_image, train_set_len, training_set.source_path, test_set_len )

	result = _ClassifyOneImageWND5( training_set, test_sig.values, feature_weights.values )

	result.path_to_source_image = test_sig.path_to_source_image
	marg_probs = np.array( result.marginal_probabilities )
	result.predicted_class_name = training_set.classnames_list[ marg_probs.argmax() ]
	# interpolated value, if applicable
	if training_set.interpolation_coefficients is not None:
		interp_val = np.sum( marg_probs * training_set.interpolation_coefficients )
		result.interpolated_value = interp_val

	return result

#=================================================================================
def ClassifyTestSet( training_set, test_set, feature_weights ):
	"""
	@remarks - all three input arguments must have the same number of features,
	and in the same order for this to work properly
	FIXME: What happens when the ground truth is not known? Currently they would all be shoved
	       into class 1, might not be a big deal since class name should be something
	       like "UNKNOWN"
	FIXME: return some python construct that contains classification results
	"""

	train_set_len = len( training_set.featurenames_list )
	test_set_len = len( test_set.featurenames_list )
	feature_weights_len = len( feature_weights.names )

	if train_set_len != test_set_len or \
	   train_set_len != feature_weights_len or \
	   test_set_len  != feature_weights_len:
		raise ValueError( "Can't classify: one or more of the inputs has a different number of" \
				"features than the others: training set={0}, test set={1}, feature weights={2}".format( \
				train_set_len, test_set_len, feature_weights_len ) + ". Perform a feature reduce." )

	print "Classifying test set '{0}' ({1} features) against training set '{2}' ({3} features)".\
	      format( test_set.source_path, test_set_len, training_set.source_path, train_set_len )

	column_header = "image\tnorm. fact.\t"
	column_header +=\
			"".join( [ "p(" + class_name + ")\t" for class_name in training_set.classnames_list ] )

	column_header += "act. class\tpred. class\tpred. val."
	print column_header

	interp_coeffs = None
	if training_set.interpolation_coefficients:
		interp_coeffs = np.array( training_set.interpolation_coefficients )

	split_statistics = TestSetClassificationResult( training_set, test_set )

	for test_class_index in range( test_set.num_classes ):
		num_class_imgs, num_class_features = test_set.data_list[ test_class_index ].shape
		for test_image_index in range( num_class_imgs ):
			one_image_features = test_set.data_list[ test_class_index ][ test_image_index,: ]
			result = _ClassifyOneImageWND5( training_set, one_image_features, feature_weights.values )

			result.path_to_source_image = test_set.imagenames_list[ test_class_index ][ test_image_index ]
			result.ground_truth_class_name = test_set.classnames_list[ test_class_index ]
			marg_probs = np.array( result.marginal_probabilities )
			result.predicted_class_name = training_set.classnames_list[ marg_probs.argmax() ]
			# interpolated value, if applicable
			if interp_coeffs is not None:
				interp_val = np.sum( marg_probs * interp_coeffs )
				result.interpolated_value = interp_val

			result.PrintToSTDOUT( line_item = True )
			split_statistics.individual_results.append( result )

	return split_statistics


#============================================================================
def UnitTest1():
	
	weights_filepath = '/home/colettace/projects/eckley_worms/feature_weights.txt'

	weights = FeatureWeights.NewFromFile( weights_filepath )
	weights.EliminateZeros()
	weights.names = FeatureNameMap.TranslateToNewStyle( weights.names )

	#big_ts = TrainingSet.NewFromFitFile( '/Users/chris/projects/josiah_worms_subset/trunk_train.fit' )
	#big_ts.PickleMe()
	big_ts = TrainingSet.NewFromPickleFile( '/Users/chris/projects/josiah_worms_subset/trunk_train.fit.pickled' )
	big_ts.featurenames_list = FeatureNameMap.TranslateToNewStyle( big_ts.featurenames_list )

	reduced_ts = big_ts.FeatureReduce( weights.names )
	reduced_ts.Normalize()
	
	ClassifyTestSet( reduced_ts, reduced_ts, weights )

#=========================================================================
def UnitTest2():

	ts = TrainingSet.NewFromDirectory( '/home/colettace/projects/eckley_worms/TimeCourse',
	                                   feature_set = "large" )
	ts.PickleMe()
	
#================================================================
def UnitTest3():

	path = "Y24-2-2_GREEN.tif"
	sigs = Signatures.LargeFeatureSet( path )
	sigs.WriteFeaturesToASCIISigFile( "pychrm_calculated.sig" )

#================================================================
def UnitTest4( max_features = 20):
	"""Generate a time curve as a function of number of features used
	to classify a single image"""

	import time
	mommy_feature_weights = FeatureWeights.NewFromPickleFile( "/home/eckleyd/RealTimeClassification/feature_weights_len_2873.weights.pickled" )
	mommy_training_set = TrainingSet.NewFromPickleFile( "/home/eckleyd/RealTimeClassification/FacingL7class_normalized.fit.pickled" )

	filename = "/home/eckleyd/ImageData/SortWorms/ScoredImages/2012-04-13/02/2012-04-13_11:57:50.tif"

	weights_to_be_tested = mommy_feature_weights.Threshold( max_features )

	count = 1
	for name in weights_to_be_tested.names:
		print "{0}\t{1}".format( count, name )
		count += 1

	import sys; sys.exit()
	timings = []

	for number_of_features_to_use in range( 1, max_features ):
		t1 = time.time()
		reduced_feature_weights = mommy_feature_weights.Threshold( number_of_features_to_use )
		sig = Signatures.NewFromFeatureNameList( filename, reduced_feature_weights.names )
		reduced_training_set = mommy_training_set.FeatureReduce( reduced_feature_weights.names )
		normalized_sig = sig.Normalize( reduced_training_set )

		result = ClassifyOneSignatureWND5( reduced_training_set, normalized_sig,\
																			 reduced_feature_weights )
		result.PrintToSTDOUT()
		t2 = time.time()
		timings.append( t2 - t1 )

	count = 1
	for timing in timings:
		print "{0}\t{1}".format( count, timing )
		count += 1

#================================================================
def UnitTest5( max_features = 20 ):
	"""Generate an accuracy curve as a function of number of features used in classification."""

	mommy_feature_weights = FeatureWeights.NewFromPickleFile( "/home/eckleyd/RealTimeClassification/feature_weights_len_2873.weights.pickled" )
	mommy_training_set = TrainingSet.NewFromPickleFile( "/home/eckleyd/RealTimeClassification/FacingL7class_normalized.fit.pickled" )

	split_results = []


	for number_of_features_to_use in range( 1, max_features ):

		reduced_feature_weights = mommy_feature_weights.Threshold( number_of_features_to_use )
		reduced_training_set = mommy_training_set.FeatureReduce( reduced_feature_weights.names )

		split_result = ClassifyTestSet( reduced_training_set, reduced_training_set,\
																			 reduced_feature_weights )
		split_result.PrintToSTDOUT()
		split_results.append( split_result )

	count = 1
	for split_result in split_results:
		print "{0}\t{1}".format( count, split_result.classification_accuracy )
		count += 1

initialize_module()

#================================================================
if __name__=="__main__":
	
	#UnitTest1()
	# UnitTest2()
	# UnitTest3()
	UnitTest4()
	#UnitTest5()
	# pass
