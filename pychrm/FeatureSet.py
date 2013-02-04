""" module for Pychrm intermediates.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2012 National Institutes of Health 

    This library is free software; you can redistribute it and/or              
    modify it under the terms of the GNU Lesser General Public                 
    License as published by the Free Software Foundation; either               
    version 2.1 of the License, or (at your option) any later version.         
                                                                               
    This library is distributed in the hope that it will be useful,            
    but WITHOUT ANY WARRANTY; without even the implied warranty of             
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          
    Lesser General Public License for more details.                            
                                                                               
    You should have received a copy of the GNU Lesser General Public           
    License along with this library; if not, write to the Free Software        
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Written by:  Christopher Coletta <christopher.coletta [at] nih [dot] gov>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Pychrm TODO (in no particular order):
* Implement multiprocessing
* Implement kernel classifier
* Implement database hookup
* OMERO integration
* Use scipy.optimize.curve_fit (http://www.scipy.org/Cookbook/FittingData) to
  further refine continuous classifiers.
* Implement dendrogram using scipy.cluster or Biopython.Phylo, including the ability
  to create dendrograms for individual images as well as entire image classes.
* Classify images against multiple classifiers, e.g., filter for confluence, then classify.
* Break out classes into their own modules.
* Generate documantation, getting started, tutorials, examples.
  Use some model-building module to map out object model (Pylint/Pyreverse)
* Migrate over to using the unittest module, and check in unit test related files 
* Implement Channels.
* Replace/optimize current slow implementation of Chebyshev coefficients
  which is the biggest pig of the algorithms
* Feature family analysis. how to get highest accuracy using least amount of feature groups:
* Leave one out analysis
* Wndchrm backend c++ cleanup - potential use of ctypes to facilitate shared memory,
  construction of new sharable ImageMatrix's.
* Import feature calculation modules only if required
* Implement ability to define continuous classification bin walls, based on which
  classification (i.e., predicted value) can be called "correct" or "incorrect"
  and consequently confusion/similarity matrices can be built.
* Add ability to remove sample/image/time by specifying name or index
"""

# ===========================================================
# The following are some Pychrm-specific module imports:

# pychrm.py has the definitions of all the SWIG-wrapped primitive
# C++ WND_CHARM objects.
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

# ============================================================
# Imports from Python core or other installed packages
import numpy as np
try:
	import cPickle as pickle
except:
	import pickle

import os
import os.path 
import itertools
import copy



# ============================================================
# BEGIN: Initialize module level globals
Algorithms = None
Transforms = None
small_featureset_featuregroup_strings = None
large_featureset_featuregroup_strings = None
small_featureset_featuregroup_list = None
large_featureset_featuregroup_list = None

error_banner = "\n*************************************************************************\n"


def initialize_module(): 
	"""If you're going to calculate any signatures, you need this stuff.
	FIXME: Rig this stuff to load only on demand."""
	
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
		fg = FeatureGroup.NewFromString( fg_str )
		if fg:
			small_featureset_featuregroup_list.append( fg )

	large_featureset_featuregroup_list = []
	for fg_str in large_featureset_featuregroup_strings.splitlines():
		fg = FeatureGroup.NewFromString( fg_str )
		if fg:
			large_featureset_featuregroup_list.append( fg )

	# while we're debugging, raise exceptions for numerical weirdness, since it all has to be dealt with somehow
	# In cases where numerical weirdness is expected and dealt with explicitly, these exceptions are
	# temporarily turned off and then restored to their previous settings.
	np.seterr (all='raise')

def output_railroad_switch( method_that_prints_output ):
	"""This is a decorator that optionally lets the user specify a file to which to redirect
	STDOUT. To use, you must use the keyword argument "output_filepath" and optionally
	the keyword argument "mode" """

	def print_method_wrapper( *args, **kwargs ):
		
		retval = None
		if "output_filepath" in kwargs:
			output_filepath = kwargs[ "output_filepath" ]
			del kwargs[ "output_filepath" ]
			if "mode" in kwargs:
				mode = kwargs[ "mode" ]
				del kwargs[ "mode" ]
			else:
				mode = 'w'
			print 'Saving output of function "{0}()" to file "{1}", mode "{2}"'.format(\
			      method_that_prints_output.__name__, output_filepath, mode )
			import sys
			backup = sys.stdout
			sys.stdout = open( output_filepath, mode )
			retval = method_that_prints_output( *args, **kwargs )
			sys.stdout.close()
			sys.stdout = backup
		else:
			retval = method_that_prints_output( *args, **kwargs )
		return retval

	return print_method_wrapper

def normalize_by_columns ( full_stack, mins = None, maxs = None ):
	"""This is a global function to normalize a matrix by columns.
	If numpy 1D arrays of mins and maxs are provided, the matrix will be normalized against these ranges
	Otherwise, the mins and maxs will be determined from the matrix, and the matrix will be normalized
	against itself. The mins and maxs will be returned as a tuple.
	Out of range matrix values will be clipped to min and max (including +/- INF)
	zero-range columns will be set to 0.
	NANs in the columns will be set to 0.
	The normalized output range is hard-coded to 0-100
	"""
# Edge cases to deal with:
#   Range determination:
#     1. features that are nan, inf, -inf
#        max and min determination must ignore invalid numbers
#        nan -> 0, inf -> max, -inf -> min
#   Normalization:
#     2. feature values outside of range
#        values clipped to range (-inf to min -> min, max to inf -> max) - leaves nan as nan
#     3. feature ranges that are 0 result in nan feature values
#     4. all nan feature values set to 0

# Turn off numpy warnings, since e're taking care of invalid values explicitly
	oldsettings = np.seterr(all='ignore')
	if (mins is None or maxs is None):
		# mask out NANs and +/-INFs to compute min/max
		full_stack_m = np.ma.masked_invalid (full_stack, copy=False)
		maxs = full_stack_m.max (axis=0)
		mins = full_stack_m.min (axis=0)

	# clip the values to the min-max range (NANs are left, but +/- INFs are taken care of)
	full_stack.clip (mins, maxs, full_stack)
	# remake a mask to account for NANs and divide-by-zero from max == min
	full_stack_m = np.ma.masked_invalid (full_stack, copy=False)

	# Normalize
	full_stack_m -= mins
	full_stack_m /= (maxs - mins)
	# Left over NANs and divide-by-zero from max == min become 0
	full_stack = full_stack_m.filled (0) * 100

	# return settings to original
	np.seterr(**oldsettings)

	return (mins,maxs)


def CheckIfClassNamesAreInterpolatable( classnames_list ):

	import re
	interp_coeffs = []
	for class_name in classnames_list:
		m = re.search( r'(\d*\.?\d+)', class_name )
		if m:
			interp_coeffs.append( float( m.group(1) ) )
		else:
			interp_coeffs = None
			break
	return interp_coeffs
# END: Initialize module level globals
#===============================================================


# BEGIN: Class definitions for WND-CHARM intermediate objects

#############################################################################
# class definition of SampleImageTiles
#############################################################################

class SampleImageTiles (object):
	"""SampleImageTiles is an image iterator wrapper (the iterator is the sample method).
	The iterator is wrapped to provide additional information such as the number of samples that will
	be extracted from the image, as well as information about each sample after calling the sample method.
	Each call to sample returns the next pychrm.SharedImageMatrix in the sample set.
	The constructor has three required parameters.
	The image parameter can be a path to an image file or a pychrm.SharedImageMatrix
	The x and y parameters can specify the number of non-overlapping samples in each dimension (is_fixed parameter is False),
	or the dimentions of each sample (is_fixed parameter is True).
	Example usage:
		image_iter = SampleImageTiles (input_image, size_x, size_y, True)
		print "Number of samples = "+str (image_iter.samples)
		for sample in image_iter.sample():
			print "({0},{1}) : ({2},{3})".format (
				image_iter.current_x, image_iter.current_y, sample.width, sample.height)
	"""

	def __init__( self, image_in, x, y, is_fixed = False):
		if isinstance (image_in, str):
			if not os.path.exists( image_in ):
				raise ValueError( "The file '{0}' doesn't exist, maybe you need to specify the full path?".format( image_in ) )
			self.image = pychrm.SharedImageMatrix()
			if 1 != self.image.OpenImage( image_in, 0, None, 0, 0 ):
				raise ValueError( 'Could not build an SharedImageMatrix from {0}, check the file.'.format( image_in ) )
		elif isinstance (image_in, pychrm.SharedImageMatrix):
			self.image = image_in
		else:
			raise ValueError("image parameter 'image_in' is not a string or a pychrm.SharedImageMatrix")

		if (is_fixed):
			self.tile_width = x
			self.tile_height = y
			self.tiles_x = int (self.image.width / x)
			self.tiles_y = int (self.image.height / y)
		else:
			self.tile_width = int (self.image.width / x)
			self.tile_height = int (self.image.height / y)
			self.tiles_x = x
			self.tiles_y = y

		self.samples = self.tiles_x * self.tiles_y

	def sample(self):
		width = self.tile_width
		height = self.tile_height
		max_x = self.image.width
		max_y = self.image.height
		original = self.image
		current_y = 0
		self.current_y = current_y
		while current_y + height <= max_y:
			current_x = 0
			self.current_x = current_x
			while current_x + width <= max_x:
				yield pychrm.SharedImageMatrix (original, current_x, current_y, current_x+width-1, current_y+height-1,0,0)
				current_x = current_x + width
				self.current_x = current_x
			current_y = current_y + height
			self.current_y = current_y


#############################################################################
# class definition of FeatureVector
#############################################################################
# FeatureVector inherits from class "object" and is thus a Python "new-style" class
class FeatureVector(object):
	"""A FeatureVector is simply a list of doubles and a list of corresponding strings (i.e., names).
	The lengths of the two lists must always be equal. A FeatureVector may be rendered 
	corrupted if the order of one list is changed but not the other. In the future,
	we might consider Python tuples instead of Python lists to store data since 
	tuples are immutable."""

	names = None
	values = None
	#================================================================
	def __init__( self, data_dict = None):
		"""@brief: constructor"""

		self.names = []
		self.values = []

		if data_dict:
			if "values" in data_dict:
				self.values = data_dict[ "values" ]
			if "names" in data_dict:
				self.names = data_dict[ "names" ]
	#================================================================
	def is_valid( self ):
		"""@brief: an instance should know all the criteria for being a valid FeatureVector"""
		if len( self.values ) != len( self.names ):
			raise RuntimeError( "Instance of {0} is invalid: ".format( self.__class__ ) + \
			  "different number of values ({0}) and names ({1}).".format( \
			  len( self.values ), len( self.names ) ) )
		return True

#############################################################################
# class definition of FeatureWeights
#############################################################################
class FeatureWeights( FeatureVector ):
	"""FeatureWeights is an abstract base class inheriting from "FeatureVector"
	that comprises one-half of a WND-CHARM classifier (the other being a FeatureSet).
	
	It is a list of strings which are the names of
	individual image descriptors (features) and a corresponding list of doubles which
	are the weights assigned to those features. Since WND-CHARM is a generalized
	pattern recognition algorithm that calculates the same image descriptors for all
	images, it is through features weights that trained classifiers can zero-in on only
	those features which provide distinctiveness across classes and ignore noisy features.
	Thus any instance of a FeatureWeights class is context-specific."""

	name = None
	associated_training_set = None

	def __init__( self, data_dict = None, name = None ):
		# call parent constructor
		super( FeatureWeights, self ).__init__( data_dict )
		self.name = name

	#================================================================
	@classmethod
	def NewFromFile( cls, weights_filepath ):
		"""Load feature weights from a C++ WND-CHARM-sytle feature weights text file,
		i.e., the one created by using the command line "wndchrm train -vw/path/to/weights.txt"""
		raise NotImplementedError

	#================================================================
	def Threshold( self, num_features_to_be_used  ):
		"""@breif Returns an instance of a FeatureWeights class with the top n relevant features in that order"""
		raise NotImplementedError

	#================================================================
	@classmethod
	def NewFromFeatureSet( cls, num_features_to_be_used  ):
		"""@breif Calculate FeatureWeights from a FeatureSet"""
		raise NotImplementedError

	#================================================================
	def PickleMe( self, pathname = None ):
		"""Creates a Python pickle of this class and saves it to a file specified by pathname"""

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

	#================================================================
	@output_railroad_switch
	def Print( self ):
		"""@breif Prints out feature values and statistics"""
		raise NotImplementedError



#############################################################################
# class definition of FisherFeatureWeights
#############################################################################
class FisherFeatureWeights( FeatureWeights ):
	"""
	FisherFeatureWeights is a concrete class using feature weighting that uses the
	Fisher Discriminant as a weight for image descriptors. Features that have a high
	Fisher score have high inter-class variability and low intra-class variability.

	This type of feature weight is one result of training an image classifier from a 
	training set that has distinct discrete classes, e.g., is a cell mono- or binucleate.
	Fisher classifiers can be used on a continuous morphology regime if the user creates
	groups or bins for the training images.
	"""

	def __init__( self, data_dict = None, name = None):
		"""Simply calls parent constructor"""
		super( FisherFeatureWeights, self ).__init__( data_dict, name )


	#================================================================
	@classmethod
	def NewFromFile( cls, weights_filepath ):
		"""Written to read in files created by wndchrm -vw/path/to/weightsfile.txt"""

		weights = cls()
		with open( weights_filepath, 'r' ) as weights_file:
			for line in weights_file:
				# split line "number <space> name"
				feature_line = line.strip().split( " ", 1 )
				weights.values.append( float( feature_line[0] ) )
				weights.names.append( feature_line[1] )

		return weights

	#================================================================
	@classmethod
	def NewFromFeatureSet( cls, training_set ):
		"""Takes a FeatureSet_Discrete as input and calculates a Fisher score for 
		each feature. Returns a newly instantiated instance of FisherFeatureWeights.

		For:
		N = number of classes
		F = number of features
		It = total number of images in training set
		Ic = number of images in a given class
		"""

		if training_set == None:
			import inspect
			form_str = 'You passed in a None as a training set to the function {0}.{1}'	
			raise ValueError( form_str.format( cls.__name__, inspect.stack()[1][3] ) )
		if not training_set.__class__.__name__ == "FeatureSet_Discrete":
			raise ValueError( "Cannot create Fisher weights from anything other than a FeatureSet_Discrete." )

		# 1D matrix 1 * F
		population_means = np.mean( training_set.data_matrix, axis = 0 )

		# WARNING, this only works in python27:
		# ====================================
		# If 'training_set' is a balanced training set (i.e., same number of images
		# in each class), you can use pure matrix calls without iteration:

		# 3D matrix N * Ic * F
		#all_images_classes_separate = np.array( training_set.data_list )

		#if len( all_images_classes_separate.shape ) == 3:

			# 2D matrix N * F
		#	intra_class_means = np.mean( all_images_classes_separate, axis = 1 )
			# 2D matrix N * F
		#	intra_class_variances = np.var( all_images_classes_separate, axis = 1 )

		#else:
		# ====================================

		# 2D matrix shape N * F
		intra_class_means = np.empty( \
				 ( training_set.num_classes, len( training_set.featurenames_list ) ) )
		# 2D matrix shape N * F
		intra_class_variances = np.empty( \
				 ( training_set.num_classes, len( training_set.featurenames_list ) ) )
		
		class_index = 0
		for class_feature_matrix in training_set.data_list:
			intra_class_means[ class_index ] = np.mean( class_feature_matrix, axis=0 )
			intra_class_variances[ class_index ] = np.var( class_feature_matrix, axis=0 )
			class_index += 1


		# 1D matrix 1 * F
		# we deal with NANs/INFs separately, so turn off numpy warnings about invalid floats.
		# for the record, in numpy:
		#     1./0. = inf, 0./inf = 0., 1./inf = 0. inf/0. = inf, inf/inf = nan
		#     0./0. = nan,  nan/0. = nan, 0/nan = nan, nan/nan = nan, nan/inf = nan, inf/nan = nan
		# We can't deal with NANs only, must also deal with pos/neg infs
		# The masked array allows for dealing with "invalid" floats, which includes nan and +/-inf
		oldsettings = np.seterr(invalid='ignore')
		feature_weights_m =  np.ma.masked_invalid (
			np.square( population_means - intra_class_means ).sum( axis = 0 ) /
		    (training_set.num_classes - 1) / np.mean( intra_class_variances, axis = 0 )
		    )
		# return numpy error settings to original
		np.seterr(**oldsettings)


		new_fw = cls()
		new_fw.names = training_set.featurenames_list[:]
		# the filled(0) method of the masked array sets all nan and infs to 0
		new_fw.values = feature_weights_m.filled(0).tolist()
		new_fw.associated_training_set = training_set
	
		return new_fw

	#================================================================
	def EliminateZeros( self ):
		"""Eliminates any features with a weight of zero, and returns a new instance of
		FisherFeatureWeights without those features."""

		new_weights = FisherFeatureWeights()
		scores = zip( self.names, self.values )
		nonzero_scores = [ (name, weight) for name, weight in scores if weight != 0 ]
		new_weights.names, new_weights.values = zip( *nonzero_scores )
		return new_weights

	#================================================================
	def Threshold( self, num_features_to_be_used = None ):
		"""Returns an instance of a FisherFeatureWeights class with the top n relevant features 
		in order."""

		# Default is top 15% of features
		if num_features_to_be_used is None:
			num_features_to_be_used = int( len( self.values ) * 0.15 )
		elif num_features_to_be_used > len( self.values ):
			raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
			                      format( len( self.values ), num_features_to_be_used ) )

		new_weights = self.__class__()
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

		new_weights.associated_training_set = self.associated_training_set

		return new_weights

	#================================================================
	def Slice( self, start_index, stop_index):
		"""Return a new instance of FisherFeatureWeights containing a chunk
		of middle-ranked features."""
		
		min_index = None
		max_index = None

		if stop_index > start_index:
			min_index = start_index
			max_index = stop_index
		else:
			min_index = stop_index
			max_index = start_index

		if (min_index < 0) or ( max_index > len( self.values ) ):
			raise ValueError( 'Cannot slice, check your start and stop indices.' )

		new_weights = self.__class__()
		raw_featureweights = zip( self.names, self.values )

		use_these_feature_weights = \
				list( itertools.islice( raw_featureweights, min_index, max_index ) )
		
		# we want lists, not tuples!
		new_weights.names, new_weights.values =\
		  [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

		new_weights.associated_training_set = self.associated_training_set

		return new_weights

	#================================================================
	@output_railroad_switch
	def Print( self ):
		"""Prints out feature weight values and statistics."""
		print "Fisher Feature Weight set {0}:".format( self.name )
		print "Rank\tValue\tName"
		print "====\t=====\t===="
		for i in range( 0, len( self.names ) ):
			#print "{0}\t{1:.5f}\t{2}".format( i+1, self.values[i], self.names[i] )
			print "{0:.6f}\t{1}".format( self.values[i], self.names[i] )
		print ""


#############################################################################
# class definition of ContinuousFeatureWeights
#############################################################################
class ContinuousFeatureWeights( FeatureWeights ):
	"""A concrete class that calculates correlation coefficients as well as 
	regression parameters for each feature. Features are weighted based on how well
	they linearly correlate (i.e., high Pearson correlation coefficient) with an experimental 
	variable.

	An example system where a continuous classifier could be used could be
	would be defining a spectrum of morphology across age or dose response."""

	slopes = None
	intercepts = None
	pearson_coeffs = None
	pearson_stderrs = None
	pearson_p_values = None
	spearman_coeffs = None
	spearman_p_values = None
	
	def __init__( self, data_dict = None ):
		"""Constructor"""
		super( ContinuousFeatureWeights, self ).__init__( data_dict )
		self.slopes = []
		self.intercepts = []
		self.pearson_coeffs = []
		self.pearson_stderrs = []
		self.pearson_p_values = []
		self.spearman_coeffs = []
		self.spearman_p_values = []


	#================================================================
	@classmethod
	def NewFromFeatureSet( cls, training_set ):
		"""Calculate regression parameters and correlation statistics that fully define
		a continuous classifier.

		At present the feature weights are proportional the Pearson correlation coefficient
		for each given feature."""
		
		from scipy import stats

		# Known issue: running stats.linregress() with np.seterr (all='raise') has caused
		# arithmetic underflow (FloatingPointError: 'underflow encountered in stdtr' )
		# I think this is something we can safely ignore in this function, and return settings
		# back to normal at the end. -CEC
		np.seterr (under='ignore')

		matrix = training_set.data_matrix
		#FIXME: maybe add some dummyproofing to constrain incoming array size

		new_fw = cls()

		new_fw.associated_training_set = training_set
		if training_set.source_path:
			new_fw.name = cls.__name__ + ' from training set "' + training_set.source_path + '"'

		r_val_squared_sum = 0
		#r_val_cubed_sum = 0

		ground_truths = np.array( [float(val) for val in training_set.ground_truths] )

		for feature_index in range( training_set.num_features ):
			feature_values = matrix[:,feature_index]

			slope, intercept, pearson_coeff, p_value, std_err = \
			             stats.linregress( ground_truths, feature_values )

			new_fw.names.append( training_set.featurenames_list[ feature_index ] )
			new_fw.pearson_coeffs.append( pearson_coeff )
			new_fw.slopes.append( slope )
			new_fw.intercepts.append( intercept )
			new_fw.pearson_stderrs.append( std_err )
			new_fw.pearson_p_values.append( p_value )

			r_val_squared_sum += pearson_coeff * pearson_coeff
			#r_val_cubed_sum += pearson_coeff * pearson_coeff * pearson_coeff

			spearman_coeff, spearman_p_val = stats.spearmanr( ground_truths, feature_values )
			new_fw.spearman_coeffs.append( spearman_coeff )
			new_fw.spearman_p_values.append( spearman_p_val )

		new_fw.values = [val*val / r_val_squared_sum for val in new_fw.pearson_coeffs ]
		#new_fw.values = [val*val*val / r_val_cubed_sum for val in new_fw.pearson_coeffs ]
		
		# Reset numpy 
		np.seterr (all='raise')

		return new_fw

	#================================================================
	def Threshold( self, num_features_to_be_used = None  ):
		"""Returns a new instance of a ContinuousFeatureWeights class derived from this
		instance where the number of features has been reduced to only the top n features,
		where n is specified by the num_features_to_be_used argument."""

		# Default is top 15% of features
		if num_features_to_be_used is None:
			num_features_to_be_used = int( len( self.values ) * 0.15 )
		elif num_features_to_be_used > len( self.values ):
			raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
			                      format( len( self.values ), num_features_to_be_used ) ) 

		new_weights = self.__class__()
		if self.name:
			if num_features_to_be_used == len( self.names ):
				new_weights.name = self.name + " (rank-ordered)"
			else:
				new_weights.name = self.name + " (top {0} features)".format( num_features_to_be_used )


		abs_val_pearson_coeffs = [ abs( val ) for val in self.pearson_coeffs ]
		raw_featureweights = zip( self.names, abs_val_pearson_coeffs, self.pearson_coeffs, \
		    self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
		    self.spearman_coeffs, self.spearman_p_values )
		
		# sort from max to min
		# sort by the second item in the tuple, i.e., index 1
		sort_func = lambda feat_a, feat_b: cmp( feat_a[1], feat_b[1] ) 

		sorted_featureweights = sorted( raw_featureweights, sort_func, reverse = True )
		
		# take most correllated features, both positive and negative
		use_these_feature_weights = list( itertools.islice( \
			sorted_featureweights, num_features_to_be_used ) )

		# we want lists, not tuples!
		new_weights.names, abs_pearson_coeffs, new_weights.pearson_coeffs, new_weights.slopes, \
		    new_weights.intercepts, new_weights.pearson_stderrs, new_weights.pearson_p_values,\
		    new_weights.spearman_coeffs, new_weights. spearman_p_values =\
		      [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

		r_val_sum = 0
		for val in abs_pearson_coeffs:
			r_val_sum += val
		new_weights.values = [ val / r_val_sum for val in abs_pearson_coeffs ]

		new_weights.associated_training_set = self.associated_training_set

		return new_weights

	#================================================================
	def Slice( self, start_index, stop_index):
		"""Return a new instance of ContinuousFeatureWeights populated with a
		chunk of middle-ranked features specified by arguments start_index and stop_index."""
		
		min_index = None
		max_index = None

		if stop_index > start_index:
			min_index = start_index
			max_index = stop_index
		else:
			min_index = stop_index
			max_index = start_index

		if (min_index < 0) or ( max_index > len( self.values ) ):
			raise ValueError( 'Cannot slice, check your start and stop indices.' )

		new_weights = self.__class__()
		if self.name:
			new_weights.name = self.name + " (sliced {0}-{1})".format( min_index, max_index )

		abs_val_pearson_coeffs = [ abs( val ) for val in self.pearson_coeffs ]
		raw_featureweights = zip( self.names, abs_val_pearson_coeffs, self.pearson_coeffs, \
		    self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
		    self.spearman_coeffs, self.spearman_p_values )

		use_these_feature_weights = \
				list( itertools.islice( raw_featureweights, min_index, max_index ) )
		
		new_weights.names, abs_pearson_coeffs, new_weights.pearson_coeffs, new_weights.slopes, \
		    new_weights.intercepts, new_weights.pearson_stderrs, new_weights.pearson_p_values,\
		    new_weights.spearman_coeffs, new_weights. spearman_p_values =\
		      [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

		r_val_sum = 0
		for val in abs_pearson_coeffs:
			r_val_sum += val
		new_weights.values = [ val / r_val_sum for val in abs_pearson_coeffs ]

		new_weights.associated_training_set = self.associated_training_set

		return new_weights

	#================================================================
	@output_railroad_switch
	def Print( self, print_legend = True ):
		"""@breif Prints out feature values and statistics"""

		header_str = "Continuous feature weight set"
		if self.name:
			header_str += ' "{0}"'.format( self.name )
		header_str += ":"
		print header_str
		print "Total num features: {0}".format( len( self.values ) )
		print "-----------------------------------"

		if print_legend:
			print "Legend:"
			print "IFW - Feature weight applied to the individual feature"
			print "IPC - Pearson correlation coefficient of feature values vs ground truth"
			print "IPE - Standard Error of IPC"
			print "IPP - P-value of IPC"
			print "ISC - Spearman correlation coefficient of feature values vs ground truth"
			print "IPP - P-value of ISC"
			print ""
		print "NUM\tIFW\tIPC\tIPE\tIPP\tISC\tIPP\tNAME"
		print "===\t===\t===\t===\t===\t===\t===\t===="	
		for i in range( len( self.values ) ):
			line_item = "{0}\t".format( i + 1 )
			line_item += "{0:2.4f}\t".format( self.values[i] )
			line_item += "{0:2.4f}\t".format( self.pearson_coeffs[i] )
			line_item += "{0:2.4f}\t".format( self.pearson_stderrs[i] )
			line_item += "{0:2.4f}\t".format( self.pearson_p_values[i] )
			line_item += "{0:2.4f}\t".format( self.spearman_coeffs[i] )
			line_item += "{0:2.4f}\t".format( self.spearman_p_values[i] )
			line_item += self.names[i]
			print line_item


#############################################################################
# class definition of Signatures
#############################################################################
class Signatures( FeatureVector ):
	"""
	Signatures is a concrete class inheriting from FeatureVector
	that contains the image descriptor values, a.k.a. "Signatures", 
	for a single image or ROI. It is the values contained in 
	Signatures that are grouped together to form FeatureSets.

	This class is mostly agnostic w.r.t. the types of image descriptors contained herein.
	Previous implementations of WND-CHARM have defined a "small feature set" of
	1025 specific image descriptors as well as a "large feature set" of 2885 features.
	These specific lists have been preserved here for convenience and interchangeability.
	"""

	# in the future, signatures can be loaded from be from a database as well.
	source_file = None
	path_to_image = None
	path_to_sigfile = None

	#: Captures some of the image descriptor meta-data, such as whether we
	#: used the large feature set to generate these features ("-l")

	options = ""

	#================================================================
	def __init__( self ):
		"""@brief: constructor"""

		# call parent constructor
		super( Signatures, self ).__init__()

	#================================================================
	@classmethod
	def NewFromTiffFile( cls, imagepath, options = "-l" ):
		"""@brief Loads precalculated sigs, or calculates new ones

		This function is meant to be the quarterback in deciding what sig
		loading function is called.
		@argument path - path to a tiff file
		@return An instance of the class Signatures for image with sigs calculated."""

		if options == "-l":
			return cls.LargeFeatureSet( imagepath, options.replace( "-l", "" ) )
		else:
			return cls.SmallFeatureSet( imagepath, options )
		# FIXME: else if a feature list or feature group list were specified...

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

		the_sigs = cls()

		# Try to find a file that corresponds to this image with signatures in it.
		# derive the name of the sig/pysig file that would exist if there were one
		# given the options requested

		sigpath = None
		root, extension = os.path.splitext( imagepath )
		options_str = options if options else ""
		if not "-l" in options_str:
			options_str += "-l"

		read_in_sig_file = False
		if os.path.exists( root + options_str + ".pysig" ):
			sigpath = root + options_str + ".pysig"
			read_in_sig_file = True
		elif os.path.exists( root + options_str + ".sig" ):
			sigpath = root + options_str + ".sig" 
			read_in_sig_file = True
		else:
			sigpath = root + options_str + ".pysig"

		if read_in_sig_file:
			print "{0}: reading in signature from file: {1}".format( imagepath, sigpath )
			the_sigs = cls.NewFromSigFile( sigpath, imagepath, options_str )


		# Check to see what, if anything was loaded from file. The file could be corrupted
		# incomplete, or calculated with different options, e.g., -S1441
		# FIXME: Here's where you'd calculate a small subset of features
		# and see if they match what was loaded from file
		if len( the_sigs.names ) <= 0:
			# All hope is lost. Calculate sigs
			print "Calculating large feature set for file: {0}".format( imagepath )
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
	def NewFromFeatureGroupList( cls, image_path_or_mat, feature_groups, options = None ):
		"""@brief calculates signatures
		@argument image_path_or_mat - path to a tiff file as a string or a pychrm.SharedImageMatrix object
		"""

		if isinstance (image_path_or_mat, str):
			path_to_image = image_path_or_mat
			if not os.path.exists( path_to_image ):
				raise ValueError( "The file '{0}' doesn't exist, maybe you need to specify the full path?".format( outfile_pathname ) )
			original = pychrm.SharedImageMatrix()
			if 1 != original.OpenImage( path_to_image, 0, None, 0, 0 ):
				raise ValueError( 'Could not build an SharedImageMatrix from {0}, check the file.'.\
					format( path_to_image ) )
		elif isinstance (image_path_or_mat, pychrm.SharedImageMatrix):
			original = image_path_or_mat
			path_to_image = "sample" # should really get this from SharedImageMatrix
		else:
			raise ValueError("image parameter 'image_path_or_mat' is not a string or a pychrm.SharedImageMatrix")

		print path_to_image
		im_cache = {}
		im_cache[ '' ] = original

		# instantiate an empty Signatures object
		signatures = cls()
		signatures.source_file = path_to_image
		signatures.path_to_image = path_to_image
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
		"""@brief calculates signatures"""

		work_order, num_features = GenerateWorkOrderFromListOfFeatureStrings( feature_names )
		sig = cls.NewFromFeatureGroupList( path_to_image, work_order, options )
		return sig.FeatureReduce( feature_names ) 

	#================================================================
	@classmethod
	def FromPickle( cls, path ):
		"""Not implemented."""
		raise NotImplementedError()

	#================================================================
	def PickleMe( self ):
		"""Not implemented."""
		raise NotImplementedError()

	#================================================================
	@classmethod
	def NewFromSigFile( cls, sigfile_path, image_path = None, options = None ):
		"""@argument sigfile_path must be a .sig or a .pysig file
		
		@return  - An instantiated signature class with feature names translated from
		           the old naming convention, if applicable.
		@remarks - old style sig files don't know about their options other than 
		           the information contained in the name. In the future pysig files
		           may keep that info within the file. Thus, for now, the argument
		           options is something set externally rather than gleaned from
		           reading the file."""

		if sigfile_path == None or sigfile_path == "":
			raise ValueError( "The path to the signature file to be opened was empty or null. "\
			                  "Check arguments to the function {0}.NewFromSigFile and try again.".\
			                  format( cls.__name__ ) )
		if not sigfile_path.endswith( (".sig", ".pysig" ) ):
			raise ValueError( "The file {0} isn't a signature file (doesn't end in .sig or .pysig)".\
			                  format( sigfile_path ) )

		#print "Loading features from sigfile {0}".format( sigfile_path )

		signatures = cls()
		# FIXME: Do we care about the .tif?

		signatures.source_file = sigfile_path
		signatures.path_to_sigfile = sigfile_path
		signatures.path_to_image = image_path
		signatures.options = options
 
		with open( sigfile_path ) as infile:
			linenum = 0
			for line in infile:
				if linenum == 0:
					# The class id here may be trash
					signatures.class_id = line.strip()
				elif linenum == 1:
					# We've already assigned self.source_file
					# the path in the sig file may be trash anyway
					#signatures.source_file = line.strip()
					pass
				else:
					value, name = line.strip().split( ' ', 1 )
					signatures.values.append( float( value ) )	
					signatures.names.append( name )
				linenum += 1
			#print "Loaded {0} features.".format( len( signatures.values ) )

		# Check if the feature name follows the old convention
		signatures.names = FeatureNameMap.TranslateToNewStyle( signatures.names ) 
		
		return signatures
	
	#================================================================
	def WriteFeaturesToASCIISigFile( self, filepath = None ):
		"""Write features to a .pysig file, in the same format as a .sig file created 
		
		If filepath is specified, you get to name it whatever you want and put it
		wherever you want. Otherwise, it's named according to convention and placed 
		next to the image file in its directory."""

		self.is_valid()

		outfile_path = ""
		if not filepath or filepath == "":
			if not self.source_file or self.source_file == "":
				raise ValueError( "Can't write sig file. No filepath specified in function call, and no path associated with this instance of Signatures." )
			outfile_path = self.source_file

			path, filename = os.path.split( outfile_path )
			
			# if the filepath is None, you don't have to check if it exists,
			# it's gonna be created right here
			if not filepath == None:
				if not os.path.exists( path ):
					raise ValueError( 'Invalid path {0}'.format( path ) )

			filename_parts = filename.rsplit( '.', 1 )
			if self.options and self.options is not "":
				outfile_path = "{0}{1}.pysig".format( filename_parts[0],\
																					self.options if self.options else "" )
			else:
				outfile_path = "{0}.pysig".format( filename_parts[0] )
			outfile_path = os.path.join( path, outfile_path )
		else:
			outfile_path = filepath

		if os.path.exists( outfile_path ):
			print "Overwriting {0}".format( outfile_path )
		else:
			print 'Writing signature file "{0}"'.format( outfile_path )

		self.path_to_sigfile = outfile_path

		with open( outfile_path, "w" ) as out_file:
			# FIXME: line 2 contains class membership, just hardcode a number for now
			out_file.write( "0\n" )
			out_file.write( "{0}\n".format( self.source_file ) )
			for i in range( 0, len( self.names ) ):
				out_file.write( "{val:0.6f} {name}\n".format( val=self.values[i], name=self.names[i] ) )

	#================================================================
	def FeatureReduce( self, requested_features ):
		"""Returns a new instance of Signatures that contains only those features specified
		by name via the argument requested_features."""

		if self.names == requested_features:
			return self

		selfs_features = set( self.names )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = error_banner + "Feature Reduction error:\n"
			err_str += "The signatures set for image file '{0}' ".format( self.source_file )
			err_str += "is missing {0}".format( len( missing_features_from_req ) )
			err_str += "/{1} features that were requested in the feature reduction list.".format(\
			             len( requested_features ) )
			raise ValueError( err_str )

		print "Performing feature reduce/reorder on signatures..."
		# The featurenames sets are either out of order or contain extraneous features.
		dictionary = dict( zip( self.names, self.values ) )
		reduced_sigs = Signatures()
		reduced_sigs.source_file = self.source_file
		reduced_sigs.path_to_image = self.path_to_image
		reduced_sigs.path_to_sigfile = self.path_to_sigfile
		reduced_sigs.options = self.options
			
		for new_name in requested_features:
			reduced_sigs.names.append( new_name )
			reduced_sigs.values.append( dictionary[ new_name ] )
		
		return reduced_sigs

	#================================================================
	def Normalize( self, training_set ):
		"""
		FIXME: should there be some flag that gets set if this sig has 
		already been normalized??
		
		@return: None
		"""

		if training_set.featurenames_list != self.names:
			raise ValueError("Can't normalize signature for {0} against training_set {1}: Features don't match."\
		  .format( self.source_file, training_set.source_path ) )
		
		print "Normalizing features"
		my_features = np.array( self.values )
		normalize_by_columns (my_features, training_set.feature_minima, training_set.feature_maxima)
		self.values = my_features.tolist()

# end definition class Signatures

#############################################################################
# class definition of FeatureGroup
#############################################################################
class FeatureGroup( object ):
	"""A FeatureGroup is a concrete class that fully defines a family of image features. 

	Its member function CalculateFeatures() is the interface which receives a pixel plane
	as input and returns a vector of doubles which are the image descriptor values 
	derived from that input pixel plane.

	In the WND-CHARM sense, a FeatureGroup is composed of: 1. A list of transforms to apply 
	sequentially to an input image/ROI/pixel plane, and 2. The algorithm to apply to 
	that transformed pixel plane which produces the image descriptor values.

	The member "algorithm" is a reference to the SWIG-wrapped C++ class FeatureAlgorithm.
	The member "transform_list" is a list of references to the SWIG-wrapped C++ class
	FeatureTransform.
	The member "name" is a string representation of the algorithm and transform list,
	following the convention "AlgorithmName ( [Transform A ( [Transform B (] )] )."""

	name = None
	algorithm = None
	transform_list = None
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

	#================================================================
	@classmethod
	def NewFromString( cls, name ):
		"""Takes a string input, parses, and returns an instance of a FeatureGroup class"""

		# The ability to comment out lines with a hashmark at the beginning of the line
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

		return cls( name, Algorithms[ alg ], tform_swig_obj_list )

#############################################################################
# Some more global functions related to Signatures and FeatureGroups:
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
def GenerateWorkOrderFromListOfFeatureStrings( feature_list ):
	"""
	Takes list of feature strings and chops off bin number at the first space on right, e.g.,
	"feature alg (transform()) [bin]" ... Returns a list of FeatureGroups.

	WARNING, RETURNS A TUPLE OF TWO THINGS!

	@return work_order - list of FeatureGroup instances
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
		fg = FeatureGroup.NewFromString( feature_group )
		output_features_count += fg.algorithm.n_features
		work_order.append( fg )

	return work_order, output_features_count

#############################################################################
# class definition of FeatureSet
#############################################################################
class FeatureSet( object ):
	"""An abstract base class inherited by FeatureSet_Discrete (used with FisherFeatureWeights)
	and FeatureSet_Continuous (used with Pearson Corellation Coefficient-weighted 
	ContinuousFeatureWeights).

	An instance of FeatureSet is one-half of a WND-CHARM classifier, the other half being the 
	FeatureWeights instance.

	The FeatureSet class is a container for sets of image descriptors, which are collected
	into Numpy matrices organized by image class or ground truth. FeatureSets are also used
	as containers for test images which have yet to be classified. FeatureSets can also be
	randomly Split() into two subset FeatureSets to be used for cross-validation of a 
	classifier, one for training the other for testing."""

	# source_path - could be the name of a .fit, or pickle file from which this
	# instance was generated, could be a directory
	#  source_path is essentially a name
	# might want to make separate name member in the future
	source_path = None
	num_features = None
	num_images = None

	#: A list of strings length M
	featurenames_list = None

	#: for discrete classes this is a list of lists of image paths
	#: for continuous it's a simple list
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

	# keep track of all the options (-l -S###, etc)
	# FIXME: expand to have all options kept track of individually
	feature_options = None

	#: These members are for discrete training sets, or for continuous training sets
	#: that have some discrete characteristics

	#: A list of strings, length C
	#: FIXME: these two belong in FeatureSet_Discrete
	classnames_list = None
	classsizes_list = None

	# A list of floats against which an image's marginal probaility values can be multiplied
	# to obtain an interpolated value.
	interpolation_coefficients = None
	
	#: A single numpy matrix N features (columns) x M images (rows)
	data_matrix = None
	data_matrix_is_contiguous = False

	#==============================================================
	def __init__( self, data_dict = None ):
		"""FeatureSet constructor"""

		self.featurenames_list = []
		self.imagenames_list = []
		#: FIXME: these two belong in FeatureSet_Discrete
		self.classnames_list = []
		self.classsizes_list = []
		self.interpolation_coefficients = []

		if data_dict != None:
			if "source_path" in data_dict:
				self.source_path = data_dict[ 'source_path' ]
			if "num_classes" in data_dict:
				self.num_classes = data_dict[ 'num_classes' ]
			if "num_features" in data_dict:
				self.num_features = data_dict[ 'num_features' ]
			if "num_images" in data_dict:
				self.num_images = data_dict[ 'num_images' ]
			if "featurenames_list" in data_dict:
				self.featurenames_list = data_dict[ 'featurenames_list' ]
			if "imagenames_list" in data_dict:
				self.imagenames_list = data_dict[ 'imagenames_list' ]
			if "feature_maxima" in data_dict:
				self.feature_maxima = data_dict[ 'feature_maxima' ]
			if "feature_minima" in data_dict:
				self.feature_minima = data_dict[ 'feature_minima' ]
			#: FIXME: these two belong in FeatureSet_Discrete
			if "classnames_list" in data_dict:
				self.classnames_list = data_dict[ 'classnames_list' ]
			if "classsizes_list" in data_dict:
				self.classsizes_list = data_dict[ 'classsizes_list' ]
			if "interpolation_coefficients" in data_dict:
				self.interpolation_coefficients = data_dict[ 'interpolation_coefficients' ]
			if "data_matrix" in data_dict:
				self.data_matrix = data_dict[ 'data_matrix' ]

	#==============================================================
	@output_railroad_switch
	def Print( self ):
		"""Prints out basic attributes about this training set, including name, path to
		source data, number and composition of image classes, number of features, etc."""

		print 'Training Set "{0}"'.format( self.source_path )
		print 'Type: {0}'.format( self.__class__.__name__ )
		print 'Total number of images: {0}'.format( self.num_images )
		print 'Number features: {0}'.format( len( self.featurenames_list ) )

	#==============================================================
	@classmethod
	def NewFromPickleFile( cls, pathname ):
		"""Returns new instance of FeatureSet build from a saved pickle file,
		with a filename ending in .fit.pickle"""

		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".fit.pickled" ):
			raise ValueError( 'Not a pickled FeatureSet file: {0}'.format( pathname ) )

		print "Loading Training Set from pickled file {0}".format( pathname )
		the_training_set = None
		with open( pathname, "rb" ) as pkled_in:
			the_training_set = cls( pickle.load( pkled_in ) )

		# re-generate data_list views from data_matrix and classsizes_list
		if ("data_list" in the_training_set.__dict__):
			the_training_set.data_list = [0] * the_training_set.num_classes
			sample_row = 0
			for i in range( the_training_set.num_classes ):
				nrows = the_training_set.classsizes_list[i]
				the_training_set.data_list[i] = the_training_set.data_matrix[sample_row : sample_row + nrows]
				sample_row += nrows

		# it might already be normalized!
		# FIXME: check for that
		# the_training_set.Normalize()

		return the_training_set

	#==============================================================
	def PickleMe( self, pathname = None ):
		"""Pickle this instance of FeatureSet and write to file whose path is optionally
		specified by argument "pathname" """

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

		# Since we may have both a data_matrix and views into it (data_list), we only want to store
		# one or the other.  Pickle is not smart enough to store numpy views as references.
		# We chose to store the data_matrix, setting data_list to [] if we have it
		# The views have to be reconstructed from the un-pickle using the classsizes_list
		self.ContiguousDataMatrix()
		data_list_copy = None
		if ("data_list" in self.__dict__):
			data_list_copy = self.data_list
			self.data_list = []
		with open( outfile_pathname, 'wb') as outfile:
			pickle.dump( self.__dict__, outfile, pickle.HIGHEST_PROTOCOL )

		# Restore the data_list
		if (data_list_copy):
			self.data_list = data_list_copy

	#==============================================================
	@classmethod
	def NewFromFitFile( cls, pathname ):
		"""Helper function which reads in a c-chrm fit file, builds a dict with the info
		Then calls the constructor and passes the dict as an argument."""
		raise NotImplementedError()

	#==============================================================
	@classmethod
	def NewFromSignature( cls, signature, ts_name = "TestSet", ):
		"""@brief Creates a new FeatureSet from a single signature."""
		raise NotImplementedError()

	#==============================================================
	@classmethod
	def NewFromDirectory( cls, top_level_dir_path, feature_set = "large", write_sig_files_todisk = True ):
		"""The equivalent of the C++ implementation of the command "wndchrm train" """
		raise NotImplementedError()

	#==============================================================
	@classmethod
	def NewFromFileOfFiles( cls, fof_path, options = None ):
		"""FIXME: add ability to specify which features are wanted if none have been calculated yet"""

		if not os.path.exists( fof_path ):
			raise ValueError( "The file '{0}' doesn't exist, maybe you need to specify the full path?".format( fof_path ) )

		print 'Loading {0} from file of files "{1}"'.format( cls.__name__, fof_path )
		
		new_ts = cls()
		new_ts.num_images = 0
		new_ts.source_path = fof_path

		classnames_set = set()

		with open( fof_path ) as fof:
			for line in fof:
				class_id_index = None
				file_path, class_name = line.strip().split( "\t" )
				if not os.path.exists( file_path ):
					raise ValueError(\
					    "The file '{0}' doesn't exist, maybe you need to specify the full path?".\
					    format( file_path ) )
				
				if not class_name in classnames_set:
					classnames_set.add( class_name )
					new_ts.classnames_list.append( class_name )
					class_id_index = len( new_ts.classnames_list ) - 1
				else:
					class_id_index = new_ts.classnames_list.index( class_name )

				if file_path.endswith( (".tif", ".tiff", ".TIF", ".TIFF" ) ): 
					sig = Signatures.NewFromTiffFile( file_path, options )
				elif file_path.endswith( (".sig", "pysig" ) ): 
					sig = Signatures.NewFromSigFile( file_path, options = options )
				else:
					raise ValueError( "File {0} isn't a .tif or a .sig file".format( file_path ) )
				new_ts.AddSignature( sig, class_id_index )

		new_ts.interpolation_coefficients = \
		                    CheckIfClassNamesAreInterpolatable( new_ts.classnames_list )
		
		if isinstance( new_ts, FeatureSet_Discrete ):
			new_ts.num_classes = len( new_ts.data_list )

		new_ts.Print()
		return new_ts
	
	#==============================================================
	@classmethod
	def NewFromSQLiteFile(cls, path):
		"""Not implemented."""
		raise NotImplementedError()

	#==============================================================
	def _ProcessSigCalculationSerially( self, feature_set = "large", write_sig_files_to_disk = True, options = None ):
		"""Virtual method."""
		raise NotImplementedError()

	#==============================================================
	def _ProcessSigCalculationParallelly( self, feature_set = "large", write_sig_files_todisk = True ):
		"""Virtual Method."""
		raise NotImplementedError()

	#==============================================================
	def ContiguousDataMatrix( self ):
		"""This method should be called to access the class data_matrix field.
		In the case where there are both views into this matrix (e.g. data_list) as well as the matrix itself,
		this method ensures that the data_matrix is a vstack of all data_lists
		After this call, the self.data_matrix field will be contiguous, and all of the views in data_list
		are consistent with it.
		In the base class, just returns the data_matrix
		"""
		return (self.data_matrix)

	#==============================================================
	def Normalize( self, training_set = None, quiet = False ):
		"""By convention, the range of values are normalized on an interval [0,100].
		Normalizing is useful in making the variation of features human readable
		and that all samples are comprable if they've been normalized against
		the same ranges of feature values.
		Raw features computed during training are normalized in the training set.
		These training ranges are then used to normalize raw features computed for test images.
		"""

		if self.normalized_against and training_set:
			# I've already been normalized, and you want to normalize me again?
			raise ValueError( "Set {0} has already been normalized against {1}.".format (
				self.source_path, self.normalized_against ) )

		elif not self.normalized_against and not training_set:
			# Normalize me against myself
			if not quiet:
				print 'Normaling set "{0}" ({1} images) against itself'.format (
					self.source_path, self.num_images )

			(self.feature_minima, self.feature_maxima) = normalize_by_columns (self.ContiguousDataMatrix())
			self.normalized_against = "itself"
		else:
			# Normalize me against the given training set
			if training_set.featurenames_list != self.featurenames_list:
				raise ValueError("Can't normalize test_set {0} against training_set {1}: Features don't match.".format (
					self.source_path, training_set.source_path ) )

			if not quiet:
				print 'Normaling set "{0}" ({1} images) against set "{2}" ({3} images)'.format(
					self.source_path, self.num_images, training_set.source_path, training_set.num_images )

			if not training_set.normalized_against:
				training_set.Normalize( )

			assert self.num_features > 0
			(self.feature_minima, self.feature_maxima) = normalize_by_columns (self.ContiguousDataMatrix(), training_set.feature_minima, training_set.feature_maxima)
			self.normalized_against = training_set.source_path
			

	#==============================================================
	def FeatureReduce( self, requested_features ):
		"""Virtual method."""
		raise NotImplementedError()

	#==============================================================
	def AddSignature( self, signature, class_id_index = None ):
		"""Virtual method."""
		raise NotImplementedError()

	#==============================================================
	def ScrambleGroundTruths( self ):
		"""Virtual method. Produce an instant negative control training set"""

	#==============================================================
	def Split( self ):
		"""Virtual method"""
		raise NotImplementedError()

# END FeatureSet class definition

#############################################################################
# class definition of FeatureSet_Discrete
#############################################################################
class FeatureSet_Discrete( FeatureSet ):
	"""One of two (thus far) concrete classes that inherit from the "FeatureSet" base class.

	The difference in structure between "FeatureSet_Discrete" and "FeatureSet_Continuous"
	is that the former has the member "data_list", a Python list, in which the image
	decompositions collected into a list of separate Numpy matrices,
	one for each discrete class. The latter has the members "data_matix" and ground_truths,
	which are single Numpy matrix into which all image descriptors are collected,
	and a list of ground truth values associated with each image, respectively."""

	num_classes = None

	# For C classes, each with Ni images and M features:
	# If the dataset is contiguous, C = 1

	# A list of numpy matrices, length C (one Ni x M matrix for each class)
	# The design is such because it's useful to be able to quickly collect feature statistics
	# across an image class excluding the other classes
	data_list = None

	#==============================================================
	def __init__( self, data_dict = None):
		"""constructor"""
		self.data_list = []
		
		super( FeatureSet_Discrete, self ).__init__( data_dict )

		if data_dict != None:
			if "data_list" in data_dict:
				self.data_list = data_dict[ 'data_list' ]

	#==============================================================
	def Print( self ):
		"""Print basic info about this FeatureSet_Discrete"""
		super( FeatureSet_Discrete, self ).Print()
		class_index = 0
		for class_name in self.classnames_list:
			print '\tClass {0} "{1}": {2} images'.format(
				class_index, class_name, len( self.imagenames_list[ class_index ] ) )
			class_index += 1

	#==============================================================
	def ContiguousDataMatrix( self ):
		"""This method should be called to access the class data_matrix field.
		In the case where there are both views into this matrix (e.g. data_list) as well as the matrix itself,
		this method ensures that the data_matrix is a vstack of all data_lists
		After this call, the self.data_matrix field will be contiguous, and all of the views in data_list
		are consistent with it.
		"""

		# If its already contiguous, or there are no data_lists, just return it
		if (self.data_matrix_is_contiguous or not self.data_list or not len (self.data_list) ):
			return self.data_matrix
		
		num_features = 0
		copy_class = 0
		copy_row = 0
		for class_mat in self.data_list:
			# Make sure all classes have the same number of features
			if (num_features and class_mat.shape[1] != num_features):
				raise ValueError ( "class index {0}:'{1}' has a different number of features than other classes ({3}).".format (
					copy_class, self.classnames_list[i], num_features) )
			else:
				num_features = class_mat.shape[1]
			# if this flag is set in the numpy, then it is not a view.
			if class_mat.flags.owndata:
				# We don't need to keep going, all we need is the startpoint for the copy
				break
			copy_row += class_mat.shape[0]
			copy_class += 1

		# resize the matrix
		if self.data_matrix is not None:
			self.data_matrix.resize (self.num_images, self.num_features)
		else:
			#print "called with empty data_matrix"
			self.data_matrix = np.empty ([ self.num_images, self.num_features ], dtype='double')
			copy_class = 0
			copy_row = 0

		# We need to start copying at the first non-view class mat to the end.
		for class_index in range (copy_class, len (self.data_list)):
			#print "copy class"+str(class_index)
			nrows = self.data_list[class_index].shape[0]
			self.data_matrix[copy_row : copy_row + nrows] = np.copy (self.data_list[class_index])
			self.data_list[class_index] = self.data_matrix[copy_row : copy_row + nrows]
			copy_row += nrows

		self.data_matrix_is_contiguous = True
		return self.data_matrix

	#==============================================================
	@classmethod
	def NewFromFitFile( cls, pathname ):
		"""Helper function which reads in a c-chrm fit file, builds a dict with the info
		Then calls the constructor and passes the dict as an argument"""

		from StringIO import StringIO
		
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
			data_dict[ 'imagenames_list' ] = []
			data_dict[ 'featurenames_list' ] = []
			data_dict[ 'classnames_list' ] = []
			data_dict[ 'classsizes_list' ] = []
			data_dict[ 'imagenames_list' ] = []
			data_dict[ 'data_matrix' ] = None
			data_dict[ 'data_list' ] = []
			tmp_string_data_list = []

			name_line = False
			line_num = 0
			feature_count = 0
			image_pathname = ""
			num_classes = 0
			num_features = 0
			num_images = 0

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
						num_images += 1
					else:
						image_pathname = line.strip()
						data_dict[ 'imagenames_list' ][ zero_indexed_class_id ].append( image_pathname )
					name_line = not name_line
				line_num += 1

		string_data = "\n"
		
		data_matrix = np.empty ([ num_images, num_features ], dtype='double')
		sample_row = 0
		data_dict[ 'data_list' ] = [0] * num_classes
		data_dict[ 'classsizes_list' ] = [0] * num_classes
		for i in range( num_classes ):
			# alt: nrows = len(tmp_string_data_list[i])
			nrows = len(data_dict['imagenames_list'][i])
			print 'generating matrix for class {0} "{1}" ({2} images)'.format(
				i, data_dict['classnames_list'][i], nrows)
			#print "{0}".format( tmp_string_data_list[i] )
			data_matrix[sample_row : sample_row + nrows] = np.genfromtxt( StringIO( string_data.join( tmp_string_data_list[i] ) ) )
			data_dict[ 'data_list' ][i] =  data_matrix[sample_row : sample_row + nrows]
			data_dict[ 'classsizes_list' ][i] = nrows
			sample_row += nrows
		data_dict[ 'data_matrix' ] = data_matrix
		data_dict[ 'interpolation_coefficients' ] = \
		                CheckIfClassNamesAreInterpolatable( data_dict[ 'classnames_list' ] )

		# Instantiate the class
		the_training_set = cls( data_dict )
		
		return the_training_set

	#==============================================================
	@classmethod
	def NewFromSignature( cls, signature, ts_name = "TestSet", ):
		"""@brief Creates a new FeatureSet from a single signature"""

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
		new_ts.classsizes_list.append( 1 )
		new_ts.featurenames_list = signature.names
		new_ts.imagenames_list.append( [ inputimage_filepath ] )
		new_ts.data_matrix = np.array( signature.values )
		data_list.append( new_ts.data_matrix[0:1] )
		data_matrix_is_contiguous = True

		return new_ts

	#==============================================================
	@classmethod
	def NewFromDirectory( cls, top_level_dir_path, feature_set = "large", write_sig_files_todisk = False ):
		"""@brief Equivalent to the "wndchrm train" command from the C++ implementation by Shamir.
		Read the the given directory, parse its structure, and populate the member
		self.imagenames_list. Then call another function to farm out the calculation of each 
		image's features to child processes (or at least it will in the future!)"""

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
		#new_ts.num_images = num_images #taken care of by AddSignatures()
		new_ts.num_images = 0 # initialize
		new_ts.num_classes = num_classes
		new_ts.classnames_list = classnames_list
		#new_ts.imagenames_list = imagenames_list  #taken care of by AddSignatures()
		new_ts.source_path = top_level_dir_path
		new_ts._ProcessSigCalculationSerially( imagenames_list, feature_set, write_sig_files_todisk )
		if feature_set == "large":
			# FIXME: add other options
			new_ts.feature_options = "-l"
		new_ts.interpolation_coefficients = CheckIfClassNamesAreInterpolatable( classnames_list )

		return new_ts


	#==============================================================
	def _ProcessSigCalculationSerially( self, imagenames_list, feature_set = "large", write_sig_files_to_disk = True, options = None ):
		"""Calculate image descriptors for all the images contained in self.imagenames_list"""

		# FIXME: check to see if any .sig, or .pysig files exist that match our
		#        Signature calculation criteria, and if so read them in and incorporate them

		sig = None
		class_id = 0
		for class_filelist in imagenames_list:
			for sourcefile in class_filelist:
				if feature_set == "large":
					sig = Signatures.LargeFeatureSet( sourcefile, options )
				elif feature_set == "small":
					sig = Signatures.SmallFeatureSet( sourcefile, options )
				else:
					raise ValueError( "You requested '{0}' feature set... sig calculation other than small and large feature set hasn't been implemented yet.".format( feature_set ) )
				# FIXME: add all the other options
				# check validity
				if not sig:
					raise ValueError( "Couldn't create a valid signature from file {0} with options {1}".format( sourcefile, options ) )
				sig.is_valid()
				if write_sig_files_to_disk:
					sig.WriteFeaturesToASCIISigFile()
				self.AddSignature( sig, class_id )
			class_id += 1

	#==============================================================
	def FeatureReduce( self, requested_features ):
		"""Returns a new FeatureSet that contains a subset of the features
		arg requested_features is a tuple of features.
		The returned FeatureSet will have features in the same order as they appear in
		requested_features"""

		# FIXME: Roll this function into parent class!!!!!!!!!

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
		reduced_ts = FeatureSet_Discrete()
		reduced_ts.source_path = self.source_path + "(feature reduced)"
		reduced_ts.num_classes = self.num_classes
		assert reduced_ts.num_classes == len( self.data_list )
		new_num_features = len( requested_features )
		reduced_ts.num_features = new_num_features
		reduced_ts.num_images = self.num_images
		reduced_ts.imagenames_list = self.imagenames_list[:] # [:] = deepcopy
		reduced_ts.classnames_list = self.classnames_list[:]
		reduced_ts.classsizes_list = self.classsizes_list[:]
		reduced_ts.featurenames_list = requested_features[:]
		if self.interpolation_coefficients:
			reduced_ts.interpolation_coefficients = self.interpolation_coefficients[:]
		reduced_ts.feature_maxima = np.empty (new_num_features)
		reduced_ts.feature_minima = np.empty (new_num_features)

		# copy features
		reduced_ts.data_matrix = np.empty ([ reduced_ts.num_images, reduced_ts.num_features ], dtype='double')
		new_index = 0
		for featurename in requested_features:
			old_index = self.featurenames_list.index( featurename )
			reduced_ts.data_matrix[:,new_index] = self.data_matrix[:,old_index]
			if self.feature_maxima is not None and self.feature_minima is not None:
				reduced_ts.feature_maxima[ new_index ] = self.feature_maxima[ old_index ]
				reduced_ts.feature_minima[ new_index ] = self.feature_minima[ old_index ]
			new_index += 1

		# regenerate the class views
		reduced_ts._RegenerateClassViews()

		return reduced_ts

	#==============================================================
	def AddSignature( self, signature, class_id_index ):
		""" AddSignature adds signatures quickly to a growing FeatureSet.
		To be quick, the signatures must be added in class index order
		i.e., they are appended to the class-ordered feature matrix, so the class_id_index must be the
		last valid class index, or the one after last.
		To add signatures in random order, use InsertSignature, which is much slower
		@argument signature is a valid signature
		@argument class_id_index identifies the class to which the signature belongs
			class_id_index is zero-indexed
		"""

		if None == class_id_index:
			raise ValueError( 'Must specify either a class_index' )
		
		if( self.data_list == None ) or ( len( self.data_list ) == 0 ) :
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
		while (len( self.data_list ) ) < class_id_index + 1:
			self.data_list.append( None )
		while (len( self.imagenames_list ) ) < class_id_index + 1:
			self.imagenames_list.append( [] )
		while (len( self.classnames_list ) ) < class_id_index + 1:
			self.classnames_list.append( "UNKNOWN"+str(class_id_index + 1) )
		while (len( self.classsizes_list ) ) < class_id_index + 1:
			self.classsizes_list.append( 0 )

		self.imagenames_list[class_id_index].append( signature.source_file )

		if self.data_list[ class_id_index ] == None:
			self.data_list[ class_id_index ] = np.array( signature.values )
		else:
			# vstack takes only one argument, a tuple, thus the extra set of parens
			self.data_list[ class_id_index ] = np.vstack( ( self.data_list[ class_id_index ] ,\
					np.array( signature.values ) ) )

		self.num_images += 1
		self.classsizes_list[class_id_index] += 1

		#print 'Added file "{0}" to class {1} "{2}" ({3} images)'.format( signature.source_file, \
		#    class_id_index, self.classnames_list[class_id_index], len( self.imagenames_list[ class_id_index ] ) ) 

	#==============================================================
	def ScrambleGroundTruths( self, fine_grained_scramble = True ):
		"""Produce an instant negative control training set"""

		import random

		new_ts = self.__class__()
		new_ts.data_list = [ None ] * self.num_classes
		new_ts.imagenames_list = [ [] for j in range( self.num_classes ) ]
		new_ts.num_classes = self.num_classes
		new_ts.classnames_list = self.classnames_list
		new_ts.classsizes_list = self.classsizes_list
		new_ts.featurenames_list = self.featurenames_list
		new_ts.num_features = len( self.featurenames_list )
		new_ts.source_path = self.source_path + " (scrambled)"
		if self.interpolation_coefficients:
			new_ts.interpolation_coefficients = self.interpolation_coefficients

		if fine_grained_scramble:
			# Dump out all the images into one big pool, and reassign them
			# back to the classes randomly
			num_images_per_class = [ len( names_list ) for names_list in self.imagenames_list ]
			num_filenames = sum( num_images_per_class )
			all_filenames = list( itertools.chain.from_iterable( self.imagenames_list ) )
			all_images = np.vstack( self.data_list )

			num_imgs, num_features = all_images.shape
			if not num_imgs == num_filenames:
				raise ValueError( "Number of filenames doesn't match number of signatures" )
			new_ts.num_images = num_imgs

			image_lottery = range( num_imgs )
			random.shuffle( image_lottery )
			image_lottery = iter( image_lottery )

			#for feature_index in range( num_features ):
			#	new_ts.featurenames_list[ feature_index ] += " (scrambled)"

			class_index = 0
			for num_images_in_this_class in num_images_per_class:
				new_matrix = np.empty( ( num_images_in_this_class, num_features ) )
				new_file_list = []

				for i in range( num_images_in_this_class ):
					index = next( image_lottery )
					new_matrix[i] = all_images[ index ]
					new_file_list.append( all_filenames[ index ] )

				new_ts.data_list[ class_index ] = new_matrix
				new_ts.imagenames_list[ class_index ] = new_file_list

				class_index += 1

		else:
			# A course-grained scramble would keep the members in their respective bins
			# but scramble just the class names.
			raise NotImplementedError

		return new_ts

	#==============================================================
	def Split( self, randomize = True, balanced_classes = False, training_set_fraction = None,\
	           i = None, j = None, training_set_only = False, quiet = False ):
		"""
Used for dividing the current FeatureSet into two subsets used for classifier
cross-validation (i.e., training set and test set).

USE CASE LOGIC CHART
--------------------

b = balanced classes
f = training_set_fraction
i = num imgs per class in training set specified
j = num samples in test set specified
o = only training set wanted, skip test set

N = not specified
Y = specified
R = made redundant by other options

b?	f?	i?	j?	o?	outcome:
--	--	--	--	--	--------
N 	N 	N 	N 	N 	unbalanced training and test sets, default 75-25 split in each class, i's & j's different for each class
N 	Y 	R 	R 	N 	unbalanced training and test sets, 75-25 default replaced, i's & j's different for each class
N 	R 	Y 	N 	N 	balanced training set, detritus goes into the unbalanced test set
N 	R 	N 	Y 	N 	balanced test set, detritus goes into the unbalanced training set
R 	R 	Y 	Y 	N 	balanced training and test sets as specified by i & j
N 	N 	N 	R 	Y 	unbalanced training set only, num samples in each class determined by fraction
R 	R 	Y 	R 	Y 	balanced training set only, num samples specified by i
Y 	N 	N 	N 	N 	balanced training and test sets, i=0.75*smallest class, j=0.25*smallest class
Y 	Y 	R 	R 	N 	balanced training and test sets, as above but with 75-25 default replaced

		"""

		# FIXME: use np.random.shuffle(arr) - shuffles first dimension (rows) of multi-D numpy, so images in our case.
		# If the class views are shuffled one-by-one, then the main matrix will be shuffled as well, but only within classes.
		# then, just do the split by slicing the classes based on train/test sizes
		# see also http://stackoverflow.com/questions/4601373/better-way-to-shuffle-two-numpy-arrays-in-unison
		# for shuffling multiple arrays in unison
		# 		def shuffle_in_unison(a, b):
		# 			rng_state = numpy.random.get_state()
		# 			numpy.random.shuffle(a)
		# 			numpy.random.set_state(rng_state)
		# 			numpy.random.shuffle(b)
		# or, using a permuted vector as an index
		#	p = numpy.random.permuation(len(a))
		#	return a[p], b[p]

		# Shuffle in unison perhaps not the best way since it may not be necessary to shuffle
		# all samples in a class? Also how is p = numpy.random.permuation(len(a)) different from
		# currently used p = random.shuffle(range(a))? - CEC


		# Step 0: General Dummyproofing
		smallest_class_size = self.NumSamplesInSmallestClass()

		if i and( i <= 0 or i > smallest_class ):
			raise ValueError( 'i must be greater than zero and less than total number of images'\
			    + ' in smallest class ({0})'.format( smallest_class ) )

		if j and( j <= 0 or j > smallest_class ):
			raise ValueError( 'j must be greater than zero and less than total number of images'\
			    + ' in smallest class ({0})'.format( smallest_class ) )

		if ( i and j ) and ( ( i + j ) > smallest_class ):
			raise ValueError( 'Values for i and j cannot add up to more than total number of images in smallest class ({0})'.format( smallest_class ) )
		
		if training_set_fraction and ( training_set_fraction < 0 or training_set_fraction > 1 ):
			raise ValueError( "Training set fraction must be a number between 0 and 1" )

		if j and ( j <= 0 ) and not training_set_only == True:
			raise UserWarning( "j value of 0 implies only training set is desired" )
			training_set_only = True

		if training_set_fraction and ( training_set_fraction == 1.0 ) and not training_set_only == True:
			raise UserWarning( "training set fraction value of 1 implies only training set is desired" )
			training_set_only = True

		if i and j and training_set_fraction:
			raise ValueError( 'Conflicting input: You specified i, j and training_set fraction, which is redundant.' )

		if j and ( j > 0 ):
			raise ValueError( 'Conflicting input: You specified both a non-zero value for j, but also training_set_only set to true.')

		# Specify defaults here instead of in method argument for the purpose of dummyproofing
		if not training_set_fraction:
			training_set_fraction = 0.75

		# Step 1: Number of samples in training set
		num_samples_in_training_set = None
		if i:
			num_samples_in_training_set = [ i ] * self.num_classes
		elif balanced_classes and j:
			num_samples_in_training_set = [ smallest_class_size - j ] * self.num_classes
		elif balanced_classes and training_set_fraction:
			num_samples_in_training_set = \
			       [ int(round( training_set_fraction * smallest_class_size )) ] * self.num_classes
		elif j:
			# you want detritus to go into training set
			num_samples_in_training_set = [ ( num - j ) for num in self.classsizes_list ]
		else:
			# unbalanced
			num_samples_in_training_set = [ int(round( training_set_fraction * num )) for num in self.classsizes_list ]

		# Step 2: Number of samples in test set
		if not training_set_only:
			num_samples_in_test_set = None
			if j:
				num_samples_in_test_set = [ j ] * self.num_classes
			elif balanced_classes and i:
				num_samples_in_test_set = [ smallest_class_size - i ] * self.num_classes
			elif balanced_classes and training_set_fraction:
				num_samples_in_test_set = [ (smallest_class_size - num) for num in num_samples_in_training_set ]
			elif i:
				# you want detritus to go into test set
				num_samples_in_test_set = [ ( num - i ) for num in self.classsizes_list ]
			else:
				# you want an unbalanced test set
				num_samples_in_test_set = [ int(round( (1.0-training_set_fraction) * num )) for num in self.classsizes_list ]


		# Say what we're gonna do:
		if not quiet:
			print "\t" + "\t".join( self.classnames_list ) + "\ttotal:"
			print "Train Set\t" + "\t".join( [ str(num) for num in num_samples_in_training_set ] ) + "\t{0}".format( sum( num_samples_in_training_set ) )
			if not training_set_only:
				print "Test Set\t" + "\t".join( [ str(num) for num in num_samples_in_test_set ] ) + "\t{0}".format( sum( num_samples_in_test_set ) )
			if randomize:
				print "Sample membership chosen at random."

		if randomize:
			import random

		training_set = None
		test_set = None

		training_set = self.__class__()
		training_set.data_list = [ None ] * self.num_classes
		training_set.num_images = 0
		training_set.num_classes = self.num_classes
		training_set.classnames_list = self.classnames_list
		training_set.classsizes_list = num_samples_in_training_set
		training_set.featurenames_list = self.featurenames_list
		training_set.num_features = len( self.featurenames_list )
		training_set.imagenames_list = [ [] for j in range( self.num_classes ) ]
		training_set.source_path = self.source_path + " (subset)"
		if self.interpolation_coefficients:
			training_set.interpolation_coefficients = self.interpolation_coefficients
	
		if not training_set_only:
			test_set = self.__class__()
			test_set.data_list = [ None ] * self.num_classes
			test_set.num_images = 0
			test_set.num_classes = self.num_classes
			test_set.classnames_list = self.classnames_list
			test_set.classsizes_list = num_samples_in_test_set
			test_set.featurenames_list = self.featurenames_list
			test_set.num_features = len( self.featurenames_list )
			test_set.imagenames_list = [ [] for j in range( self.num_classes ) ]
			test_set.source_path = self.source_path + " (subset)"
			if self.interpolation_coefficients:
				test_set.interpolation_coefficients = self.interpolation_coefficients

		# assemble training and test sets
		for class_index in range( self.num_classes ):

			# If randomize, choose samples at random, but once they're chosen, pack them into
			# the FeatureSets in alphanumeric order.
			sort_func = lambda A, B: cmp( self.imagenames_list[ class_index ][ A ], self.imagenames_list[ class_index ][ B ] )

			num_images = self.classsizes_list[ class_index ]
			sample_lottery = range( num_images )
			if randomize:
				random.shuffle( sample_lottery )

			# training set:
			training_matrix = np.empty( [ training_set.classsizes_list[ class_index ], self.num_features ], dtype='double' )

			training_samples = sample_lottery[ 0 : training_set.classsizes_list[ class_index ] ]
			training_samples = sorted( training_samples, sort_func )

			train_samp_count = 0
			for sample_index in training_samples:
				sample_name = self.imagenames_list[ class_index ][ sample_index ]
				training_matrix[ train_samp_count,: ] = self.data_list[ class_index ][ sample_index ]
				training_set.imagenames_list[ class_index ].append( sample_name )
				training_set.num_images += 1
				train_samp_count += 1
			training_set.data_list[ class_index ] = training_matrix

			# test set:
			if not training_set_only:
				test_samp_count = 0
				test_matrix = np.empty( [ test_set.classsizes_list[ class_index ], self.num_features ], dtype='double' )

				test_samples = sample_lottery[ training_set.classsizes_list[ class_index ] : \
				   training_set.classsizes_list[ class_index ] + test_set.classsizes_list[ class_index ] ]
				test_samples = sorted( test_samples, sort_func )

				for sample_index in test_samples:
					sample_name = self.imagenames_list[ class_index ][ sample_index ]
					test_matrix[ test_samp_count,: ] = self.data_list[ class_index ][ sample_index ]
					test_set.imagenames_list[ class_index ].append( sample_name )
					test_set.num_images += 1
					test_samp_count += 1
				test_set.data_list[ class_index ] = test_matrix

		if training_set_only:
			return training_set

		return training_set, test_set

	#==============================================================
	def NumSamplesInSmallestClass( self ):
		"""Method name says it all."""
		return np.amin( np.array( [ len(samplenames) for samplenames in self.imagenames_list ] ) )
	#==============================================================
	def _RegenerateClassViews( self ):
		"""Rebuild the views in self.data_list for convenient access."""

		sample_row = 0
		self.data_list = [0] * self.num_classes
		for class_index in range( self.num_classes ):
			nrows = self.classsizes_list[ class_index ]
			self.data_list[class_index] = self.data_matrix[sample_row : sample_row + nrows]
			sample_row += nrows

# END FeatureSet_Discrete class definition


#############################################################################
# class definition of FeatureSet_Continuous
#############################################################################
class FeatureSet_Continuous( FeatureSet ):
	"""One of two (thus far) concrete classes that inherit from the "FeatureSet" base class.

	The difference in structure between "FeatureSet_Discrete" and "FeatureSet_Continuous"
	is that the former has the member "data_list", a Python list, in which the image
	decompositions collected into a list of separate Numpy matrices,
	one for each discrete class. The latter has the members "data_matix" and ground_truths,
	which are single Numpy matrix into which all image descriptors are collected,
	and a list of ground truth values associated with each image, respectively."""

	#: Ground truth numerical values accociated with each image
	ground_truths = None

	def __init__( self, data_dict = None ):

		# call parent constructor
		self.ground_truths = []

		super( FeatureSet_Continuous, self ).__init__( data_dict )

		if data_dict != None:
			if "data_matrix" in data_dict:
				self.data_matrix = data_dict[ 'data_matrix' ]
			if "ground_truths" in data_dict:
				self.ground_truths = data_dict[ 'ground_truths' ]

	#==============================================================
	def Print( self ):
		"""Calls parent class method"""
		super( FeatureSet_Continuous, self ).Print()

	#==============================================================
	@classmethod
	def NewFromFitFile( cls, pathname ):
		"""Construct a FeatureSet_Continuous from a ".fit" training set file created by 
		the C++ implentation of WND-CHARM."""

		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".fit" ):
			raise ValueError( 'Not a .fit file: {0}'.format( pathname ) )

		pickled_pathname = pathname + ".pychrm"

		print "Creating Continuous Training Set from legacy WND-CHARM text file file {0}".format( pathname )
		with open( pathname ) as fitfile:
			new_ts = cls()

			new_ts.source_path = pathname
			tmp_string_data_list = []

			name_line = False
			line_num = 0
			feature_count = 0
			image_pathname = ""
			num_classes = 0
			num_features = 0

			import re

			for line in fitfile:
				if line_num is 0:
					num_classes = int( line )
				elif line_num is 1:
					num_features = int( line )
					new_ts.num_features = num_features
				elif line_num is 2:
					new_ts.num_images = int( line )
				elif line_num <= ( num_features + 2 ):
					new_ts.featurenames_list.append( line.strip() )
					feature_count += 1
				elif line_num == ( num_features + 3 ):
					pass # skip a line
				elif line_num <= ( num_features + 3 + num_classes ):
					line = line.strip()
					new_ts.classnames_list.append( line )
					new_ts.classsizes_list.append( 0 )
					m = re.search( r'(\d*\.?\d+)', line )
					if m:
						new_ts.interpolation_coefficients.append( float( m.group(1) ) )
					else:
						raise ValueError( "Can't create continuous training set, one of the class names " \
								"'{0}' is not able to be interpreted as a number.".format( line ) )
				else:
					# Read in features
					# Comes in alternating lines of data, then tile name
					if not name_line:
						# strip off the class identity value, which is the last in the array
						split_line = line.strip().rsplit( " ", 1)
						zero_indexed_class_id = int( split_line[1] ) - 1
						tmp_string_data_list.append( split_line[0] )
						new_ts.ground_truths.append( new_ts.interpolation_coefficients[ zero_indexed_class_id ] )
						new_ts.classsizes_list[ zero_indexed_class_id ] += 1
					else:
						image_pathname = line.strip()
						new_ts.imagenames_list.append( image_pathname )
					name_line = not name_line
				line_num += 1

		string_data = "\n"
		print "parsing text into a numpy matrix"

		from StringIO import StringIO
		new_ts.data_matrix = np.genfromtxt( StringIO( string_data.join( tmp_string_data_list ) ) )
		
		return new_ts

	#==============================================================
	def _ProcessSigCalculationSerially( self, imagenames_list, feature_set = "large", write_sig_files_to_disk = True, options = None ):
		"""Begin calculation of image features for the images listed in self.imagenames_list,
		but do not use parallel computing/ create child processes to do so.
		
		This function is to be deprecated ASAP in favor of parallel feature calculation"""
		# FIXME: check to see if any .sig, or .pysig files exist that match our
		#        Signature calculation criteria, and if so read them in and incorporate them

		sig = None
		class_id = 0
		for class_filelist in imagenames_list:
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
			

	#==============================================================
	def FeatureReduce( self, requested_features ):
		"""Returns a new FeatureSet that contains a subset of the features
		arg requested_features is a tuple of features
		the returned FeatureSet will have features in the same order as they appear in
		     requested_features"""

		# FIXME: Roll this function into parent class
		# Check that self's faturelist contains all the features in requested_features

		selfs_features = set( self.featurenames_list )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = error_banner + "Feature Reduction error:\n"
			err_str += "The training set '{0}' is missing ".format( self.source_path )
			err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
					len( missing_features_from_req ), len( requested_features ) )
			err_str += "\nDid you forget to convert the feature names from the C++ implementation style into their modern counterparts using FeatureNameMap.TranslateToNewStyle()?"

			raise ValueError( err_str )

		# copy everything but the signature data
		reduced_ts = FeatureSet_Continuous()
		new_num_features = len( requested_features )
		reduced_ts.source_path = self.source_path + "({0} features)".format( new_num_features )
		reduced_ts.num_features = new_num_features
		reduced_ts.num_images = self.num_images
		reduced_ts.imagenames_list = self.imagenames_list[:] # [:] = deepcopy
		reduced_ts.featurenames_list = requested_features[:]
		if self.interpolation_coefficients:
			reduced_ts.interpolation_coefficients= self.interpolation_coefficients[:]
		if self.classnames_list:
			reduced_ts.classnames_list = self.classnames_list[:]
		if self.classsizes_list:
			reduced_ts.classsizes_list = self.classsizes_list[:]
		if self.ground_truths:
			reduced_ts.ground_truths = self.ground_truths[:]
		reduced_ts.feature_maxima = np.empty (new_num_features)
		reduced_ts.feature_minima = np.empty (new_num_features)

		# copy feature minima/maxima
		# Have to explicitly check for "is not None" due to exception:
		# ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
		if (not self.feature_maxima == None) and (not self.feature_minima == None):
			new_index = 0
			for featurename in requested_features:
				old_index = self.featurenames_list.index( featurename )
				reduced_ts.feature_maxima[ new_index ] = self.feature_maxima[ old_index ]
				reduced_ts.feature_minima[ new_index ] = self.feature_minima[ old_index ]
				new_index += 1

		# feature reduce
		num_imgs_in_class, num_old_features = self.data_matrix.shape
		# NB: double parentheses required when calling numpy.zeros(), i guess it's a tuple thing
		reduced_ts.data_matrix = np.zeros( ( num_imgs_in_class, new_num_features ) )
		new_column_index = 0
		for featurename in requested_features:
			fat_column_index = self.featurenames_list.index( featurename )
			reduced_ts.data_matrix[:,new_column_index] = self.data_matrix[:,fat_column_index]
			new_column_index += 1

		return reduced_ts

	#==============================================================
	def AddSignature( self, signature, ground_truth = None ):
		"""@argument signature is a valid signature
		@argument class_id_index identifies the class to which the signature belongs"""
		
		if (self.data_matrix == None) or ( len( self.data_matrix ) == 0 ) :
			self.data_matrix = np.array( signature.values )
			self.featurenames_list = signature.names
			self.num_features = len( signature.names )
		else:
			# Make sure features are in order.
			# Feature Reduce will throw an exception if they can't be placed in order
			signature = signature.FeatureReduce( self.featurenames_list )
			# vstack takes only one argument, a tuple, thus the extra set of parens
			self.data_matrix = np.vstack( (self.data_matrix, np.array( signature.values )) )

		self.num_images += 1
		self.imagenames_list.append( signature.source_file )

		self.ground_truths.append( ground_truth )

	#==============================================================
	def ScrambleGroundTruths( self ):
		"""Produce an instant negative control training set"""

		import random
		random.shuffle( self.ground_truths )
		self.source_path += " (scrambled)"

	#==============================================================
	def Split( self, randomize = True, training_set_fraction = 0.75,\
	           i = None, j = None, training_set_only = False, quiet = False ):
		"""Used for dividing the current FeatureSet into two subsets used for classifier
		cross-validation (i.e., training set and test set).
		
		Number of images in training and test sets are allocated by i and j, respectively
		Otherwise they are given by training_set fraction.
		"""

		# FIXME: use np.random.shuffle(arr) - shuffles first dimension (rows) of multi-D numpy, so images in our case.
		# see also http://stackoverflow.com/questions/4601373/better-way-to-shuffle-two-numpy-arrays-in-unison
		# for shuffling multiple arrays in unison
		# 		def shuffle_in_unison(a, b):
		# 			rng_state = numpy.random.get_state()
		# 			numpy.random.shuffle(a)
		# 			numpy.random.set_state(rng_state)
		# 			numpy.random.shuffle(b)
		# or, using a permuted vector as an index
		#	p = numpy.random.permuation(len(a))
		#	return a[p], b[p]
		# Figure out how many images will be in which class
		if i and j:
			if (i + j) > self.num_images:
				raise ValueError( 'Values for i and j cannot add up to more than total number of images in parent set ({0})'.format( self.num_images ) )

		if i:
			if i <= 0 or i > self.num_images:
				raise ValueError( 'i must be greater than zero and less than total number of images'\
				    + ' in parent set ({0})'.format( self.num_images ) )
				num_images_in_training_set = i
		else:
			if training_set_fraction <= 0 or training_set_fraction >= 1:
				raise ValueError( "Training set fraction must be a number between 0 and 1" )
			num_images_in_training_set = int( round( training_set_fraction * self.num_images ) )

		if j:
			if j <= 0:
				training_set_only = True
			elif j > ( self.num_images - num_images_in_training_set ):
				raise ValueError( 'j must be less than total number of images in parent set ({0}) minus # images in training set ({1})'.format( self.num_images, num_images_in_training_set ) )
			training_set_only = False
			num_images_in_test_set = j
		else:
			num_images_in_test_set = self.num_images - num_images_in_training_set

		# Say what we're gonna do:
		if not quiet:
			out_str = ''
			if randomize:
				out_str += 'Randomly splitting '
			else:
				out_str += 'Splitting '
			out_str += '{0} "{1}" ({2} images) into '.format( type( self ).__name__, \
				 self.source_path, self.num_images )
			out_str += "training set ({0} images)".format( num_images_in_training_set )
			if not training_set_only:
				out_str += " and test set ({0} images)".format( num_images_in_test_set )
			print out_str

		# initialize everything
		training_set = None
		test_set = None
		training_set = self.__class__()
		training_set.num_images = num_images_in_training_set
		training_set.featurenames_list = self.featurenames_list
		training_set.num_features = len( self.featurenames_list )
		training_set.imagenames_list = []
		training_set.source_path = self.source_path + " (subset)"
		if self.classnames_list:
			training_set.classnames_list = self.classnames_list
		if self.classsizes_list:
			training_set.classsizes_list = [ 0 ] * self.num_classes
		if self.interpolation_coefficients:
			training_set.interpolation_coefficients = self.interpolation_coefficients
	
		if not training_set_only:
			test_set = self.__class__()
			test_set.num_images = num_images_in_test_set
			test_set.featurenames_list = self.featurenames_list
			test_set.num_features = len( self.featurenames_list )
			test_set.imagenames_list = []
			test_set.source_path = self.source_path + " (subset)"
			if self.classnames_list:
				test_set.classnames_list = self.classnames_list
			if self.classsizes_list:
				test_set.classsizes_list = [ 0 ] * self.num_classes
			if self.interpolation_coefficients:
				test_set.interpolation_coefficients = self.interpolation_coefficients

		image_lottery = range( self.num_images )
		if randomize:
			import random
			random.shuffle( image_lottery )

		train_image_lottery = image_lottery[:num_images_in_training_set]
		test_image_lottery = image_lottery[num_images_in_training_set: \
		                                   num_images_in_training_set + num_images_in_test_set ]

		# build the training and test sets such that their ground truths are sorted in 
		# ascending order.
		sort_func = lambda A, B: cmp( self.ground_truths[A], self.ground_truths[B] ) if self.ground_truths[A] != self.ground_truths[B] else cmp( self.imagenames_list[A], self.imagenames_list[B] )
		# I guess we really don't have to sort the training set, do we?
		# train_image_lottery = sorted( train_image_lottery, sort_func )
		test_image_lottery = sorted( test_image_lottery, sort_func )

		# Training Set first
		train_image_count = 0
		training_matrix = np.empty( ( num_images_in_training_set, len( self.featurenames_list ) ) ) 
		while train_image_count < num_images_in_training_set:
			image_index = train_image_lottery[ train_image_count ]
			image_name = self.imagenames_list[ image_index ]
			training_set.imagenames_list.append( image_name )
			training_matrix[ train_image_count ] = self.data_matrix[ image_index ]
			training_set.ground_truths.append( self.ground_truths[ image_index ] )
			train_image_count += 1
		training_set.data_matrix = training_matrix

		# Now Test Set, if applicable
		test_matrix = None
		if not training_set_only:
			test_matrix = np.empty( ( num_images_in_test_set, len( self.featurenames_list ) ) )
			test_image_count = 0
			while test_image_count < num_images_in_test_set:
				image_index = test_image_lottery[ test_image_count ]
				image_name = self.imagenames_list[ image_index ]
				test_set.imagenames_list.append( image_name )
				test_matrix[ test_image_count ] = self.data_matrix[ image_index ]
				test_set.ground_truths.append( self.ground_truths[ image_index ] )
				test_image_count += 1
			test_set.data_matrix = test_matrix
		
		return training_set, test_set

# END FeatureSet_Continuous class definition

#=================================================================================
class ClassificationResult( object ):
	"""Interface class which all Result classes inherit"""

	def GenerateStats( self ):
		"""Virtual method"""
		raise NotImplementedError
	
	@output_railroad_switch
	def Print( self ):
		"""Virtual method"""
		raise NotImplementedError

	def GenerateHTML( self ):
		"""Virtual method"""
		raise NotImplementedError


#=================================================================================
class ImageClassificationResult( ClassificationResult ):
	"""Abstract base class that is meant to contain the information from a single
	classification of a single image/ROI."""

	name = None
	source_file = None
	ground_truth_value = None
	predicted_value = None
	batch_number = None

	#: For the future: a member to indicate the position of the ROI within an image
	kernel_location = None

	@output_railroad_switch
	def Print():
		raise NotImplementedError

#=================================================================================
class DiscreteImageClassificationResult( ImageClassificationResult ):
	"""Concrete class that contains the result of a classification for a single image/ROI.
	which includes predicted class, marginal probabilities, etc."""

	normalization_factor = None
	marginal_probabilities = None
	#: predicted_class_name will always be a string
	#: the interpolated value, if applicable, gets stored in self.predicted_vlaue
	predicted_class_name = None
	ground_truth_class_name = None

	def __init__( self ):
		self.marginal_probabilities = []

	#==============================================================
	@output_railroad_switch
	def Print( self, line_item = False ):
		"""Output classification line item data, including predicted class and marginal
		probabilities"""
		
		if line_item:
			# img name:
			output_str = self.source_file
			# normalization factor:
			output_str += "\t{val:0.3g}\t".format( val=self.normalization_factor )
			# marginal probabilities:
			output_str += "\t".join(\
					[ "{val:0.3f}".format( val=prob ) for prob in self.marginal_probabilities ] )
			output_str += "\t"
			# actual class:
			if self.ground_truth_class_name:
				output_str += "{0}\t".format( self.ground_truth_class_name )
			else:
				output_str += "*\t"
			# predicted class:
			output_str += self.predicted_class_name + "\t"
			# interpolated value, if applicable
			if self.predicted_value:
				output_str += "{val:0.3f}".format( val=self.predicted_value )
			print output_str
		else:
			print "Image:             \t{0}".format( self.source_file )
			print "Normalization Factor:\t{0}".format( self.normalization_factor )
			print "Marg. Probabilities:\t" + "\t".join(\
					[ "{val:0.3f}".format( val=prob ) for prob in self.marginal_probabilities ] )
			print "Ground truth class:\t {0}".format( self.ground_truth_class_name ) 
			print "Predicted class:\t {0}".format( self.predicted_class_name ) 
			if self.predicted_value:
				print "Interpolated Value:\t{0}".format( self.predicted_value )
				#print "Interpolated Value:\t{val=0.3f}".format( val=self.predicted_value )

	#==============================================================
	@classmethod
	def _WND5( cls, trainingset, testimg, feature_weights ):
		"""
		Don't call this function directly, use the wrapper functions 
		DiscreteBatchClassificationResult.New()  (for test sets) or
		DiscreteImageClassificationResult.NewWND5() (for single images)
		Both of these functions have dummyproofing.

		If you're using this function, your training set data is not continuous
		for N images and M features:
			trainingset is list of length L of N x M numpy matrices
			testtile is a 1 x M list of feature values
		NOTE: the trainingset and test image must have the same number of features!!!
		AND: the features must be in the same order!!
		Returns an instance of the class DiscreteImageClassificationResult
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
				# epsilon checking for each feature is too expensive
				# FIXME: Do quick and dirty summation check until we can figure something else out
				dists = np.absolute( sig_matrix[ tile_index ] - testimg )
				w_dist = np.sum( dists )
#				print "train img {0} dist : {1}".format( tile_index, w_dist )
#				if (np.isinf(w_dist)):
#					print "dists: "+str(dists)
				if w_dist < epsilon:
					num_collisions += 1
					continue
				dists = np.multiply( weights_squared, np.square( dists ) )
				w_dist = np.sum( dists )
				# The exponent -5 is the "5" in "WND5"
				class_similarities[ class_index ] += w_dist ** -5
			#print "\n"

			class_similarities[ class_index ] /= ( num_tiles - num_collisions )
#			print "class_similarities: "+str(class_similarities)

		result = cls()
		norm_factor = sum( class_similarities )
		result.normalization_factor = norm_factor 
		result.marginal_probabilities = [ x / norm_factor for x in class_similarities ]

		return result

	#=================================================================================
	@classmethod
	def NewWND5( cls, training_set, feature_weights, test_sig, quiet = False ):
		"""@brief: A wrapper function for _ClassifyOneImageWND5 that does dummyproofing
		@return: An instance of a DiscreteBatchClassificationResult"""

		if not isinstance( training_set, FeatureSet_Discrete ):
			raise ValueError( 'First argument to NewWND5 must be of type "FeatureSet_Discrete", you gave a {0}'.format( type( training_set ).__name__ ) )
		
		if not isinstance( feature_weights, FeatureWeights ):
			raise ValueError( 'Second argument to NewWND5 must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )

		if not isinstance( test_sig, Signatures ):
			raise ValueError( 'Third argument to NewWND5 must be of type "Signatures", you gave a {0}'.format( type( test_sig ).__name__ ) )

		# check to see if sig is valid
		test_sig.is_valid()

		train_set_len = len( training_set.featurenames_list )
		test_set_len = len( test_sig.names )
		feature_weights_len = len( feature_weights.names )

		if test_sig.names != feature_weights.names:
			raise ValueError("Can't classify, features in signature don't match features in weights." )

		if test_sig.names != training_set.featurenames_list:
			raise ValueError("Can't classify, features in signature don't match features in training_set." )

		print "Classifying image '{0}' ({1} features) against test set '{2}' ({3} features)".\
			 format( test_sig.source_file, train_set_len, training_set.source_path, test_set_len )

		result = cls._WND5( training_set, test_sig.values, feature_weights.values )

		result.source_file = test_sig.source_file
		marg_probs = np.array( result.marginal_probabilities )
		result.predicted_class_name = training_set.classnames_list[ marg_probs.argmax() ]
		# interpolated value, if applicable
		if len (training_set.interpolation_coefficients) == len (marg_probs):
			interp_val = np.sum( marg_probs * training_set.interpolation_coefficients )
			result.predicted_value = interp_val

		if not quiet:
			column_header = "image\tnorm. fact.\t"
			column_header +=\
			 "".join( [ "p(" + class_name + ")\t" for class_name in training_set.classnames_list ] )
			column_header += "act. class\tpred. class\tpred. val."
			print column_header
			result.Print( line_item = True )

		return result


#=================================================================================
class ContinuousImageClassificationResult( ImageClassificationResult ):
	"""Concrete class that contains the result of a classification for a single image/ROI.
	The primary result is a predicted value for the given image based on the regression
	parameters obtained when generating a continuous classifier."""

	#==============================================================
	@output_railroad_switch
	def Print( self, line_item = False ):
		"""Output results."""

		if line_item:
			# img name:
			output_str = self.source_file
			output_str += "\t"
			# actual class:
			if self.ground_truth_value is not None:
				output_str += str( self.ground_truth_value ) + "\t"
			else:
				output_str += "*\t"
			# predicted class:
			output_str += str( self.predicted_value )
			print output_str
		else:
			print "Image:             \t{0}".format( self.source_file )
			print "Ground truth class:\t {0}".format( self.ground_truth_value ) 
			print "Predicted class:\t {0}".format( self.predicted_value ) 

	#==============================================================
	@classmethod
	def _LinearRegression( cls, one_image_features, feature_weights ):
		"""Produce a predicted value for a single image based on the regression parameters
		contained in the ContinuousFeatureWeights argument "feature_weights".
		
		Don't call this function directly, but instead use the member function
		New() on the ContinuousBatchClassificationResult class which has dummyproofing."""

		per_feature_predicted_vals = []

		for i in range( len( one_image_features ) ):
			feat_val = one_image_features[i]
			weight = feature_weights.values[i]
			slope = feature_weights.slopes[i]
			intercept = feature_weights.intercepts[i]

			# y = mx+b
			# feature value = slope * age + intercept
			# solve for the abscissa:
			per_feature_predicted_vals.append( weight * ( feat_val - intercept ) / slope )

		result = cls()
		result.predicted_value = sum( per_feature_predicted_vals )

		return result

#=================================================================================
class BatchClassificationResult( ClassificationResult ):
	"""An abstract base class which serves as container for individual 
	ImageClassificationResult instances.
	
	The equivalent concept to a BatchClassificationResult in the C++ implementation
	of WND-CHARM is the split, i.e., command line arg of -n10 would generate 10 
	train/test splits."""

	name = None
	training_set = None
	test_set = None
	feature_weights = None
	figure_of_merit = None
	individual_results = None
	predicted_values = None
	ground_truth_values = None

	num_classifications = None


	# This stuff is for correllation analysis
	pearson_coeff = None
	pearson_p_value = None
	pearson_std_err = None
	spearman_coeff = None
	spearman_p_value = None
	ground_truth_values = None
	predicted_values = None


	#==============================================================
	def __init__( self, training_set = None, test_set = None, feature_weights = None ):
		"""BatchClassificationResult constructor"""

		self.training_set = training_set
		self.test_set = test_set
		self.feature_weights = feature_weights
		self.individual_results = []

		self.num_classifications = 0


	#==============================================================
	def GenerateStats( self ):
		"""Base method that's called by daughter classes.
		
		If this BatchClassificationResult contains ImageClassificationResults that
		has known ground truth as well as predicted values, this method will calculate
		statistics describing how well predicted values correlate with ground truth.

		Requires scipy.stats package to be installed"""
		#FIXME: Calculate p-value of our standard error figure of merit.

		if len( self.individual_results ) == 0:
			raise ValueError( 'No individual results in this batch with which to generate statistics.' )

		self.num_classifications = len( self.individual_results )

		if self.ground_truth_values and self.predicted_values and \
		        len( self.ground_truth_values ) == len( self.predicted_values ):

			gt = np.array( self.ground_truth_values )
			pv = np.array( self.predicted_values )

			diffs = gt - pv
			diffs = np.square( diffs )
			err_sum = np.sum( diffs )

			import math; from scipy import stats
			self.figure_of_merit = math.sqrt( err_sum / self.num_classifications )

			# For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
			np.seterr (under='ignore')
			slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
			             stats.linregress( self.ground_truth_values, self.predicted_values )

			self.spearman_coeff, self.spearman_p_value =\
			       stats.spearmanr( self.ground_truth_values, self.predicted_values )
			np.seterr (all='raise')


	#==============================================================
	def RankOrderSort( self ):
		"""Rank-order sorts ground truth/predicted value data points for purposes of 
		being graphed."""

		value_pairs = zip( self.ground_truth_values, self.predicted_values )

		# sort by ground_truth value first, predicted value second
		sort_func = lambda A, B: cmp( A[0], B[0] ) if A[0] != B[0] else cmp( A[1], B[1] ) 

		sorted_pairs = sorted( value_pairs, sort_func )
		
		# we want lists, not tuples!
		self.ground_truth_values, self.predicted_values =\
			[ list( unzipped_tuple ) for unzipped_tuple in zip( *sorted_pairs ) ]	

	#==============================================================
	def PickleMe( self ):
		"""Not Implemented"""
		raise NotImplementedError
	def NewFromPickleFile( self ):
		"""Not Implemented"""
		raise NotImplementedError


#=================================================================================
class DiscreteBatchClassificationResult( BatchClassificationResult ):
	"""Concrete class that is a container for DiscreteImageClassificationResult instances.
	Use this class to run classification experiments and to generate statistics on their results.
	For discrete (Fisher) classifiers, the figure_of_merit is classification accuracy."""

	num_correct_classifications = None
	classification_accuracy = None
	
	num_classifications_per_class = None
	num_correct_classifications_per_class = None

	confusion_matrix = None
	similarity_matrix = None
	average_class_probability_matrix = None

	#==============================================================
	def __init__( self, training_set = None, test_set = None, feature_weights = None ):
		"""Simply calls parent constructor."""
		super( DiscreteBatchClassificationResult, self ).__init__( training_set, test_set, feature_weights)

	#==============================================================
	def GenerateStats( self ):
		"""Fills out the confusion matrix, similarity matrix, and class probability matrix."""

		# Run a standard error analysis if ground_truth/predicted vals exist (call base method ):
		# This also sets self.num_classifications
		super( DiscreteBatchClassificationResult, self ).GenerateStats()

		num_classes = self.training_set.num_classes

		# initialize the matrices

		# Remember! Dicts are not guaranteed to maintain key order but lists are
		# When cycling through the matrix, iterate over the lists, and not the keys of the dict
		self.confusion_matrix = {}
		self.average_class_probability_matrix = {}

		# initialize the rows
		# Does the test set have ground truth?
		if self.test_set.classnames_list:
			for test_set_class_name in self.test_set.classnames_list:
				self.confusion_matrix[ test_set_class_name ] = {}
				self.average_class_probability_matrix[ test_set_class_name ] = {}
		else:
			self.confusion_matrix[ "UNKNOWN" ] = {}
			self.average_class_probability_matrix[ "UNKNOWN" ] = {}
		# now the columns:
		for training_set_class_name in self.training_set.classnames_list:
			for test_set_class_name in self.test_set.classnames_list:
				self.confusion_matrix[ test_set_class_name ][ training_set_class_name ] = 0
				self.average_class_probability_matrix[ test_set_class_name ][ training_set_class_name ] = 0

		self.num_correct_classifications = 0

		from collections import defaultdict
		self.num_classifications_per_class = defaultdict( int )
		self.num_correct_classifications_per_class = defaultdict( int )

		for indiv_result in self.individual_results:
			self.num_classifications_per_class[ indiv_result.ground_truth_class_name ] += 1

			if indiv_result.ground_truth_class_name == indiv_result.predicted_class_name:
				self.num_correct_classifications += 1
				self.num_correct_classifications_per_class[ indiv_result.ground_truth_class_name ] += 1

			test_set_class = indiv_result.ground_truth_class_name
			if test_set_class == None: test_set_class = "UNKNOWN"

			self.confusion_matrix[ test_set_class ][ indiv_result.predicted_class_name ] += 1

			# FIXME: is there any possibility that the order of the values in the marginal
			# probability array don't correspond with the order of the training set classes?
			index = 0
			for training_set_class_name in self.training_set.classnames_list:
				val = indiv_result.marginal_probabilities[ index ]
				self.average_class_probability_matrix[ test_set_class ][ training_set_class_name ] += val
				index += 1

		# Do the averaging for the similarity and avg prob matrices.
		for row in self.test_set.classnames_list:
			for col in self.training_set.classnames_list:
				self.average_class_probability_matrix[ row ][ col ] /= \
				   self.num_classifications_per_class[ row ]

		# The similarity matrix is just the average class probability matrix
		# normalized to have 1's in the diagonal.
		# Doesn't make sense to do this unless the matrix is square
		# if row labels == column labels:
		if self.test_set.classnames_list == self.training_set.classnames_list:
			from copy import deepcopy
			self.similarity_matrix = deepcopy( self.average_class_probability_matrix )
			for row in self.test_set.classnames_list:
				for col in self.training_set.classnames_list:
					self.similarity_matrix[ row ][ col ] /= self.similarity_matrix[ row ][ row ]
					self.similarity_matrix[ col ][ row ] /= self.similarity_matrix[ row ][ row ]


		self.classification_accuracy = float( self.num_correct_classifications) / float( self.num_classifications )

	#==============================================================
	@output_railroad_switch
	def Print( self ):
		"""Prints out statistics from this batch's image classifications, which include 
		classification accuracy, confusion matrix, similarity matrix, and average class 
		probability matrix."""

		if self.classification_accuracy == None:
			self.GenerateStats()

		print "==========================================="
		print "Batch summary:"
		print "Total number of classifications: {0}".format( self.num_classifications )
		print "Total number of CORRECT classifications: {0}".format( self.num_correct_classifications )
		print "Total classification accuracy: {0:0.4f}".format( self.classification_accuracy )
		if not self.figure_of_merit == None:
			print "Standard Error: {0:0.4f}".format( self.figure_of_merit )
			print "Pearson Coefficient: {0:0.4f}".format( self.pearson_coeff )
			print "Spearman Coefficient: {0:0.4f}".format( self.spearman_coeff )

		print ""

		# Remember: iterate over the sorted list of class names, not the keys in the dict,
		# because the dict keys aren't guaranteed to be ordered, nor is test_set.classnames_list.
		# Also remember: there may be different numbers of classes in train and test set
		# or they may be named differently.

		train_set_class_names = sorted( self.training_set.classnames_list )
		test_set_class_names = sorted( self.test_set.classnames_list )

		column_headers = "\t".join( train_set_class_names )
		column_headers += "\n"
		column_headers += "\t".join( [ '-'*len(name) for name in train_set_class_names ] )

		print "Confusion Matrix:"
		print column_headers

		# See how the row labels are test set class names
		# and the column labels are training set class names?

		for row_name in test_set_class_names:
			line = ""
			for col_name in train_set_class_names:
				line += '{0}\t'.format( self.confusion_matrix[ row_name ][ col_name ] )
			print line

		print ""

		if self.similarity_matrix:
			print "Similarity Matrix:"
			print column_headers
			for row_name in test_set_class_names:
				line = ""
				for col_name in train_set_class_names:
					line += '{0:0.4f}\t'.format( self.similarity_matrix[ row_name ][ col_name ] )
				print line
			print ""

		print "Average Class Probability Matrix:"
		print column_headers
		for row_name in test_set_class_names:
			line = ""
			for col_name in train_set_class_names:
				line += '{0:0.4f}\t'.format( self.average_class_probability_matrix[ row_name ][ col_name ] )
			print line
		print ""
	
	#==============================================================
	@output_railroad_switch
	def PrintDistanceMatrix( self, method = 'mps' ):
		"""Generate a distance matrix indicating similarity of image classes to one another.
		The output can be fed into the PHYLIP or other program to create a dendrogram.

		Methods can be {'max', 'mean', 'top', 'bottom', 'mps'}.
		
		To generate a distance matrix using marginal probability space ('mps'), it is legal for
		the test set and training set to have a different number of classes, since the coefficients
		in the average class probability matrix are used as coordinates, but otherwise
		the the input (similarity) matrix must be square."""

		# The source matrix is a numpy matrix of coefficients from which the 
		# distances will be derived
		source_matrix = np.empty( (self.test_set.num_classes, self.test_set.num_classes) ) #initialize
		row_index = 0
		for row in self.test_set.classnames_list:
			col_index = 0
			for col in self.training_set.classnames_list:
				source_matrix[ row_index ][ col_index ] = self.average_class_probability_matrix[row][col]
				col_index += 1
			row_index += 1

		output_matrix = {}

		# initialize the rows
		for test_set_class_name in self.test_set.classnames_list:
			output_matrix[ test_set_class_name ] = {}
		# now the columns:
		for training_set_class_name in self.training_set.classnames_list:
			for test_set_class_name in self.test_set.classnames_list:
				output_matrix[ test_set_class_name ][ training_set_class_name ] = 0

		if method == 'mps':
			if self.average_class_probability_matrix == None:
				self.GenerateStats()

			for row in range( self.test_set.num_classes ):
				for col in range( self.test_set.num_classes ):
					if row == col:
						output_matrix[ row, col ] = 0
					else:
						output_matrix[ row, col ] = np.sum( np.square( source_matrix[row]-source_matrix[col] ) )

		else:
			if self.similarity_matrix == None:
				self.GenerateStats()
			# Have to check again, since for some experiments it doesn't make sense to have a 
			# similarity matrix, e.g., when test set ground truth is not known.
			if self.similarity_matrix == None:
				raise ValueError( 'Similarity matrix is not possible with this data set. ' + \
				                  'Possible reason for this is that your test set does not have ' + \
				                  'ground truth defined.' )
			

			if method == 'max':
				raise NotImplementedError
			elif method == 'mean':
				raise NotImplementedError
			elif method == 'top':
				raise NotImplementedError
			elif method == 'bottom':
				raise NotImplementedError
			else:
				raise ValueError( "Distance matrix method {0} not recognized. Valid options are ".format(\
													method ) + "'max', 'mean', 'top', 'bottom', 'mps'" )

		#column_headers = "\t".join( self.test_set.classnames_list )
		#column_headers += "\n"
		#column_headers += "\t".join( [ '-'*len(name) for name in self.test_set.classnames_list ] )
		#print "Distance Matrix (method = '{0}'):".format( method )
		#print column_headers
		print self.test_set.num_classes
		for row in range( self.test_set.num_classes ):
			line = "{0}\t".format( self.test_set.classnames_list[ row ] )
			for col in range( self.test_set.num_classes ):
				line += '{0:0.4f}\t'.format( output_matrix[ row, col ] )
			print line
		#print ""


	#==============================================================
	@classmethod
	@output_railroad_switch
	def New( cls, training_set, test_set, feature_weights, batch_number = None, batch_name = None, quiet = False, norm_factor_threshold = None):
		"""The equivalent of the "wndcharm classify" command in the command line implementation
		of WND-CHARM. Input a training set, a test set, and feature weights, and returns a
		new instance of a DiscreteBatchClassificationResult, with self.individual_results
		filled with a new instances of DiscreteImageClassificationResult.

		FIXME: What happens when the ground truth is not known? Currently they would all be shoved
					 into class 1, might not be a big deal since class name should be something
					 like "UNKNOWN"
		"""

		# type checking
		if not isinstance( training_set, FeatureSet_Discrete ):
			raise ValueError( 'First argument to New must be of type "FeatureSet_Discrete", you gave a {0}'.format( type( test_set ).__name__ ) )	
		if not isinstance( test_set, FeatureSet_Discrete ):
			raise ValueError( 'Second argument to New must be of type "FeatureSet_Discrete", you gave a {0}'.format( type( test_set ).__name__ ) )	
		if not isinstance( feature_weights, FeatureWeights ):
			raise ValueError( 'Third argument to New must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )
	
		# feature comparison
		if test_set.featurenames_list != feature_weights.names:
			raise ValueError( "Can't classify, features in test set don't match features in weights. Try translating feature names from old style to new, or performing a FeatureReduce()" )
		if test_set.featurenames_list != training_set.featurenames_list:
			raise ValueError( "Can't classify, features in test set don't match features in training set. Try translating feature names from old style to new, or performing a FeatureReduce()" )

		train_set_len = len( training_set.featurenames_list )
		test_set_len = len( test_set.featurenames_list )
		feature_weights_len = len( feature_weights.names )

		print "Classifying test set '{0}' ({1} features) against training set '{2}' ({3} features)".\
					format( test_set.source_path, test_set_len, training_set.source_path, train_set_len )

		if not quiet:
			column_header = "image\tnorm. fact.\t"
			column_header +=\
				"".join( [ "p(" + class_name + ")\t" for class_name in training_set.classnames_list ] )
			column_header += "act. class\tpred. class\tpred. val."
			print column_header

		batch_result = cls( training_set, test_set, feature_weights )
		batch_result.name = batch_name

		train_set_interp_coeffs = None
		if training_set.interpolation_coefficients:
			train_set_interp_coeffs = np.array( training_set.interpolation_coefficients )
			batch_result.predicted_values = []

		test_set_interp_coeffs = None
		if test_set.interpolation_coefficients:
			test_set_interp_coeffs = np.array( test_set.interpolation_coefficients )
			batch_result.ground_truth_values = []

		for test_class_index in range( test_set.num_classes ):
			num_class_imgs, num_class_features = test_set.data_list[ test_class_index ].shape
			for test_image_index in range( num_class_imgs ):
				one_image_features = test_set.data_list[ test_class_index ][ test_image_index,: ]
				result = DiscreteImageClassificationResult._WND5( training_set, one_image_features, feature_weights.values )
				
				if norm_factor_threshold and (result.normalization_factor > norm_factor_threshold):
					continue
				result.source_file = test_set.imagenames_list[ test_class_index ][ test_image_index ]
				result.ground_truth_class_name = test_set.classnames_list[ test_class_index ]
				result.batch_number = batch_number
				result.name = batch_name
				marg_probs = np.array( result.marginal_probabilities )
				result.predicted_class_name = training_set.classnames_list[ marg_probs.argmax() ]

				# interpolated value, if applicable
				if train_set_interp_coeffs is not None:
					interp_val = np.sum( marg_probs * train_set_interp_coeffs )
					result.predicted_value = interp_val
					batch_result.predicted_values.append( interp_val )

				if test_set_interp_coeffs is not None:
					result.ground_truth_value = test_set_interp_coeffs[ test_class_index ]
					batch_result.ground_truth_values.append( test_set_interp_coeffs[ test_class_index ] )

				if not quiet:
					result.Print( line_item = True )
				batch_result.individual_results.append( result )

		return batch_result

#=================================================================================
class ContinuousBatchClassificationResult( BatchClassificationResult ):
	"""Concrete class that is a container for ContinuousImageClassificationResult instances.
	Use this class to run classification experiments and to generate statistics on their results.
	For continuous (Pearson) classifiers, the figure_of_merit is the standard error between
	predicted and ground truth."""

	def __init__( self, training_set, test_set, feature_weights ):
		"""constructor"""
		super( ContinuousBatchClassificationResult, self ).__init__( training_set, test_set, feature_weights )
		self.predicted_values = []

	#=====================================================================
	def GenerateStats( self ):
		"""Calls base class method to run a standard error analysis if ground_truth/
		predicted vals exist."""
		super( ContinuousBatchClassificationResult, self ).GenerateStats()

	#=====================================================================
	@output_railroad_switch
	def Print( self ):
		"""Calculates and outputs batch-level statistics based on the
		ContinuousImageClassificationResults contained in self.individualresults."""
		if self.figure_of_merit == None:
			self.GenerateStats()

		print "==========================================="
		print "Number of observations: {0}".format( self.num_classifications )
		if self.figure_of_merit != None:
			print "Standard error of predicted vs. ground truth values: {0}".format( self.figure_of_merit )
		#print "p-value for this split: {0}".format( self.p_value )
		#print "Standard error for this split: {0}".format( self.std_err )


	#=====================================================================
	@classmethod
	def New( cls, test_set, feature_weights, quiet = False, batch_number = None, batch_name = None):
		"""The equivalent of "wndchrm classify -C" in the command line implenentation of
		WND-CHARM. """

		# type checking
		if not isinstance( test_set, FeatureSet_Continuous ):
			raise ValueError( 'First argument to New must be of type "FeatureSet_Continuous", you gave a {0}'.format( type( test_set ).__name__ ) )	
		if not isinstance( feature_weights, ContinuousFeatureWeights ):
			raise ValueError( 'Second argument to New must be of type "ContinuousFeatureWeights", you gave a {0}'.format( type( feature_weights ).__name__ ) )

		# feature comparison
		if test_set.featurenames_list != feature_weights.names:
			raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )

		# say what we're gonna do
		if not quiet:
			out_str = 'Classifying test set "{0}" ({1} images, {2} features)\n\tagainst training set "{3}" ({4} images)'
			print out_str.format( test_set.source_path, \
			                      test_set.num_images, \
			                      len( test_set.featurenames_list ), \
			                      feature_weights.associated_training_set.source_path, \
			                      feature_weights.associated_training_set.num_images )

		if not quiet:
			column_header = "image\tground truth\tpred. val."
			print column_header

		batch_result = cls( feature_weights.associated_training_set, test_set, feature_weights )
		batch_result.name = batch_name

		if test_set.ground_truths is not None and len( test_set.ground_truths ) != 0:
			batch_result.ground_truth_values = test_set.ground_truths

		for test_image_index in range( test_set.num_images ):
			one_image_features = test_set.data_matrix[ test_image_index,: ]
			result = ContinuousImageClassificationResult._LinearRegression( one_image_features, feature_weights )
			result.batch_number = batch_number
			result.name = batch_name
			result.source_file = test_set.imagenames_list[ test_image_index ]
			result.ground_truth_value = test_set.ground_truths[ test_image_index ]
			batch_result.predicted_values.append( result.predicted_value )

			if not quiet:
				result.Print( line_item = True )
			batch_result.individual_results.append( result )

		batch_result.GenerateStats()
		return batch_result


#============================================================================
class ClassificationExperimentResult( BatchClassificationResult ):
	"""Abstract base class which serves as a container for BatchClassificationResult instances
	and their associated statistics.
	
	N.B. This is a container for batches, not for individual results, i.e.,
	the list self.individual_results contains instances of type BatchClassificationResult."""

	#: A dictionary where the name is the key, and the value is a list of individual results
	accumulated_individual_results = None

	feature_weight_statistics = None

	#: keep track of stats related to predicted values for reporting purposes
	individual_stats = None
	#=====================================================================
	def GenerateStats( self ):
		"""Not yet implemented.
		
		Considerations for future implementation.
		1. The test set may or may not have ground truth (discrete and continuous)
		2. The results may not have a predicted value (discrete only)
		3. Continuous classifications do not have marginal probabilities
		4. Hybrid test sets (discrete test sets loaded into a continuous test set)
		   have "pseudo-classes," i.e., easily binnable ground truth values."""

		# Step 1: feature weight statistics:

		feature_weight_lists = {}
		for batch_result in self.individual_results:
			weight_names_and_values = zip( batch_result.feature_weights.names, batch_result.feature_weights.values)
			for name, weight in weight_names_and_values:
				if not name in feature_weight_lists:
					feature_weight_lists[ name ] = []
				feature_weight_lists[ name ].append( weight )

		for feature_name in feature_weight_lists:
			feature_weight_lists[ feature_name ] = np.array( feature_weight_lists[ feature_name ] )

		fwl = feature_weight_lists
		feature_weight_stats = [ (np.mean(fwl[fname]), len(fwl[fname]), np.std(fwl[fname]), np.min(fwl[fname]), np.max(fwl[fname]), fname) for fname in fwl ]
		# Sort on mean values, i.e. index 0 of tuple created above
		sort_func = lambda A, B: cmp( A[0], B[0] )
		self.feature_weight_statistics = sorted( feature_weight_stats, sort_func, reverse = True )

		#Step 2: aggregate all results across splits for individual images
		# FIXME: Do it!


	#=====================================================================
	@output_railroad_switch
	def PredictedValueAnalysis( self ):
		"""For use in only those classification experiments where predicted values are generated
		as part of the classification.
		
		This function is meant to elucidate the amount of variability of classifications 
		across batches. ImageClassificationResult imformation is aggregated for each individual
		image/ROI encountered, and statistics are generated for each image/ROI and printed out."""

		if self.individual_results == 0:
			raise ValueError( 'No batch results to analyze' )

		self.predicted_values = []
		self.ground_truth_values = []

		self.accumulated_individual_results = {}
		self.individual_stats = {}

		for batch in self.individual_results:
			for result in batch.individual_results:
				if not result.source_file in self.accumulated_individual_results:
					# initialize list of individual results for this file
					self.accumulated_individual_results[ result.source_file ] = []
				self.accumulated_individual_results[ result.source_file ].append( result )
	
		for filename in self.accumulated_individual_results:
			vals = np.array( [result.predicted_value for result in self.accumulated_individual_results[filename] ])
			self.ground_truth_values.append( self.accumulated_individual_results[filename][0].ground_truth_value )
			self.predicted_values.append( np.mean(vals) )
			self.individual_stats[filename] = ( len(vals), np.min(vals), np.mean(vals), \
			                                    np.max(vals), np.std(vals) ) 

		print "\n\nIndividual results\n========================\n"
		mp = "  "
		discrlineoutstr = "\tsplit {split_num:02d} '{batch_name}': predicted: {pred_class}, actual: {actual_class}. Norm dists: ( {norm_dists} ) Interp val: {pred_val:0.3f}"
		contlineoutstr = "\tsplit {split_num:02d} '{batch_name}': actual: {actual_class}. Predicted val: {pred_val:0.3f}"
		outstr = "\t---> Tested {0} times, low {1:0.3f}, mean {2:0.3f}, high {3:0.3f}, std dev {4:0.3f}"

		#create view
		res_dict = self.accumulated_individual_results

		# sort by ground truth, then alphanum
		sort_func = lambda A, B: cmp( A, B ) if res_dict[A][0].ground_truth_value == res_dict[B][0].ground_truth_value else cmp( res_dict[A][0].ground_truth_value, res_dict[B][0].ground_truth_value  ) 
		sorted_images = sorted( self.accumulated_individual_results.iterkeys(), sort_func )

		for samplename in sorted_images:
			print 'File "' + samplename + '"'
			for result in self.accumulated_individual_results[ samplename ]:

				if isinstance( result, DiscreteImageClassificationResult ):
					marg_probs = [ "{0:0.3f}".format( num ) for num in result.marginal_probabilities ]
					print discrlineoutstr.format( split_num = result.batch_number, \
				                         batch_name = result.name, \
				                         pred_class = result.predicted_class_name, \
				                         actual_class = result.ground_truth_value, \
				                         norm_dists = mp.join( marg_probs ), \
				                         pred_val = result.predicted_value )
				elif isinstance( result, ContinuousImageClassificationResult ):
					print contlineoutstr.format( split_num = result.batch_number, \
				                         batch_name = result.name, \
				                         actual_class = result.ground_truth_value, \
				                         pred_val = result.predicted_value )
				else:
					raise ValueError( 'expected an ImageClassification result but got a {0} class'.\
				  	       format( type( result ).__name__ ) ) 
			print outstr.format( *self.individual_stats[ samplename ] )

	#=====================================================================
	@output_railroad_switch
	def WeightsOptimizationAnalysis( self ):
		"""For use in only those classification experiments where predicted values are generated
		as part of the classification.
		
		This function is used when the user has varied the number/composition of image features
		in each of the batches and wants to figure out which classifier results in the best
		figure of merit. This function evaluates which of the batches was created using
		the optimal classifier configuration and outputs that."""

		print "Calculating optimized continuous classifier for training set {0}".\
				format( self.training_set.source_path )

		best_classification_result = None
		best_weights = None
		opt_number_features = None
		min_std_err = float( "inf" )

		for batch_result in self.individual_results:
			weights = batch_result.feature_weights
			num_features = len( weights.values )
			if batch_result.figure_of_merit < min_std_err:
				min_std_err = batch_result.figure_of_merit
				best_classification_result = batch_result
				best_weights = weights
				opt_number_features = num_features

		print "Optimum number of features: {0}".format( opt_number_features )
		best_classification_result.Print()
		best_weights.Print()
		print ""
		print "Aggregate feature weight analysis:"
		print "-----------------------------------"
		print "Legend:"
		print "NUM - Number of features used in aggregate / Individual feature rank"
		print "ASE - Standard Error of Final Predicted Value (using aggregated feature) vs ground truth"
		print "APC - Pearson correlation coefficient of Final Predicted Values vs ground truth"
		print "APE - Standard Error of APC"
		print "APP - P-value of APC"
		print "ASC - Spearman correlation coefficient of Final Predicted Values vs ground truth"
		print "APP - P-value of ASC"
		print ""

		print "NUM\tASE\tAPC\tAPE\tAPP\tASC\tAPP"
		print "===\t===\t===\t===\t===\t===\t==="
		for result in self.individual_results:
			line_item = "{0}\t".format( len( result.feature_weights.values ) ) # NUM
			line_item += "{0:.4f}\t".format( result.figure_of_merit ) # ASE
			line_item += "{0:.4f}\t".format( result.pearson_coeff ) # APC
			line_item += "{0:.4f}\t".format( result.pearson_std_err ) # APE
			line_item += "{0:.4f}\t".format( result.pearson_p_value ) # APP
			line_item += "{0:.4f}\t".format( result.spearman_coeff ) # ASC
			line_item += "{0:.4f}\t".format( result.spearman_p_value ) # ASP
			print line_item


#============================================================================
class DiscreteClassificationExperimentResult( ClassificationExperimentResult ):
	"""Concrete class which serves as a container for DiscreteBatchClassificationResults
	and their associated statistics. The information contained here comprises everything
	that would appear in a HTML file generated by the C++ implementation of WND-CHARM.

	In this subclass, the figure of merit is self.classification_accuracy"""

	num_correct_classifications = None

	confusion_matrix = None
	average_similarity_matrix = None
	average_class_probability_matrix = None

	#=====================================================================
	def GenerateStats( self ):
		"""Not fully implemented yet. Need to implement generation of confusion, similarity, and
		average class probability matrices from constituent batches."""

		# Base class does feature weight analysis
		super( DiscreteClassificationExperimentResult, self ).GenerateStats()
		
		self.num_correct_classifications = 0
		for batch_result in self.individual_results:
			for indiv_result in batch_result.individual_results:
				self.num_classifications += 1
				if indiv_result.ground_truth_class_name == indiv_result.predicted_class_name:
					self.num_correct_classifications += 1

		self.figure_of_merit = float( self.num_correct_classifications) / float( self.num_classifications )

		#FIXME: Create confusion, similarity, and class probability matrices


	#=====================================================================
	@output_railroad_switch
	def Print( self ):
		"""Generate and output statistics across all batches, as well as the figures of merit
		for each individual batch."""
		
		if self.figure_of_merit == None:
			self.GenerateStats()

		print "==========================================="
		print "Experiment name: {0}".format( self.name )
		print "Summary:"
		print "Number of batches: {0}".format( len( self.individual_results ) )
		print "Total number of classifications: {0}".format( self.num_classifications )
		print "Total number of CORRECT classifications: {0}".format( self.num_correct_classifications )
		print "Total classification accuracy: {0:0.4f}\n".format( self.figure_of_merit )

		print "Batch Accuracies:"
		print "#\tname\tclassification accuracy"
		print "------------------------------------"

		count = 1
		for batch_result in self.individual_results:
			if batch_result.figure_of_merit == None:
				batch_result.GenerateStats()
			print "{0}\t\"{1}\"\t{2:0.4f}".format( count, batch_result.name, batch_result.classification_accuracy )
			count += 1

		outstr = "{0}\t{1:0.3f}\t{2}\t{3:0.3f}\t{4:0.3f}\t{5:0.3f}\t{6}"
		print "Feature Weight Analysis:"
		print "Rank\tmean\tcount\tStdDev\tMin\tMax\tName"
		print "----\t----\t-----\t------\t---\t---\t----"
		count = 1
		for fw_stat in self.feature_weight_statistics:
			print outstr.format( count, *fw_stat )
			count += 1



#============================================================================
class ContinuousClassificationExperimentResult( ClassificationExperimentResult ):
	"""Concrete class which serves as a container for ContinuousBatchClassificationResult instances
	and their associated statistics. The information contained here comprises everything
	that would appear in a HTML file generated by the C++ implementation of WND-CHARM.

	In this subclass, the figure of merit is the average standard error arcoss batches."""

	def __init__( self, name = None ):
		self.name = name
		super( ContinuousClassificationExperimentResult, self ).__init__()

	#=====================================================================
	def GenerateStats( self ):
		"""Calculates statistics describing how well predicted values
		correlate with ground truth across all batches.

		Requires scipy.stats package to be installed"""

		from itertools import chain

		lists_of_ground_truths = []
		lists_of_predicted_values = []

		self.figure_of_merit = 0
		for batch_result in self.individual_results:
			self.num_classifications += len( batch_result.individual_results )
			if batch_result.figure_of_merit == None:
				batch_result.GenerateStats()
			lists_of_ground_truths.append( batch_result.ground_truth_values )
			lists_of_predicted_values.append( batch_result.predicted_values )

		self.ground_truth_values = list( chain( *lists_of_ground_truths ) )
		self.predicted_values = list( chain( *lists_of_predicted_values ) )
	
		gt = np.array( self.ground_truth_values )
		pv = np.array( self.predicted_values )
		
		diffs = gt - pv
		diffs = np.square( diffs )
		err_sum = np.sum( diffs )

		import math; from scipy import stats
		self.figure_of_merit = math.sqrt( err_sum / self.num_classifications )

		# For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
		np.seterr (under='ignore')
		slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
								 stats.linregress( self.ground_truth_values, self.predicted_values )

		self.spearman_coeff, self.spearman_p_value =\
					 stats.spearmanr( self.ground_truth_values, self.predicted_values )
		np.seterr (all='raise')


	#=====================================================================
	@output_railroad_switch
	def Print( self ):
		"""Output statistics from this experiment."""
		if self.figure_of_merit == None:
			self.GenerateStats()

		print "==========================================="
		print "Experiment name: {0}".format( self.name )
		print "Summary:"
		print "Number of batches: {0}".format( len( self.individual_results ) )
		print "Total number of classifications: {0}".format( self.num_classifications )
		print "Total standard error: {0}".format( self.figure_of_merit )
		print "Pearson corellation coefficient: {0}".format( self.pearson_coeff )
		print "Pearson p-value: {0}".format( self.pearson_p_value )		


#============================================================================
class BaseGraph( object ):
	"""An abstract base class that is supposed to hold onto objects on which to call
	matplotlib.pyplot API methods."""

	# general stuff:
	chart_title = None
	file_name = None
	batch_result = None

	# pyplot-specific stuff
	figure = None
	main_axes = None

	def SaveToFile( self, filepath ):
	
		if self.figure == None:
			raise ValueError( 'No figure to save!' )
		self.figure.savefig( filepath )
		print 'Wrote chart "{0}" to file "{1}.png"'.format( self.chart_title, filepath )
			
#============================================================================
class PredictedValuesGraph( BaseGraph ):
	"""This is a concrete class that can produce two types of graphs that are produced
	from ImageClassificationResult data stored in a BatchClassificationResult."""
	
	# members concerned with class-dependent figures 
	grouped_coords = None
	num_classes = None
	classnames_list = None
	class_values = None

	#=================================================================
	def __init__( self, result ):
		"""Constructor sorts ground truth values contained in BatchClassificationResult
		and loads them into self.grouped_coords"""

		self.batch_result = result
		self.grouped_coords = {}

		self.classnames_list = result.test_set.classnames_list
		self.class_values = result.test_set.interpolation_coefficients

		#FIXME: implement user-definable bin edges

		result.RankOrderSort()
		whole_list = zip( result.ground_truth_values, result.predicted_values )

		class_index = 0
		for class_name in self.classnames_list:
			self.grouped_coords[ class_name ] = []
			for coords in whole_list:
				if coords[0] == self.class_values[ class_index ]:
					self.grouped_coords[ class_name ].append( coords )
			class_index += 1

	#=====================================================================
	def SaveToFile( self, filepath ):
		"""Calls base class function"""
		super( PredictedValuesGraph, self ).SaveToFile( filepath )

	#=====================================================================
	def RankOrderedPredictedValuesGraph( self, chart_title = None ):
		"""This graph visualizes the distribution of predicted values generated by classification.
		For each individual ImageClassificationResult with ground truth value (i.e., class id) and
		predicted value, all results are grouped within their class, sorted by predicted value
		in ascending order, then ploted side-by-side.
		
		Required the package matplotlib to be installed."""

		from matplotlib import pyplot

		color = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

		self.figure = pyplot.figure()
		self.main_axes = self.figure.add_subplot(111)
		if chart_title:
			self.chart_title = chart_title
			self.main_axes.set_title( chart_title )
		self.main_axes.set_title( chart_title )
		self.main_axes.set_xlabel( 'count' )
		self.main_axes.set_ylabel( 'Predicted Value Scores' )

		abscissa_index = 1

		for class_name in self.classnames_list:
			ground_truth_vals, predicted_vals = zip( *self.grouped_coords[ class_name ] )
			x_vals = [ i + abscissa_index for i in range( len( ground_truth_vals ) ) ]
			self.main_axes.scatter( x_vals, predicted_vals, color = next( color ), marker = 'o', label = class_name )
			abscissa_index += len( ground_truth_vals )

		# FIXME: put legend in bottom left
		self.main_axes.legend( loc = 'lower right')
		
	#=====================================================================
	def KernelSmoothedDensityGraph( self, chart_title = None ):
		"""This graph visualizes the distribution of predicted values generated by classification.
		A kernel-smoothed probability density function is plotted for each image class on
		the same chart allowing comparison of distribution of predicted values amoung image class.
		
		Requires the packages matplotlib and scipy. Uses scipy.stats.gaussian_kde to
		generate kernel-smoothed probability density functions."""

		from matplotlib import pyplot

		color = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

		self.figure = pyplot.figure()
		self.main_axes = self.figure.add_subplot(111)
		if chart_title:
			self.chart_title = chart_title
			self.main_axes.set_title( chart_title )
		self.main_axes.set_xlabel( 'Age score' )
		self.main_axes.set_ylabel( 'Probability density' )

		from scipy import stats

		for class_name in self.classnames_list:
			ground_truth_vals, predicted_vals = zip( *self.grouped_coords[ class_name ] )

			pred_vals = np.array( predicted_vals )
			lobound = pred_vals.min()
			hibound = pred_vals.max()
			kernel_smoother = stats.gaussian_kde( pred_vals )
			intervals = np.mgrid[ lobound:hibound:100j ]
			density_estimates = kernel_smoother.evaluate( intervals )
			self.main_axes.plot( intervals, density_estimates, color = next( color ), linewidth = 3, label = class_name )

		self.main_axes.legend()

#============================================================================
class FeatureTimingVersusAccuracyGraph( BaseGraph ):
	"""A cost/benefit analysis of the number of features used and the time it takes to calculate
	that number of features for a single image"""

	#FIXME: Add ability to do the first 50 or 100 features, make the graph, then
	#       ability to resume from where it left off to do the next 50.

	timing_axes = None

	def __init__( self, training_set, feature_weights, test_image_path, chart_title = None, max_num_features = 300):
		import time
		timings = []
	
		experiment = DiscreteClassificationExperimentResult( training_set, training_set, feature_weights)
		for number_of_features_to_use in range( 1, max_num_features + 1 ):

			reduced_ts = None
			reduced_fw = None
			three_timings = []
			# Take the best of 3
			for timing in range( 3 ):
				# Time the creation and classification of a single signature
				t1 = time.time()
				reduced_fw = feature_weights.Threshold( number_of_features_to_use )
				sig = Signatures.NewFromFeatureNameList( test_image_path, reduced_fw.names )
				reduced_ts = training_set.FeatureReduce( reduced_fw.names )
				sig.Normalize( reduced_ts )
		
				result = DiscreteImageClassificationResult.NewWND5( reduced_ts, reduced_fw, sig )
				result.Print()
				# FIXME: save intermediates just in case of interruption or parallization
				# result.PickleMe()
				t2 = time.time()
				three_timings.append( t2 - t1 )

			timings.append( min( three_timings ) )

			# now, do a fit-on-fit test to measure classification accuracy
			batch_result = DiscreteBatchClassificationResult.New( reduced_ts, reduced_ts, reduced_fw )
			batch_result.Print()
			experiment.individual_results.append( batch_result )

		from matplotlib import pyplot

		x_vals = list( range( 1, max_num_features + 1 ) )

		self.figure = pyplot.figure()
		self.main_axes = self.figure.add_subplot(111)
		if chart_title == None:
			self.chart_title = "Feature timing v. classification accuracy"	
		else:
			self.chart_title = chart_title
		self.main_axes.set_title( self.chart_title )
		self.main_axes.set_xlabel( 'Number of features' )
		self.main_axes.set_ylabel( 'Classification accuracy (%)', color='b' )
		classification_accuracies = \
		  [ batch_result.classification_accuracy * 100 for batch_result in experiment.individual_results ]

		self.main_axes.plot( x_vals, classification_accuracies, color='b', linewidth=2 )
		for tl in self.main_axes.get_yticklabels():
			tl.set_color('b')	

		self.timing_axes = self.main_axes.twinx()
		self.timing_axes.set_ylabel( 'Time to calculate features (s)', color='r' )
		self.timing_axes.plot( x_vals, timings, color='r' )
		for tl in self.timing_axes.get_yticklabels():
			tl.set_color('r')	

	def SaveToFile( self, filepath ):
		super( FeatureTimingVersusAccuracyGraph, self ).SaveToFile( filepath )

#============================================================================
class Dendrogram( object ):
	"""Not implemented. In the future might use scipy.cluster (no unrooted dendrograms though!)
	or Biopython.Phylo to visualize. Perhaps could continue C++ implementation's use of PHYLIP's
	Fitch-Margoliash program "fitch" to generate Newick phylogeny, and visualize using
	native python tools."""
	pass

#================================================================

initialize_module()

