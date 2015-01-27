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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

# wndcharm.py has the definitions of all the SWIG-wrapped primitive C++ WND_CHARM objects.
import wndcharm
import numpy as np

# ============================================================
# BEGIN: Initialize module level globals
Algorithms = []
Transforms = []
# The numbers *must* be consistent with what's defined in wndchrm C-codebase.
# Feature Vector revved to v3.x on 20141201 due to issue 39
feature_vector_major_version = 3
# Feature vector lengths in current version
# #define NUM_LC_FEATURES  4059
# #define NUM_L_FEATURES   2919
# #define NUM_C_FEATURES   2199
# #define NUM_DEF_FEATURES 1059
# These are definitions for Version 2 and 3 features.
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

def initialize_module(): 
	"""If you're going to calculate any features, you need this stuff.
	FIXME: Rig this stuff to load only on demand."""
	
	global Algorithms
	global Transforms
	# The verbosity is set by the environment variable WNDCHRM_VERBOSITY, in wndchrm_error.cpp
	# wndcharm.cvar.verbosity = 7

	# These are the auto-registered ComputationTasks, which come in different flavors
	all_tasks = wndcharm.ComputationTaskInstances.getInstances()
	for task in all_tasks:
		if task.type == task.ImageTransformTask:
			Transforms.append (task)
#			print task.name + " added to Transforms"
		elif task.type == task.FeatureAlgorithmTask:
			Algorithms.append (task)
#			print task.name + " added to Algorithms"

	# The standard feature plans get loaded as needed from C++ statics, so they don't need to be initialized.
	# Standard sets (other plan "parts" are in Tasks.h under StdFeatureComputationPlans)
	# 	getFeatureSet();
	# 	getFeatureSetColor();
	# 	getFeatureSetLong();
	# 	getFeatureSetLongColor();

	# e. g.:
# 	small_feature_plan = wndcharm.StdFeatureComputationPlans.getFeatureSet()
# 	print "small feature set groups:"
# 	last_feature_group = None;
# 	for i in range( 0, small_feature_plan.n_features):
# 		feature_group = small_feature_plan.getFeatureGroupByIndex(i)
# 		feature_name = small_feature_plan.getFeatureNameByIndex(i)
# 		if feature_group.name != last_feature_group:
# 			print "feature_group "+feature_group.name
# 		last_feature_group = feature_group.name
# 		print "  feature_name "+feature_name

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
		elif "output_stream" in kwargs:
			output_stream = kwargs[ "output_stream" ]
			del kwargs[ "output_stream" ]
			print 'Saving output of function "{0}()" to stream'.format(\
			      method_that_prints_output.__name__)
			import sys
			backup = sys.stdout
			try:
				sys.stdout = output_stream
				retval = method_that_prints_output( *args, **kwargs )
			finally:
				sys.stdout = backup
		else:
			retval = method_that_prints_output( *args, **kwargs )
		return retval

	return print_method_wrapper

def normalize_by_columns( full_stack, mins=None, maxs=None ):
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

# Turn off numpy warnings, since we're taking care of invalid values explicitly
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
	# Note the deep copy to change the numpy parameter in-place.
	full_stack[:] = full_stack_m.filled (0) * 100.0

	# return settings to original
	np.seterr(**oldsettings)

	return (mins,maxs)


def CheckIfClassNamesAreInterpolatable( classnames_list ):
	"""N.B., this method takes only the first number it finds in the class label."""

	import re
	p = re.compile( r'(-?\d*\.?\d+)' )
	interp_coeffs = []
	for class_name in classnames_list:
		m = p.search( class_name )
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
	Each call to sample returns the next wndcharm.ImageMatrix in the sample set.
	The constructor has three required parameters.
	The image parameter can be a path to an image file or a wndcharm.ImageMatrix
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
			self.image = wndcharm.ImageMatrix()
			if 1 != self.image.OpenImage( image_in, 0, None, 0, 0 ):
				raise ValueError( 'Could not build an ImageMatrix from {0}, check the file.'.format( image_in ) )
		elif isinstance (image_in, wndcharm.ImageMatrix):
			self.image = image_in
		else:
			raise ValueError("image parameter 'image_in' is not a string or a wndcharm.ImageMatrix")

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
				yield wndcharm.ImageMatrix (original, current_x, current_y, current_x+width-1, current_y+height-1,0,0)
				current_x = current_x + width
				self.current_x = current_x
			current_y = current_y + height
			self.current_y = current_y


#############################################################################
# class definition of FeatureWeights
#############################################################################
class FeatureWeights( object ):
	"""FeatureWeights is an abstract base class
	that comprises one-half of a WND-CHARM classifier (the other being a FeatureSpace).
	
	It is a list of strings which are the names of
	individual image descriptors (features) and a corresponding list of doubles which
	are the weights assigned to those features. Since WND-CHARM is a generalized
	pattern recognition algorithm that calculates the same image descriptors for all
	images, it is through features weights that trained classifiers can zero-in on only
	those features which provide distinctiveness across classes and ignore noisy features.
	Thus any instance of a FeatureWeights class is context-specific."""

	def __init__( self, name=None, size=None ):
		self.name = name
		self.associated_training_set = None
		if size is not None:
			self.featurenames_list = [None] * size
			self.values = [None] * size
		else:
			self.featurenames_list = None
			self.values = None

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
	def NewFromFeatureSpace( cls, num_features_to_be_used  ):
		"""@breif Calculate FeatureWeights from a FeatureSpace"""
		raise NotImplementedError

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

	def __init__( self, name=None, size=None ):
		"""Simply calls parent constructor"""
		super( FisherFeatureWeights, self ).__init__( name=name, size=size )

	#================================================================
	@classmethod
	def NewFromFile( cls, weights_filepath ):
		"""Written to read in files created by wndchrm -vw/path/to/weightsfile.txt"""

		weights = cls( name=weights_filepath )
		with open( weights_filepath ) as weights_file:
			# split line "number <space> name"
			raw_vals, raw_names = \
			  zip( *[ line.strip().split( None, 1 ) for line in weights_file.read().splitlines() ] )

		weights.values = [ float( val ) for val in raw_vals ]
		weights.featurenames_list = [None] * len( raw_names )

		for i, name_raw_str in enumerate( raw_names ):
			# getFeatureInfoByName does some checking, returns a None if it can't parse it
			retval = wndcharm.FeatureNames.getFeatureInfoByName( name_raw_str )
			if retval:
				weights.featurenames_list[i] = retval.name
			else:
				weights.featurenames_list[i] = name_raw_str

		return weights

	#================================================================
	@classmethod
	def NewFromFeatureSpace( cls, training_set ):
		"""Takes a FeatureSpace as input and calculates a Fisher score for 
		each feature. Returns a newly instantiated instance of FisherFeatureWeights.

		For:
		N = number of classes
		F = number of features
		It = total number of images in training set
		Ic = number of images in a given class
		"""

		# we deal with NANs/INFs separately, so turn off numpy warnings about invalid floats.
		oldsettings = np.seterr(all='ignore')

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
		intra_class_means = np.empty(
			( training_set.num_classes, len( training_set.featurenames_list ) ) )
		# 2D matrix shape N * F
		intra_class_variances = np.empty(
			( training_set.num_classes, len( training_set.featurenames_list ) ) )
		
		class_index = 0
		for class_feature_matrix in training_set.data_list:
			intra_class_means[ class_index ] = np.mean( class_feature_matrix, axis=0 )
			# Note that by default, numpy divides by N instead of the more common N-1, hence ddof=1.
			intra_class_variances[ class_index ] = np.var( class_feature_matrix, axis=0, ddof=1 )
			class_index += 1

		# 1D matrix 1 * F
		# we deal with NANs/INFs separately, so turn off numpy warnings about invalid floats.
		# for the record, in numpy:
		#     1./0. = inf, 0./inf = 0., 1./inf = 0. inf/0. = inf, inf/inf = nan
		#     0./0. = nan,  nan/0. = nan, 0/nan = nan, nan/nan = nan, nan/inf = nan, inf/nan = nan
		# We can't deal with NANs only, must also deal with pos/neg infs
		# The masked array allows for dealing with "invalid" floats, which includes nan and +/-inf
		denom = np.mean( intra_class_variances, axis=0 )
		denom[denom == 0] = np.nan
		feature_weights_m =  np.ma.masked_invalid (
			( np.square( population_means - intra_class_means ).sum( axis = 0 ) /
		    (training_set.num_classes - 1) ) / denom
		    )
		# return numpy error settings to original
		np.seterr(**oldsettings)

		new_fw = cls()
		new_fw.featurenames_list = training_set.featurenames_list[:]
		# the filled(0) method of the masked array sets all nan and infs to 0
		new_fw.values = feature_weights_m.filled(0).tolist()
		new_fw.associated_training_set = training_set
	
		return new_fw

	#================================================================
	def EliminateZeros( self ):
		"""Eliminates any features with a weight of zero, and returns a new instance of
		FisherFeatureWeights without those features."""

		new_weights = FisherFeatureWeights()
		scores = zip( self.featurenames_list, self.values )
		nonzero_scores = [ (name, weight) for name, weight in scores if weight != 0 ]
		new_weights.featurenames_list, new_weights.values = zip( *nonzero_scores )
		return new_weights

	#================================================================
	def Threshold( self, num_features_to_be_used=None, _all=False ):
		"""Returns an instance of a FisherFeatureWeights class with the top n relevant features 
		in order.

		if _all == True: simple reorder by feature rank returning even 0-weighted features
		If _all == 'nonzero': returns non-zero weighted features ranked by weight."""

		if _all:
			num_features_to_be_used = len( self.values )
		# Default is top 15% of features
		elif num_features_to_be_used is None:
			num_features_to_be_used = int( len( self.values ) * 0.15 )
		elif num_features_to_be_used > len( self.values ):
			raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
			                      format( len( self.values ), num_features_to_be_used ) )

		new_weights = self.__class__()
		raw_featureweights = zip( self.values, self.featurenames_list )
		# raw_featureweights is now a list of tuples := [ (value1, name1), (value2, name2), ... ] 

		sorted_featureweights = sorted( raw_featureweights, key=lambda a: a[0], reverse = True )
		# take top N features
		from itertools import islice
		use_these_feature_weights = list( islice( sorted_featureweights, num_features_to_be_used ) )

		if _all is not True:
			# You have a problem if any of the features have corellation coefficients of 0
			for i, feat_fig_of_merit in enumerate( [ line[0] for line in use_these_feature_weights ] ):
				if feat_fig_of_merit == 0:
					if _all == 'non-zero':
						print "Features rank index {0} and below have a correllation coefficient of 0. ".format( i )
						print 'Using {0} features'.format( i )
						use_these_feature_weights = use_these_feature_weights[ : i ]
						break
					else:
						err_msg = "Can't reduce feature weights \"{0}\" to {1} features. ".format( self.name, num_features_to_be_used )
						err_msg += "Features ranked {0} and below have a Fisher score of 0. ".format( i )
						err_msg += "Request less features. "
						raise ValueError( err_msg )

		# we want lists, not tuples!
		new_weights.values, new_weights.featurenames_list =\
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
		raw_featureweights = zip( self.featurenames_list, self.values )

		use_these_feature_weights = \
				list( itertools.islice( raw_featureweights, min_index, max_index ) )
		
		# we want lists, not tuples!
		new_weights.featurenames_list, new_weights.values =\
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
		for i, (val, name) in enumerate( zip( self.values, self.featurenames_list ), start=1 ):
			print "{0}\t{1:.6f}\t{2}".format( i, val, name )
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
	
	def __init__( self, name=None, size=None ):
		"""Constructor"""
		super( ContinuousFeatureWeights, self ).__init__( name=name, size=size )
		if size is not None:
			self.slopes = [None] * size
			self.intercepts = [None] * size
			self.pearson_coeffs = [None] * size
			self.pearson_stderrs = [None] * size
			self.pearson_p_values = [None] * size
			self.spearman_coeffs = [None] * size
			self.spearman_p_values = [None] * size
		else:
			self.slopes = None
			self.intercepts = None
			self.pearson_coeffs = None
			self.pearson_stderrs = None
			self.pearson_p_values = None
			self.spearman_coeffs = None
			self.spearman_p_values = None

	#================================================================
	@classmethod
	def NewFromFeatureSpace( cls, training_set ):
		"""Calculate regression parameters and correlation statistics that fully define
		a continuous classifier.

		At present the feature weights are proportional the Pearson correlation coefficient
		for each given feature."""
		
		from scipy.stats import linregress, spearmanr

		# Known issue: running stats.linregress() with np.seterr (all='raise') has caused
		# arithmetic underflow (FloatingPointError: 'underflow encountered in stdtr' )
		# I think this is something we can safely ignore in this function, and return settings
		# back to normal at the end. -CEC
		np.seterr (under='ignore')

		if training_set.name:
			name = cls.__name__ + ' from training set "' + training_set.name + '"'
		else:
			name = None

		new_fw = cls( name=name, size=training_set.num_features )
		new_fw.associated_training_set = training_set

		#r_val_sum = 0
		r_val_squared_sum = 0
		#r_val_cubed_sum = 0

		ground_truths = np.array( [ float(val) for val in training_set.ground_truths ] )

		new_fw.featurenames_list = training_set.featurenames_list[:]

		for feature_index in range( training_set.num_features ):

			feature_values = training_set.data_matrix[ :, feature_index ]

			slope, intercept, pearson_coeff, p_value, std_err = \
			             linregress( ground_truths, feature_values )

			new_fw.pearson_coeffs[ feature_index ] = pearson_coeff
			new_fw.slopes[ feature_index ] = slope
			new_fw.intercepts[ feature_index ] = intercept
			new_fw.pearson_stderrs[ feature_index ] = std_err
			new_fw.pearson_p_values[ feature_index ] = p_value

			#from math import fabs
			#r_val_sum += fabs(pearson_coeff)
			r_val_squared_sum += pearson_coeff * pearson_coeff
			#r_val_cubed_sum += pearson_coeff * pearson_coeff * pearson_coeff

			try:
				spearman_coeff, spearman_p_val = spearmanr( ground_truths, feature_values )
			except FloatingPointError:
				# to avoid: "FloatingPointError: invalid value encountered in true_divide"
				spearman_coeff, spearman_p_val = ( 0, 1 )

			new_fw.spearman_coeffs[ feature_index ] = spearman_coeff
			new_fw.spearman_p_values[ feature_index ] = spearman_p_val

		#new_fw.values = [ fabs(val) / r_val_sum for val in new_fw.pearson_coeffs ]
		new_fw.values = [val*val / r_val_squared_sum for val in new_fw.pearson_coeffs ]
		#new_fw.values = [val*val*val / r_val_cubed_sum for val in new_fw.pearson_coeffs ]
		
		# Reset numpy 
		np.seterr (all='raise')

		return new_fw

	#================================================================
	def Threshold( self, num_features_to_be_used=None, _all=False, use_spearman=False,
                 min_corr_coeff=None ):
		"""Returns a new instance of a ContinuousFeatureWeights class derived from this
		instance where the number of features has been reduced to only the top n features,
		where n is specified by the num_features_to_be_used argument.

		If min_corr_coeff is specified, the argument num_features_to_be_used is ignored.

		if _all == True implies simple reorder by feature rank, including 0-weighted features
		If _all == 'nonzero', returns non-zero weighted features ranked by weight."""

		if _all:
			num_features_to_be_used = len( self.values )
		elif num_features_to_be_used is None:
			# Default is top 15% of features
			num_features_to_be_used = int( len( self.values ) * 0.15 )
		elif num_features_to_be_used < 1 or num_features_to_be_used > len( self.values ):
			raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
			                      format( len( self.values ), num_features_to_be_used ) ) 

		new_weights = self.__class__()
		if self.name:
			if num_features_to_be_used == len( self.featurenames_list ):
				new_weights.name = self.name + " (rank-ordered)"
			else:
				new_weights.name = self.name + " (top {0} features)".format( num_features_to_be_used )

		if use_spearman:
			abs_corr_coeffs = [ abs( val ) for val in self.spearman_coeffs ]
		else:
			abs_corr_coeffs = [ abs( val ) for val in self.pearson_coeffs ]

		raw_featureweights = zip( abs_corr_coeffs, self.featurenames_list, self.pearson_coeffs, \
		    self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
		    self.spearman_coeffs, self.spearman_p_values )
		
		# sort from max to min
		sorted_featureweights = sorted( raw_featureweights, key=lambda r: r[0], reverse = True )
		
		# take most correlated features, both positive and negative
		if min_corr_coeff is not None:
			try:
				val = abs( float( min_corr_coeff ) )
			except:
				raise ValueError( 'Cannot convert {0} to a float.'.format( min_corr_coeff ) )
			if val <= 0 or val > 1:
				raise ValueError( 'Abs val of min correlation coefficient must be between 0 and 1.' )

			from itertools import takewhile
			use_these_feature_weights = list( takewhile( lambda x: x[0]>min_corr_coeff, sorted_featureweights ) )
		else:
			from itertools import islice
			use_these_feature_weights = list( islice( sorted_featureweights, num_features_to_be_used ) )
			if _all is not True:
				# You have a problem if any of the features have corellation coefficients of 0
				for i, feat_fig_of_merit in enumerate( [ line[0] for line in use_these_feature_weights ] ):
					if feat_fig_of_merit == 0:
						if _all == 'nonzero':
							print "Features rank index {0} and below have a correllation coefficient of 0. ".format( i )
							print 'Using {0} features'.format( i )
							use_these_feature_weights = use_these_feature_weights[ : i ]
							break
						else:
							err_msg = "Can't reduce feature weights \"{0}\" to {1} features. ".format( self.name, num_features_to_be_used )
							err_msg += "Features ranked {0} and below have a correllation coefficient of 0. ".format( i )
							err_msg += "Request less features. "
							raise ValueError( err_msg )

		# we want lists, not tuples!
		abs_corr_coeffs, new_weights.featurenames_list, new_weights.pearson_coeffs, new_weights.slopes, \
		    new_weights.intercepts, new_weights.pearson_stderrs, new_weights.pearson_p_values,\
		    new_weights.spearman_coeffs, new_weights. spearman_p_values =\
		      [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

		r_val_sum = 0
		for val in abs_corr_coeffs:
			r_val_sum += val * val
		new_weights.values = [ ( (val*val) / r_val_sum ) for val in abs_corr_coeffs ]

		new_weights.associated_training_set = self.associated_training_set

		return new_weights

	#================================================================
	def Slice( self, start_index, stop_index ):
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
		raw_featureweights = zip( self.featurenames_list, abs_val_pearson_coeffs, self.pearson_coeffs, \
		    self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
		    self.spearman_coeffs, self.spearman_p_values )

		use_these_feature_weights = \
				list( itertools.islice( raw_featureweights, min_index, max_index ) )
		
		new_weights.featurenames_list, abs_pearson_coeffs, new_weights.pearson_coeffs, new_weights.slopes, \
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
	def Print( self, print_legend=True ):
		"""@brief Prints out feature values and statistics"""

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
			if len( self.featurenames_list[i] ) < 50:
				line_item += self.featurenames_list[i]
			else:
				line_item += self.featurenames_list[i][:50] + '... (truncated)'
			print line_item


#############################################################################
# class definition of FeatureGroup
# N.B.: Now implemented in C++ with the FeatureGroup class (originally in FeatureNames.h/.cpp)
#############################################################################
# 
# 	A FeatureGroup is composed of: 1. A list of transforms to apply
# 	sequentially to an input image/ROI/pixel plane, and 2. The algorithm to apply to
# 	that transformed pixel plane which produces the image descriptor values.
#
# 	The member "algorithm" is a reference to the SWIG-wrapped C++ class FeatureAlgorithm.
# 	The member "transform_list" is a list of references to the SWIG-wrapped C++ class
# 	FeatureTransform.
# 	The member "name" is a string representation of the algorithm and transform list,
# 	following the convention "AlgorithmName ([Transform A ([Transform B (])])."""
#   e.g.: "Tamura Textures (Wavelet (Edge ()))"
#   These FeatureGroup names are keys to a static cache of FeatureGroups.
#
#================================================================
def GenerateComputationPlanFromListOfFeatureStrings( feature_list ):
	"""Takes list of feature strings and chops off bin number at the first space on right, e.g.,
	"feature alg (transform()) [bin]" ... Returns a FeatureComputationPlan

	@return work_order - a FeatureComputationPlan
	"""
	feature_plan = wndcharm.FeatureComputationPlan ('custom')

	feature_groups = set()
	for feature in feature_list:
		split_line = feature.rsplit( " ", 1 )
		# add to set to ensure uniqueness
		if split_line[0] not in feature_groups:
			feature_plan.add( split_line[0] )
			feature_groups.add( split_line[0] )
	feature_plan.finalize()

	return feature_plan

#############################################################################
# class definition of FeatureVector
#############################################################################
class FeatureVector( object ):
	"""
	FeatureVector contains features for a single image or ROI, as well as all the sampling
	options. It is the values contained here that are grouped together to form FeatureSpaces."""

	import re
	sig_filename_parser = re.compile( ''.join( [
	# basename:
	r'(?P<basename>.+?)',
	# ROI:
	r'(?P<roi>-B(?P<x>\d+)_(?P<y>\d+)_(?P<w>\d+)_(?P<h>\d+))?',
	# downsampling:
	r'(?:-d(?P<downsample>\d+))?',
	# pixel intensity adjustment:
	r'(?:-S(?P<pixel_intensity_mean>\d+)(?:_(?P<pixel_intensity_stddev>\d+))?)?',
	# rotations:
	r'(?:-R_(?P<rot>\d))?',
	# tiling info:
	r'(?P<tiling_scheme>-t(?P<tile_num_rows>\d+)(?:x(?P<tile_num_cols>\d+))?_(?P<tile_col_index>\d+)_(?P<tile_row_index>\d+))?',
	# color features:
	r'(?P<color>-c)?',
	# long feature set:
	r'(?P<long>-l)?',
	# extension: .sig or .pysig
	r'\.(?:py)?sig$' ] ) )

	#==============================================================
	def __init__( self, **kwargs ):

		#: Row name, common across all channels
		self.name = None

		#: feature names
		self.featurenames_list = None
		#: feature vector
		self.values = None
		#: feature_maxima, feature_minima, and normalized_against are members
		#: this object shares with the FeatureSpaceObject, as well as implementations
		#: of FeatureReduce() and Normalize(), and CompatibleFeatureSetVersion()
		self.feature_maxima = None
		self.feature_minima = None
		self.normalized_against = None
		#: the prefix string to which sampling options will be appended to form .sig filepath
		self.basename = None
		#: Can also be a reference to a wndchrm.ImageMatrix object
		self.source_filepath = None
		#: Path to .sig file, in future hdf/sql file
		self.auxiliary_feature_storage = None
		#: label is stringified ground truth
		self.label = None
		self.ground_truth = None
		self.samplegroupid = None
		self.channel = None
		self.time_index = None
		self.tiling_scheme = None
		#: If no ROI image subsample, whole image is tile 1 of 1 in this sample group.
		self.tile_num_rows = 1
		self.tile_num_cols = 1
		#: indices count from 0
		self.tile_row_index = 0
		self.tile_col_index = 0
		self.samplesequenceid = 0
		#: downsample (in percents)
		self.downsample = 0
		self.pixel_intensity_mean = None
		self.pixel_intensity_stddev = None
		self.roi = None
		self.h = None
		self.w = None
		self.x = None
		self.y = None
		self.z = None
		self.z_delta = None
		self.rot = None
		self.fs_col = 0

		#: self.num_features should always be == len( self.featurenames_list ) == len( self.values )
		self.num_features = None

		# WND-CHARM feature bank-specific params:
		self.color = None
		self.long = None
		self.feature_set_version = None
		self.feature_computation_plan = None

		self.Update( **kwargs )

	#==============================================================
	def __repr__( self ):
		return '"{0}" ({1}) grp {2}, idx {3}, col {4}'.format( self.name, self.label,
				self.samplegroupid, self.samplesequenceid, self.fs_col )

	#==============================================================
	def Update( self, force=False, **kwargs ):
		"""force - if True then kwargs coming in with None values will overwrite
		members in self, even if those members are non-None before calling Update."""

		self_namespace = vars( self )
		for key, val in kwargs.iteritems():
			# don't bother setting a val if it's None unless force
			if val is not None or force:  
				if key in self_namespace:
					#if self_namespace[ key ] != None and self_namespace[ key ] != val:
					#	from warnings import warn
					#	warn( "Overwriting attrib {0} old val {1} new val {2}".format( key, self_namespace[ key ], val ) )
					self_namespace[ key ] = val
				else:
					raise AttributeError( 'No instance variable named "{0}" in class {1}'.format(
						key, self.__class__.__name__ ) )

		# FIXME: feature_set_version refers to WND-CHARM specific feature set specifications.
		# Want to be able to handle other feature sets from other labs in the future.
		if self.feature_set_version == None:
			major = feature_vector_major_version

			# The feature_set_version helps describe what features contained in the set.
			# Major version has to do with fixing bugs in the WND_CHARM algorithm code base.
			# The minor version describes the composition of features in the feature set.
			# Minor versions 1-4 have specific combination of WND-CHARM features.
			# Minor version 0 refers to user-defined combination of features.

			# Check to see if there is a user-defined set of features for this feature vector:
			if self.featurenames_list:
				if len( self.featurenames_list ) not in feature_vector_minor_version_from_num_features:
					minor = 0
				else:
					# FIXME: If features are out of order, should have a minor version of 0
					minor = feature_vector_minor_version_from_num_features[ len( self.featurenames_list ) ]
			else:
				if not self.long:
					if not self.color:
						minor = 1
					else:
						minor = 3
				else:
					if not self.color:
						minor = 2
					else:
						minor = 4
			self.feature_set_version = '{0}.{1}'.format( major, minor )
		else:
			major, minor = [ int( val ) for val in self.feature_set_version.split('.') ]

		if self.num_features == None:
			if self.featurenames_list:
				self.num_features = len( self.featurenames_list )
			else:
				if major == 1:
					num_feats_dict = feature_vector_minor_version_from_num_features_v1
				else:
					num_feats_dict = feature_vector_minor_version_from_num_features
				for num_feats, version in num_feats_dict.iteritems():
					if version == minor:
						self.num_features = num_feats

		# When reading in sampling opts from the path, they get pulled out as strings
		# instead of ints:
		if self.tile_row_index and type( self.tile_row_index ) != int:
			self.tile_row_index = int( self.tile_row_index )
		if self.tile_col_index and type( self.tile_col_index ) != int:
			self.tile_col_index = int( self.tile_col_index )
		if self.tile_num_rows and type( self.tile_num_rows ) != int:
			self.tile_num_rows = int( self.tile_num_rows )
		if self.tile_num_cols and type( self.tile_num_cols ) != int:
			self.tile_num_cols = int( self.tile_num_cols )

		# sequence order:
		# index 0 = position row 0, col 0
		# index 1 = position row 0, col 1
		# index 2 = position row 1, col 0, etc...
		self.samplesequenceid = self.tile_row_index + ( self.tile_num_cols * self.tile_col_index )
		return self

	#==============================================================
	def Derive( self, **kwargs ):
		"""Make a copy of this FeatureVector, except members passed as kwargs"""

		from copy import deepcopy
		new_obj = self.__class__()
		self_namespace = vars( self )
		new_obj_namespace = vars( new_obj )

		# skip these (if any):
		convenience_view_members = []

		# Are all keys in kwargs valid instance attribute names?
		invalid_kwargs = set( kwargs.keys() ) - set( self_namespace.keys() )
		if len( invalid_kwargs ) > 0:
			raise ValueError( "Invalid keyword arg(s) to Derive: {0}".format( invalid_kwargs ) )

		# Go through all of self's members and copy them to new_fs
		# unless a key-val pair was passed in as kwargs
		for key in self_namespace:
			if key in convenience_view_members:
				continue
			if key in kwargs:
				new_obj_namespace[key] = kwargs[key]
			else:
				new_obj_namespace[key] = deepcopy( self_namespace[key] )
		return new_obj

	#==============================================================
	def __deepcopy__( self, memo ):
		"""Make a deepcopy of this FeatureVector"""
		return self.Derive()

	#==============================================================
	def GenerateSigFilepath( self ):
		"""The C implementation of wndchrm placed feature metadata
		in the filename in a specific order, recreated here."""

		# FIXME: sigpaths for FeatureVectors with different channels
		# may have sig file names that will collide/overwrite each other.
		if self.basename:
			base = self.basename
		elif self.source_filepath:
			from os.path import splitext
			base, ext = splitext( self.source_filepath )
		elif self.name:
			base = self.name
			self.basename = base
		else:
			raise ValueError( 'Need for "basename" or "source_filepath" or "name" attribute in FeatureVector object to be set to generate sig filepath.')

		self_namespace = vars(self)
		# ROI:
		roi_params = 'x', 'y', 'w', 'h',
		bb = [self_namespace[key] for key in roi_params]

		if all( [ val is not None for val in bb ] ):
			base += "-B{0}_{1}_{2}_{3}".format( *bb )
		if self.downsample:
			base += "-d" + str(self.downsample)
		if self.pixel_intensity_mean is not None:
			base += "-S" + str(self.pixel_intensity_mean)
			if self.pixel_intensity_stddev is not None:
				base += "-S" + str(self.pixel_intensity_stddev)
		if self.rot is not None:
			base += "-R_" + str(self.rot)
		if self.tile_num_rows and self.tile_num_rows != 1:
			base += "-t" + str(self.tile_num_rows)
			if self.tile_num_cols and self.tile_num_cols != 1:
				base +="x" + str(self.tile_num_cols)
			if self.tile_row_index is not None and self.tile_col_index is not None:
				base += "_{0}_{1}".format( self.tile_row_index, self.tile_col_index )
			else:
				raise ValueError('Need to specify tile_row_index and tile_col_index in self for tiling params')
		if self.color:
			base += '-c'
		if self.long:
			base += '-l'

		return base + '.sig'

	#================================================================
	def GenerateFeatures( self, write_sig_files_to_disk=True, quiet=True ):
		"""@brief Loads precalculated sigs, or calculates new ones, based on which instance
		attributes have been set, and what their values are.
		
		Returns self for convenience."""

		# 0: What features does the user want?
		# 1: are there features already calculated somewhere?
		# 2: if so are they current/complete/expected/correct?
		# 3: if not, what's left to calculate?
		# 4: Calculate the rest
		# 5: Reduce the features down to what the user asked for

		if self.values is not None and len( self.values ) != 0:
			return self

		# Make sure Feature Vector version string is correct, etc:
		self.Update()

		try:
			self.LoadSigFile( quiet=quiet )
			# FIXME: Here's where you'd calculate a small subset of features
			# and see if they match what was loaded from file. The file could be corrupted
			# incomplete, or calculated with different options, e.g., -S1441
			return self
		except IOError:
			pass

		# All hope is lost, calculate features.

		# Use user-assigned feature computation plan, if provided:
		if self.feature_computation_plan != None:
			comp_plan = self.feature_computation_plan
			self.feature_set_version = comp_plan.feature_vec_type
		else:
			major, minor = self.feature_set_version.split('.')
			if minor == '0':
				comp_plan = GenerateComputationPlanFromListOfFeatureStrings( self.featurenames_list )
			elif minor == '1':
				comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSet()
			elif minor == '2':
				comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetLong()
			elif minor == '3':
				comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetColor()
			elif minor == '4':
				comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetColorLong()
			else:
				raise ValueError( "Not sure which features you want." )
			self.feature_computation_plan = comp_plan

		# Here are the ImageMatrix API calls:
		# void normalize(double min, double max, long range, double mean, double stddev);
		# int OpenImage(char *image_file_name, int downsample, rect *bounding_rect, double mean, double stddev);
		# void Rotate (const ImageMatrix &matrix_IN, double angle);

		if self.rot is not None:
			raise NotImplementedError( "FIXME: Implement rotations." )

		self_namespace = vars( self )
		bb_members = 'x', 'y', 'w', 'h',
		bb_vals = tuple( [ self_namespace[ var ] for var in bb_members ] )
		if all( [ val is not None for val in bb_vals ] ):
			bb = wndcharm.rect( *bb_vals )
		else:
			bb = None

		if self.pixel_intensity_mean:
			mean = self.pixel_intensity_mean
		else:
			mean = 0

		if self.pixel_intensity_stddev:
			stddev = self.pixel_intensity_stddev
		else:
			stddev = 0

		if isinstance( self.source_filepath, str ):
			the_tiff = wndcharm.ImageMatrix()
			if 1 != the_tiff.OpenImage( self.source_filepath, self.downsample, bb, mean, stddev ):
				raise ValueError( 'Could not build an ImageMatrix from {0}, check the path.'.\
					format( self.source_filepath ) )
		elif isinstance( self.source_filepath, wndcharm.ImageMatrix ):
			the_tiff = self.source_filepath
			raise NotImplementedError( "FIXME: Still haven't implemented downsample, bounding box, pixel intensity mean and stddev for an already open instance of ImageMatrix." )
		else:
			raise ValueError("image parameter 'image_path_or_mat' is not a string or a wndcharm.ImageMatrix")

		# pre-allocate space where the features will be stored (C++ std::vector<double>)
		tmp_vec = wndcharm.DoubleVector( comp_plan.n_features )

		# Get an executor for this plan and run it
		plan_exec = wndcharm.FeatureComputationPlanExecutor( comp_plan )
		plan_exec.run( the_tiff, tmp_vec, 0 )

		# get the feature names from the plan
		comp_names = [ comp_plan.getFeatureNameByIndex(i) for i in xrange( comp_plan.n_features ) ]

		# convert std::vector<double> to native python list of floats
		comp_vals = list( tmp_vec )

		# Feature Reduction/Reorder step:
		# Feature computation may give more features than are asked for by user, or out of order.
		if self.featurenames_list:
			if self.featurenames_list != comp_names:
				self.values = [ comp_vals[ comp_names.index( name ) ] for name in self.featurenames_list ]
		else:
			self.featurenames_list = comp_names
			self.values = comp_vals

		if write_sig_files_to_disk:
			self.ToSigFile( quiet=quiet )

		# Feature names need to be modified for their sampling options.
		# Base case is that channel goes in the innermost parentheses, but really it's not
		# just channel, but all sampling options.
		# For now, let the FeatureSpace constructor code handle the modification of feature names
		# for its own self.featurenames_list
		return self

	#==============================================================
	def CompatibleFeatureSetVersion( self, version ):
		"""Note that if either minor version is 0 (i.e not a standard feature vector)
		we return true while in fact, the compatibility is unknown"""

		try:
			version = version.feature_set_version
		except AttributeError:
			# Assume version is the version string then.
			pass

		if self.feature_set_version is None or version is None:
			err_str = "Can't tell if FeatureSpace {0} is compatible with version {1} because "
			if self.feature_set_version is None and version is None:
				err_str += "both are null."
			elif self.feature_set_version is None:
				err_str += "the FS instance's version string is null."
			else:
				err_str += "input version string is null."
			raise AttributeError( err_str.format( self.name, version ) )

		their_major, their_minor = [ int(v) for v in version.split('.',1) ]
		our_major, our_minor = [ int(v) for v in self.feature_set_version.split('.',1) ]

		if( their_major != our_major ):
			return False
		if our_minor and their_minor and our_minor != their_minor:
			return False

		return True

	#==============================================================
	def Normalize( self, input_feat_container=None, inplace=True, quiet=False ):
		"""By convention, the range of feature values in the WND-CHARM algorithm are
		normalized on the interval [0,100]. Normalizing is useful in making the variation 
		of features human readable. Normalized samples are only comprable if they've been 
		normalized against the same feature maxima/minima."""

		if self.normalized_against:
			# I've already been normalized, and you want to normalize me again?
			raise ValueError( "{0} \"{1}\" has already been normalized against {2}.".format (
				self.__class__.__name__, self.name, self.normalized_against ) )

		newdata = {}

		if not input_feat_container:
			# Specific to FeatureVector implementation:
			# Doesn't make sense to Normalize a 1-D FeatureVector against itself
			# The FeatureSpace implementation of this function has stuff in this block
			err = "Can't normalize {0} \"{1}\" against itself (Normalize() called with blank arg)."
			raise ValueError( err.format( self.__class__.__name__, self.name ) )
		else:
			# Recalculate my feature space according to maxima/minima in input_feat_container
			if input_feat_container.featurenames_list != self.featurenames_list:
				err_str = "Can't normalize {0} \"{1}\" against {2} \"{3}\": Features don't match.".format(
				  self.__class__.__name__, self.name,
					input_feat_container.__class__.__name__, input_feat_container.name )
				raise ValueError( err_str )
			if not self.CompatibleFeatureSetVersion( input_feat_container ):
				err_str = 'Incompatible feature versions: "{0}" ({1}) and "{2}" ({3})'
				raise ValueError( err_str.format( self.name, self.feature_set_version,
				    input_feat_container.name, input_feat_container.feature_set_version ) )

			if not quiet:
				# Specific to FeatureVector implementation:
				# no num_samples member:
				print 'Normalizing {0} "{1}" ({2} features) against {3} "{4}"'.format(
					self.__class__.__name__, self.name, len( self.featurenames_list),
					input_feat_container.__class__.__name__, input_feat_container.name )

			# Need to make sure there are feature minima/maxima to normalize against:
			if not input_feat_container.normalized_against:
				input_feat_container.Normalize( quiet=quiet )

			mins = input_feat_container.feature_minima
			maxs = input_feat_container.feature_maxima
			newdata['normalized_against'] = input_feat_container

		newdata['values'] = np.copy( self.values )
		newdata['feature_minima'], newdata['feature_maxima'] = \
		    normalize_by_columns( newdata['values'], mins, maxs )

		if inplace:
			self.Update( **newdata )
		else:
			return self.Derive( **newdata )

	#==============================================================
	def FeatureReduce( self, requested_features, inplace=False ):
		"""Returns a new FeatureVector that contains a subset of the data by dropping
		features (columns), and/or rearranging columns.

		requested_features := an object with a "featurenames_list" member
		    (FeatureVector/FeatureSpace/FeatureWeights) or an iterable containing
		    strings that are feature names.

		Implementation detail: compares input "requested_features" to self.featurenames_list,
		and "requested_features" becomes the self.featurenames_list of the returned FeatureVector."""

		try:
			requested_features = requested_features.featurenames_list
		except AttributeError:
			# assume it's already a list then
			pass

		# Check that self's faturelist contains all the features in requested_features
		selfs_features = set( self.featurenames_list )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = "Feature Reduction error:\n"
			err_str += '{0} "{1}" is missing '.format( self.__class__.__name__, self.name )
			err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
					len( missing_features_from_req ), len( requested_features ) )
			err_str += "\nDid you forget to convert the feature names into their modern counterparts?"
			raise ValueError( err_str )

		# The implementation of FeatureReduce here is similar to FeatureSpace.FeatureReduce
		# Here is where the implementations diverge"
		num_features = len( requested_features )

		newdata = {}
		newdata[ 'name' ] = self.name + "(feature reduced)"
		newdata[ 'featurenames_list' ] = requested_features
		newdata[ 'num_features' ] = num_features

		new_order = [ self.featurenames_list.index( name ) for name in requested_features ]

		# N.B. 1-D version used here, contrast with FeatureSpace.FeatureReduce() implementation.
		newdata[ 'values' ] = self.values[ new_order ]

		if self.feature_maxima is not None:
			newdata[ 'feature_maxima' ] = self.feature_maxima[ new_order ]
		if self.feature_minima is not None:
			newdata[ 'feature_minima' ] = self.feature_minima[ new_order ]

		# If the feature vectors sizes changed then they are no longer standard feature vectors.
		if self.feature_set_version is not None and num_features != self.num_features:
			newdata[ 'feature_set_version' ] = \
					"{0}.0".format( self.feature_set_version.split('.',1)[0] )

		if inplace:
			return self.Update( **newdata )
		return self.Derive( **newdata )
	#================================================================
	def LoadSigFile( self, sigfile_path=None, quiet=False ):

		if sigfile_path:
			path = sigfile_path
		elif self.auxiliary_feature_storage:
			path = self.auxiliary_feature_storage
		else:
			path = self.GenerateSigFilepath()

		with open( path ) as infile:
			lines = infile.read().splitlines()

		import re
		# First line is metadata
		self.class_id, self.feature_set_version = \
		        re.match( '^(\S+)\s*(\S+)?$' , lines[0] ).group( 1, 2 )
		if self.feature_set_version is None:
			self.feature_set_version = "1.0"

		# 2nd line is path to original tiff file, just skip

		# Process rest of lines
		values, names = zip( *[ line.split( None, 1 ) for line in lines[2:] ] )

		# np.fromstring is a 3x PIG:
		# %timeit out = np.array( [ float(val) for val in thing ] )
		# 10 loops, best of 3: 38.3 ms per loop
		# %timeit out = np.fromstring( " ".join( thing ), sep=" " )
		# 10 loops, best of 3: 98.1 ms per loop

		self.values = np.array( [ float( val ) for val in values ] )

		# We would know by know if there was a sigfile processing error,
		# e.g., file doesn't exist.
		# Safe to set this member now if not already set
		if not self.auxiliary_feature_storage:
			self.auxiliary_feature_storage = path

		# Subtract path so that path part doesn't become part of name
		from os.path import basename
		# Pull sampling options from filename
		path_removed = basename( path )
		self.name = path_removed
		result = self.sig_filename_parser.search( path_removed )
		if result:
			self.Update( **result.groupdict() )

		# This is really slow:
		#for i, name in enumerate( names ):
			#retval = wndcharm.FeatureNames.getFeatureInfoByName( name )
			#if retval:
			#	self.featurenames_list[i] = retval.name
			#else:
			# self.featurenames_list[i] = name

		# Use pure Python for old-style name translation
		#from wndcharm import FeatureNameMap
		#self.featurenames_list = FeatureNameMap.TranslateToNewStyle( featurenames_list )

		# Deprecate old-style naming support anyway, those features are pretty buggy
		# -CEC 20150104
		self.featurenames_list = list( names )

		# Cleanup for legacy edge case:
		# Set the minor version to the vector type based on # of features
		# The minor versions should always specify vector types, but for version 1 vectors,
		# the version is not written to the file, so it gets read as 0.
		if( self.feature_set_version == "1.0" ):
			self.version = "1." + str(
				feature_vector_minor_version_from_num_features_v1.get( len( self.values ),0 ) )

		if not quiet:
			print "Loaded features from file {0}".format( path )
		return self

	#================================================================
	@classmethod
	def NewFromSigFile( cls, sigfile_path, image_path=None, quiet=False ):
		"""@return  - An instantiated FeatureVector class with feature names translated from
		           the old naming convention, if applicable."""
		return cls( source_filepath=image_path ).LoadSigFile( sigfile_path, quiet )

	#================================================================
	def ToSigFile( self, path=None, quiet=False ):
		"""Write features C-WND-CHARM .sig file format

		If filepath is specified, you get to name it whatever you want and put it
		wherever you want. Otherwise, it's named according to convention and placed 
		next to the image file in its directory."""
		from os.path import exists
		if not path:
			path = self.GenerateSigFilepath()
		if not quiet:
			if exists( path ):
				print "Overwriting {0}".format( path )
			else:
				print 'Writing signature file "{0}"'.format( path )
		self.auxiliary_feature_storage = path
		with open( path, "w" ) as out:
			# FIXME: line 1 contains class membership and version
			# Just hardcode the class membership for now.
			out.write( "0\t{0}\n".format( self.feature_set_version ) )
			out.write( "{0}\n".format( self.source_filepath ) )
			for val, name in zip( self.values, self.featurenames_list ):
				out.write( "{0:0.6g} {1}\n".format( val, name ) )

# end definition class FeatureVector

#############################################################################
# class definition of FeatureSpace
#############################################################################
class FeatureSpace( object ):
	"""An instance of FeatureSpace is one-half of a WND-CHARM classifier, the other half being the
	FeatureWeights instance.

	The FeatureSpace class is a container for sets of image descriptors, which are collected
	into Numpy matrices organized by image class or ground truth. FeatureSpaces are also used
	as containers for test images which have yet to be classified. FeatureSpaces can also be
	randomly Split() into two subset FeatureSpaces to be used for cross-validation of a
	classifier, one for training the other for testing.
	"""

	# Used for parsing "File of Files" definition of FeatureSpace
	import re
	channel_col_finder = re.compile(r'(?P<path>.*?)?\{(?P<opts>.*?)?\}')
	channel_opt_finder = re.compile(r'(?:(?P<key>.+?)=)?(?P<value>.+)')

	#==============================================================
	def __init__( self, name=None, source_filepath=None, num_samples=None, num_samples_per_group=1,
					num_features=None, discrete=True, feature_set_version=None ):
		"""FeatureSpace constructor"""
		
		# Let F = # features for a given 5D ROI.
		# Let S = total # samples (rows) in a feature set.
		# Let C = # of discrete classes, if any, for a classification problem.
		# Let Si = # of samples in a given class whose index is i.
		# Let G = # tiles/ROIS in a sample group. Is 1 if no tiling.

		# BASIC DATA MEMBERS
		# -------------------------------------------
		#: type: string
		self.name = name

		#: type: string
		#: Path to FeatureSpace source/definition/origination file or directory.
		self.source_filepath = source_filepath

		#: type: FeatureSpace, or string containing 'self'
		#: By convention, the range of values are normalized on an interval [0,100].
		#: Reference to self or another FeatureSpace, indicating source of feature
		#: maxima/minima to transform feature space to normalized interval.
		self.normalized_against = None

		#: type: numpy.ndarray
		#: 2D numpy matrix with shape=(F,S) that contains all features.
		self.data_matrix = None
		#: If classification, per-class views into the feature matrix
		self.data_list = None

		#: @type: boolean
		#: Set to True when features packed into single matrix via internal
		self.data_matrix_is_contiguous = False
		
		#: The feature vector version contained in this FeatureSpace
		#: The major version must match for all feature vectors in the FeatureSpace
		#: The minor version must match also if it is one of the standard feature vectors (i.e. non-0)
		self.feature_set_version = feature_set_version

		#: A string keep track of all the options (-l -S###, -t etc)
		#: FIXME: expand to have all options kept track of individually
		self.feature_options = None
		self.tile_rows = None
		self.tile_cols = None

		#: Do the samples belong to discrete classes for supervised learning, or not
		#: (regression, clustering)
		self.discrete = discrete

		#: shape - supposed to work like numpy.ndarray.shape
		self.shape = None

		# SAMPLE METADATA DATA MEMBERS
		# N.B. Here data members are grouped into sets of two, the first being the "_contiguous"
		# version which is a simple list of length S, and the second being compound lists of lists,
		# the outer list being length C, inner lists of length Si, which are per-class convenience
		# views into their _contiguous data member counterpart.
		# -------------------------------------------

		#: FIXME: Eliminate in favor of len( self.samplenames_list )
		self.num_samples = num_samples

		#: A list of sample names in same row order as their samples appear in self.data_matrix.
		#: Corresponding lists of lists of sample names grouped by view.

		self._contiguous_samplenames_list = None
		self.samplenames_list = None

		#: By default, samples are independent/not grouped for splitting purposes
		self.num_samples_per_group = num_samples_per_group

		#: Keeps track of which samples are grouped together and must not be separated
		#: when FeatureSpace is split for cross-validation purposes
		self._contiguous_samplegroupid_list = None
		self.samplegroupid_list = None

		#: An intra-sample group tile/ROI index, max value = G
		self._contiguous_samplesequenceid_list = None
		self.samplesequenceid_list = None

		#: A list of floats with is the "target" vector for regression, interpolation, etc.
		self._contiguous_ground_truths = None
		self.ground_truths = None

		# List data members for discrete data whose len() is the number of classes
		#: List of strings which are the class names
		self.classnames_list = None
		#: float-ified versions of class labels, if applicable
		self.interpolation_coefficients = None
		#: Number of samples in each class
		self.classsizes_list = None


		# FEATURE METADATA DATA MEMBERS
		# -------------------------------------------
		#: FIXME: Eliminate in favor of len( self.featurenames_list )
		self.num_features = num_features

		#: block out some features for purposes of feature contribution analysis, et al.
		self.feature_mask = None

		#: A list of strings length M
		self.featurenames_list = None

		#: Contains pre-normalized feature maxima so feature space of this or other
		#: FeatureSpaces can be transformed.
		self.feature_maxima = None

		#: Contains pre-normalized feature minima so feature space of this or other
		#: FeatureSpaces can be transformed.
		self.feature_minima = None

		### Now initialize array-like members if possible:
		if self.num_samples and self.num_features:
			self.shape = ( self.num_samples, self.num_features )
			self.data_matrix = np.empty( self.shape, dtype='double' )

		if self.num_samples:
			self._contiguous_samplenames_list = [None] * self.num_samples
			self._contiguous_samplegroupid_list = [None] * self.num_samples
			self._contiguous_samplesequenceid_list = [None] * self.num_samples
			self._contiguous_ground_truths = [None] * self.num_samples

		if self.num_features:
			self.featurenames_list = [None] * self.num_features

	#==============================================================
	def Derive( self, **kwargs ):
		"""Make a copy of this FeatureSpace, except members passed as kwargs"""

		from copy import deepcopy
		new_obj = self.__class__()
		self_namespace = vars( self )
		new_obj_namespace = vars( new_obj )

		# Don't bother copying these "view" members which are rebuilt by self._RebuildViews()
		convenience_view_members = [ 'data_list', 'samplenames_list', 'samplegroupid_list',\
		    'samplesequenceid_list', 'ground_truths' ]

		# Are all keys in kwargs valid instance attribute names?
		invalid_kwargs = set( kwargs.keys() ) - set( self_namespace.keys() )
		if len( invalid_kwargs ) > 0:
			raise ValueError( "Invalid keyword arg(s) to Derive: {0}".format( invalid_kwargs ) )

		# Go through all of self's members and copy them to new_fs
		# unless a key-val pair was passed in as kwargs
		for key in self_namespace:
			if key in convenience_view_members:
				continue
			if key in kwargs:
				new_obj_namespace[key] = kwargs[key]
			else:
				new_obj_namespace[key] = deepcopy( self_namespace[key] )
		new_obj._RebuildViews()
		return new_obj

	#==============================================================
	def Update( self, **kwargs ):
		"""Replace instance attribute values with the ones passed as kwargs."""

		# FIXME: Should we call self._RebuildViews() at the end every time?
		self_namespace = vars( self )
		for key, val in kwargs.iteritems():
			if key in self_namespace:
				#if self_namespace[ key ] != None and self_namespace[ key ] != val:
				#	from warnings import warn
				#	warn( "Overwriting attrib {0} old val {1} new val {2}".format( key, self_namespace[ key ], val ) )
				self_namespace[ key ] = val
			else:
				raise AttributeError( 'No instance variable named "{0}" in class {1}'.format(
				  key, self.__class__.__name__ ) )
		return self

	#==============================================================
	def __deepcopy__( self, memo ):
		"""Make a deepcopy of this FeatureSpace"""
		return self.Derive()

	#==============================================================
	@output_railroad_switch
	def Print( self, verbose=False ):
		"""Prints out basic attributes about this training set, including name, path to
		source data, number and composition of image classes, number of features, etc."""

		print 'Summary of {0} "{1}":'.format( self.__class__.__name__ , self.name )
		if self.name != self.source_filepath:
			print 'source: "{0}"'.format( self.source_filepath )
		print 'Total samples: {0} ({1} groups, {2} samples/group)'.format( self.num_samples,
		  len( set( self._contiguous_samplegroupid_list ) ), self.num_samples_per_group )
		print 'Total num features: {0}'.format( len( self.featurenames_list ) )
		print 'Feature Set Version: {0}'.format( self.feature_set_version )

		if self.discrete:
			rpt_str = '\tClass {0} "{1}": {2} samples ({3} groups)'
			if self.classnames_list is not None:
				for i, class_name in enumerate( self.classnames_list ):
					print rpt_str.format( i, class_name, len( self.samplenames_list[i] ),
							len( set( self.samplegroupid_list[i] ) ) )

		if verbose: # verbose implies print info for each sample
			if self.num_samples_per_group == 1:
				sample_metadata = \
				  zip( self._contiguous_samplenames_list, self._contiguous_ground_truths )
				header_str = "SAMP NAME\tGROUND TRUTH\n==============================================================="
				format_str = "{0}\t{1}"
			else:
				sample_metadata = zip( self._contiguous_samplenames_list, 
							self._contiguous_samplegroupid_list, self._contiguous_samplesequenceid_list,
							self._contiguous_ground_truths )
				header_str = "SAMP NAME\tGROUP INDEX\tTILE INDEX\tGROUND TRUTH\n===================================================================="
				format_str = "{0}\t{1:03d}\t{2:02d}\t{3}"

			print header_str
			for line_item in sample_metadata:
				print format_str.format( *line_item )
		print ""

	#==============================================================
	def __repr__( self ):
		"""Prints out basic attributes about this training set, including name, path to
		source data, number and composition of image classes, number of features, etc."""

		outstr = '<' +str( self.__class__.__name__ ) + ' '
		outstr += '"{0}"'.format( self.name ) + ' '
		outstr += 'n_features=' + str( self.num_features ) + ' '
		outstr += 'n_total_samples=' + str( self.num_samples )
		if self.discrete:
			outstr += ' n_classes=' + str( self.num_classes ) + ' '
			outstr += 'samples_per_class=(' + ', '.join( [ '"{0}": {1}'.format( name, quant ) \
							for name, quant in zip( self.classnames_list, self.classsizes_list ) ] ) + ')'
		outstr += '>'
		return outstr

	#==============================================================
	def CompatibleFeatureSetVersion( self, version ):
		"""Note that if either minor version is 0 (i.e not a standard feature vector)
		we return true while in fact, the compatibility is unknown"""

		try:
			version = version.feature_set_version
		except AttributeError:
			# Assume version is the version string then.
			pass

		if self.feature_set_version is None or version is None:
			err_str = "Can't tell if FeatureSpace {0} is compatible with version {1} because "
			if self.feature_set_version is None and version is None:
				err_str += "both are null."
			elif self.feature_set_version is None:
				err_str += "the FS instance's version string is null."
			else:
				err_str += "input version string is null."
			raise AttributeError( err_str.format( self.name, version ) )

		their_major, their_minor = [ int(v) for v in version.split('.',1) ]
		our_major, our_minor = [ int(v) for v in self.feature_set_version.split('.',1) ]

		if( their_major != our_major ):
			return False
		if our_minor and their_minor and our_minor != their_minor:
			return False

		return True

	#==============================================================
	@classmethod
	def NewFromPickleFile( cls, pathname ):
		"""Returns new instance of FeatureSpace build from a saved pickle file,
		with a filename ending in .fit.pickle"""

		path, filename = os.path.split( pathname )
		if filename == "":
			raise ValueError( 'Invalid pathname: {0}'.format( pathname ) )

		if not filename.endswith( ".fit.pickled" ):
			raise ValueError( 'Not a pickled FeatureSpace file: {0}'.format( pathname ) )

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

		if (the_training_set.feature_set_version is None):
			the_training_set.feature_set_version = "1." + str(
				feature_vector_minor_version_from_num_features_v1.get( 
					len(the_training_set.featurenames_list), 0 ) )

		return the_training_set

	#==============================================================
	def PickleMe( self, pathname=None ):
		"""Pickle this instance of FeatureSpace and write to file whose path is optionally
		specified by argument "pathname" """

		outfile_pathname = ""
		if pathname != None:
			outfile_pathname = pathname
		else:
			# try to generate a path based on member source_filepath
			if self.source_filepath == None or self.source_filepath == "":
				raise ValueError( "Can't pickle this training set: its 'source_filepath' member"\
						"is not defined, and you did not specify a file path for the pickle file." )
			if os.path.isdir( self.source_filepath ):
				# this trainingset was generated from a directory
				# naming convention is /path/to/topleveldir/topleveldir-options.fit.pickled
				root, top_level_dir = os.path.split( self.source_filepath )
				if self.feature_options != None and self.feature_options != "":
					outfile_pathname = os.path.join( self.source_filepath, \
							                  top_level_dir + self.feature_options + ".fit.pickled" )
				else:
					outfile_pathname = os.path.join( self.source_filepath, \
					                      top_level_dir + ".fit.pickled" )
			else:
				# was genearated from a file, could have already been a pickled file
				if self.source_filepath.endswith( "fit.pickled" ):
					outfile_pathname = self.source_filepath
				elif self.source_filepath.endswith( ".fit" ):
					outfile_pathname = self.source_filepath + ".pickled"
				else:
					outfile_pathname = self.source_filepath + ".fit.pickled"	

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
	def NewFromFitFile( cls, pathname, discrete=True, quiet=False,
		    global_sampling_options=None, **kwargs ):
		"""Helper function which reads in a c-chrm fit file.

		tile_options - an integer N -> NxN tile scheme, or a tuple (N,M) -> NxM tile scheme
		discrete_data - if false, try to interpret classes as continuous variable."""

		if not global_sampling_options:
			global_sampling_options = FeatureVector( **kwargs )

		import re
		from os.path import basename, split

		path, filename = split( pathname )
		if not filename.endswith( ".fit" ):
			raise ValueError( 'Not a .fit file: {0}'.format( pathname ) )

		if not quiet:
			print "Creating Training Set from legacy WND-CHARM text file file {0}".format( pathname )
		new_fs = cls()

		new_fs.source_filepath = filename
		new_fs.name = basename( filename )
		new_fs.discrete = discrete

		fitfile = open( pathname )

		name_line = False
		line_num = 0
		sample_count = 0

		new_fs.tile_rows = global_sampling_options.tile_num_rows
		new_fs.tile_cols = global_sampling_options.tile_num_cols
		new_fs.num_samples_per_group = new_fs.tile_rows * new_fs.tile_cols
		new_fs.global_sampling_options = global_sampling_options

		for line in fitfile:
			if line_num is 0:
				# 1st line: number of classes and feature vector version
				num_classes, feature_set_version = re.match('^(\S+)\s*(\S+)?$', line.strip()).group(1, 2)
				if feature_set_version is None:
					feature_set_version = "1.0"
				new_fs.feature_set_version = feature_set_version
				num_classes = int( num_classes )
				new_fs.num_classes = num_classes
				new_fs.classsizes_list = [0] * num_classes
				new_fs.classnames_list = [0] * num_classes

			elif line_num is 1:
				# 2nd line: num features
				num_features = int( line )
				new_fs.num_features = num_features
				new_fs.featurenames_list = [None] * num_features
				if( feature_set_version == "1.0" ):
					feature_set_version = "1." + str(
						feature_vector_minor_version_from_num_features_v1.get ( num_features,0 ) )
					new_fs.feature_set_version = feature_set_version

			elif line_num is 2:
				# 3rd line: number of samples
				num_samples = int( line )
				new_fs.num_samples = num_samples
				new_fs.shape = ( num_samples, num_features )
				new_fs.data_matrix = np.empty( new_fs.shape, dtype='double' )
				new_fs._contiguous_samplenames_list = [None] * num_samples

			elif line_num < ( num_features + 3 ):
				# Lines 4 through num_features contains the feature names
				#retval = wndchrm.FeatureNames.getFeatureInfoByName( line.strip() )
				#name = retval.name if retval else line.strip()
				name = line.strip()
				new_fs.featurenames_list[ line_num - 3 ] = name

			elif line_num == ( num_features + 3 ):
				# The line after the block of feature names is blank
				pass

			elif line_num < ( num_features + 4 + num_classes ):
				# Class labels
				class_index = line_num - num_features - 4
				new_fs.classnames_list[ class_index ] = line.strip()

			else:
				# Everything else after is a feature or a sample name
				# Comes in alternating lines of data, then path to sample original file (tif or sig)
				if not name_line:
					# strip off the class identity value, which is the last in the array
					features_string, class_index_string  = line.strip().rsplit( " ", 1 )
					class_index = int( class_index_string ) - 1
					new_fs.classsizes_list[ class_index ] += 1
					# np.fromstring is a PIG, see timeit data elsewhere in code

					new_fs.data_matrix[ sample_count ] = \
					  np.array( [ float(val) for val in features_string.split() ] )
				else:
					new_fs._contiguous_samplenames_list[ sample_count ] = line.strip()
					sample_count += 1
				name_line = not name_line

			line_num += 1

		fitfile.close()

		_retval = CheckIfClassNamesAreInterpolatable( new_fs.classnames_list )
		if _retval:
			# Numeric ground truth/target vector
			new_fs.interpolation_coefficients = _retval
			new_fs._contiguous_ground_truths = [ _retval[ class_index ] \
			    for class_index in xrange( num_classes ) \
			      for i in xrange( new_fs.classsizes_list[ class_index ] ) ]
		else:
			# Just a string label ground truth
			new_fs._contiguous_ground_truths = [ new_fs.classnames_list[ class_index ] \
			  for class_index in xrange( num_classes ) \
			    for i in xrange( new_fs.classsizes_list[ class_index ] ) ]

		if new_fs.num_samples_per_group != 1:
			# sample sequence id = tile id
			# goes: [ 1, 2, 3, 4, 1, 2, 3, 4, ... ]
			new_fs._contiguous_samplesequenceid_list = [ i \
			  for j in xrange( num_samples/new_fs.num_samples_per_group ) \
			    for i in xrange( new_fs.num_samples_per_group ) ]
			# samples with same group id can't be split
			# goes: [ 1, 1, 1, 1, 2, 2, 2, 2, ... ]
			new_fs._contiguous_samplegroupid_list = [ j \
			  for j in xrange( num_samples/new_fs.num_samples_per_group ) \
			    for i in xrange( new_fs.num_samples_per_group ) ]
		else:
			new_fs._contiguous_samplesequenceid_list = [1] * num_samples
			new_fs._contiguous_samplegroupid_list = range( num_samples )

		print "Features version from .fit file: {0}".format( new_fs.feature_set_version )
		new_fs._RebuildViews()
		return new_fs

	#==============================================================
	def ToFitFile( self, path ):
		"""Writes features to ASCII text file which can be read by classic wnd-charm."""
		fit = open( path, 'w' )

		# 1st line: number of classes and feature vector version
		fit.write( str(self.num_classes) + ' ' + self.feature_set_version + '\n' )
		# 2nd line: num features
		fit.write( str(self.num_features) + '\n' )
		# 3rd line: number of samples
		fit.write( str(self.num_samples) + '\n' )
		# Lines 4 through num_features contains the feature names
		for name in self.featurenames_list:
			fit.write( name + '\n' )
		# The line after the block of feature names is blank
		fit.write('\n')
		# Then all the Class labels
		for label in self.classnames_list:
			fit.write( label + '\n' )

		# In the fit file format, a sample's class membership is denoted by the final int
		# at the end of the line of features. A class index of 0 implies it belongs
		# to the UNKNOWN CLASS so in practical terms, fit file indexing starts at 1.
		if self.interpolation_coefficients is None and self.classnames_list is None:
			# For classless data, assign all samples to the UNKNOWN CLASS.
			class_indices = [0] * self.num_samples
		else:
			if self.interpolation_coefficients:
				ground_truth_class_vals = self.interpolation_coefficients
			else:
				ground_truth_class_vals = self.classnames_list

			class_indices = [None] * self.num_samples
			for i in xrange( self.num_samples ):
				try:
					val = str( 1 + ground_truth_class_vals.index( self._contiguous_ground_truths[i] ) )
				except ValueError:
					val = str( 0 )
				class_indices[i] = val

		# Finally, alternating lines of features and paths to sample original file (tif or sig)
		for i, sample_name in enumerate( self._contiguous_samplenames_list ):
			self.data_matrix[i].tofile( fit, sep=' ', format='%g' )
			# add class index of sample to end of features line
			fit.write( ' ' + class_indices[i] + '\n' )
			fit.write( sample_name + '\n' )

		fit.close()

	#==============================================================
	def Save( self, pathname ):
		"""Writes features to HDF5 file."""
		raise NotImplementedError( "FIXME" )

	#==============================================================
	def _RebuildViews( self, reorder=False ):
		"""Construct self's data members into either A) lists of per-class lists of 
		features/meature metadata which are optimized for classification-style machine 
		learning problems or B) single contiguous lists of data for regression-style problems.

		reorder - a sample was added out of order, reorder by class membership"""

		if reorder:
			raise NotImplementedError( "Sorry, this method doesn't handle sorting by classes yet" )

		if self.discrete is None:
			errmsg = 'FeatureSpace {0} "discrete" member hasn\'t been set. '.format( self )
			errmsg += 'Please set the flag on the object indicating classification vs. regression/clustering.'
			raise ValueError( errmsg )

		if self.discrete == True:
			self.data_list = [None] * self.num_classes
			self.samplenames_list = [None] * self.num_classes
			self.samplegroupid_list = [None] * self.num_classes
			self.samplesequenceid_list = [None] * self.num_classes
			if self._contiguous_ground_truths:
				self.ground_truths = [None] * self.num_classes

			class_bndry_index = 0
			for class_index in xrange( self.num_classes ):
				n_class_samples = self.classsizes_list[ class_index ]
				self.data_list[ class_index ] = \
					self.data_matrix[ class_bndry_index : class_bndry_index + n_class_samples ]
				self.samplenames_list[ class_index ] = \
					self._contiguous_samplenames_list[ class_bndry_index : class_bndry_index + n_class_samples ]
				self.samplegroupid_list[ class_index ] = \
					self._contiguous_samplegroupid_list[ class_bndry_index : class_bndry_index + n_class_samples ]
				self.samplesequenceid_list[ class_index ] = \
					self._contiguous_samplesequenceid_list[ class_bndry_index : class_bndry_index + n_class_samples ]
				if self._contiguous_ground_truths:
					self.ground_truths[ class_index ] = \
						self._contiguous_ground_truths[ class_bndry_index : class_bndry_index + n_class_samples ]

				class_bndry_index += n_class_samples

		else:
			self.data_list = self.data_matrix
			self.samplenames_list = self._contiguous_samplenames_list
			self.samplegroupid_list = self._contiguous_samplegroupid_list 
			self.samplesequenceid_list = self._contiguous_samplesequenceid_list
			self.ground_truths = self._contiguous_ground_truths

		self.data_matrix_is_contiguous = True
		return self

	#==============================================================
	@classmethod
	def NewFromDirectory( cls, top_level_dir_path, discrete=True, num_samples_per_group=1,
	  quiet=False, global_sampling_options=None, write_sig_files_to_disk=True, **kwargs ):
		"""@brief Equivalent to the "wndchrm train" command from the C++ implementation by Shamir.
		Read the the given directory and parse its structure for class membership.
		Populate a list of FeatureVector instances, then call helper functions to
		load/calculate features and populate this object."""

		if not global_sampling_options:
			global_sampling_options = FeatureVector( **kwargs )

		if not quiet:
			print "Creating Training Set from directories of images {0}".format( top_level_dir_path )

		samples = []
		tile_num_rows = global_sampling_options.tile_num_rows
		tile_num_cols = global_sampling_options.tile_num_cols

		from copy import deepcopy
		from os import walk
		from os.path import join, basename
		sample_group_count = 0
		for root, dirs, files in walk( top_level_dir_path ):
			if root == top_level_dir_path:
				if len( dirs ) <= 0:
					# no class structure
					filelist = [ join( root, _file ) for _file in files \
					                  if _file.endswith(('.tif','.tiff','.TIF','.TIFF')) ]
					if len( filelist ) <= 0:
						raise ValueError( 'No tiff files in directory {0}'.format( root ) )
					for _file in filelist:
						tile_position_index = col_index * tile_num_rows + row_index
						for col_index in xrange( tile_num_cols ):
							for row_index in xrange( tile_num_rows ):
								fv = deepcopy( global_sampling_options )
								fv.source_filepath = _file
								fv.label = "UNKNOWN"
								fv.tile_row_index = row_index
								fv.tile_col_index = col_index
								fv.Update()
								samples.append( fv )
						sample_group_count += 1
					break
			else:
				# If we're here, we're down in one of the sub directories
				# This class's name will be "subdir" in /path/to/topleveldir/subdir
				filelist = [ join( root, _file ) for _file in files \
													if _file.endswith(('.tif','.tiff','.TIF','.TIFF')) ]
				if len( filelist ) <= 0:
					# Maybe the user hid them?
					#raise ValueError( 'No tiff files in directory {0}'.format( root ) )
					continue
				class_name = basename( root )
				for _file in filelist:
					tile_position_index = col_index * tile_num_rows + row_index
					for col_index in xrange( tile_num_cols ):
						for row_index in xrange( tile_num_rows ):
							fv = deepcopy( global_sampling_options )
							fv.source_filepath = _file
							fv.label = class_name
							fv.tile_row_index = row_index
							fv.tile_col_index = col_index
							fv.Update()
							samples.append( fv )
					sample_group_count += 1

		# FIXME: Here's where the parallization magic can (will!) happen.
		[ fv.GenerateFeatures( write_sig_files_to_disk, quiet ) for fv in samples ]

		name = basename( top_level_dir_path )
		return cls._NewFromListOfFeatureVectors( samples, name=name,
		       source_filepath=top_level_dir_path,
					 num_samples=None,
		       num_samples_per_group=(tile_num_rows*tile_num_cols),
					 num_features=global_sampling_options.num_features,
		       discrete=discrete, quiet=quiet )

	#==============================================================
	@classmethod
	def NewFromFileOfFiles( cls, pathname, discrete=True, quiet=False,
		     global_sampling_options=None, write_sig_files_to_disk=True, **kwargs ):
		"""Create a FeatureSpace from a file of files.

		The original FOF format (pre-2015) was just two columns, a path and a ground truth
		separated by a tab character. The extention to this format supports additional optional
		columns specifying additional paths and preprocessing options for a more complex
		feature space."""

		from os import getcwd
		from os.path import split, splitext, isfile, join
		from copy import deepcopy

		if not global_sampling_options:
			global_sampling_options = FeatureVector( **kwargs )

		if not quiet:
			print 'Loading {0} from file of files "{1}"'.format( cls.__name__, pathname )

		basepath, ext = splitext( pathname )
		dir_containing_fof, file_name = split( basepath )
		cwd = getcwd()

		num_fs_columns = None
		num_features = None
		feature_set_version = None

		# A FeatureSpace's samples in feature space need to be grouped by ground truth.
		# These variables help to determine if the samples listed in the input FOF are grouped.
		seen_ground_truths = set()
		previous_ground_truth = None
		current_ground_truth = None
		samples_grouped_by_ground_truth = True

		tile_num_rows = global_sampling_options.tile_num_rows
		tile_num_cols = global_sampling_options.tile_num_cols
		num_samples_per_group = tile_num_rows * tile_num_cols

		# Keeps track of the sample names to help organize like
		# samples into sample groups
		samp_name_to_samp_group_id_dict = {}

		def ReturnSampleGroupID( name ):
			if name not in samp_name_to_samp_group_id_dict:
				samp_name_to_samp_group_id_dict[ name ] = len(samp_name_to_samp_group_id_dict)
			return samp_name_to_samp_group_id_dict[ name ]

		samples = []
		fof = open( pathname )

		# There are a number of places to look for the existence of the files
		# specified via relative path.

		for line_num, line in enumerate( fof ):
			cols = line.strip().split('\t', 2)

			if samples_grouped_by_ground_truth:
				if current_ground_truth == None:
					# first line
					current_ground_truth = cols[1]
				else:
					previous_ground_truth = current_ground_truth
					current_ground_truth = cols[1]
					if current_ground_truth != previous_ground_truth:
						# crossed a class boundary, can't use previous ground truth anymore
						seen_ground_truths.add( previous_ground_truth )
					if current_ground_truth in seen_ground_truths:
						# Samples not grouped, will require sample sorting.
						# No need to keep checking.
						samples_grouped_by_ground_truth = False

			# Classic two-column FOF format
			if len( cols ) < 3:

				# If first time through, set a flag to indicate this is an classic version FOF
				# Set number of sample columns = -1 so as not to confuse with columns 0, 1, 2, ...
				# in multichannel FOF format below.
				if num_fs_columns == None:
					num_fs_columns = -1
				elif num_fs_columns != -1:
					err_smg = "File {0}, line {1} has old-style two-column FOF format, while the other lines use the new-style format with {3} columns"
					raise ValueError( err_msg.format( pathname, line_num, num_fs_columns + 3 ) )

				# Create a sampling opts template for this line in the FOF
				base_sample_opts = deepcopy( global_sampling_options )
				base_sample_opts.name = cols[0]
				base_sample_opts.label = cols[1]
				if not discrete:
					base_sample_opts.ground_truth = float(cols[1])
				# Note only difference with base_sampling_opts in the 3+ col version code below
				# is the fs_col is always 0
				base_sample_opts.fs_col = 0

				if not isfile( cols[0] ):
					# Try prepending current working directory
					path_to_sample = join( cwd, cols[0] )
					if not isfile( path_to_sample ):
						# Try prepending path to this FOF
						path_to_sample = join( dir_containing_fof, cols[0] )
						if not isfile( path_to_sample ):
							# Don't know what else to tell ya, pal.
							raise ValueError( "Can't find sample \"{0}\"".format( cols[0] ) )
				else:
					path_to_sample = cols[0]

				if path_to_sample.endswith('sig'):
					# Not possible to tile over a feature vector, so it's just one and done here.
					# Need to actually load the sig file here to tell how many features are contained
					base_sample_opts.LoadSigFile( path_to_sample )
					# Classic 2-col FOF listing sig path is an edgecase where we assign sample group id
					# based on base prefix resulting from parsing out the sampling options ("-l", etc).
					base_sample_opts.samplegroupid = ReturnSampleGroupID( base_sample_opts.basename )
					samples.append( base_sample_opts )
				else:
					base_sample_opts.source_filepath = path_to_sample
					base_sample_opts.samplegroupid = ReturnSampleGroupID( cols[0] )
					for col_index in xrange( tile_num_cols ):
						for row_index in xrange( tile_num_rows ):
							fv = deepcopy( base_sample_opts )
							fv.Update( tile_row_index=row_index, tile_col_index=col_index )
							samples.append( fv )
				# By now (after perhaps needing to load sig file) we know how many features in this sample
				num_feats_in_this_row = base_sample_opts.num_features
				if feature_set_version is None:
					feature_set_version = base_sample_opts.feature_set_version
				elif feature_set_version != base_sample_opts.feature_set_version:
					err_str = 'FOF line {0} has feature set version "{1}" that differs from rest "{1}"'
					raise ValueError( err_str.format( line_num, base_sample_opts.feature_set_version,
					    feature_set_version ) )
			# NEW THREE-column FOF format:
			# Within "column 3" are sub-columns indicating channel options.
			# Generates a dict of FeatureSpace processing options for each column of
			# channel options found, which takes the form:
			# [optional path to tiff or sig file] {[optname1=val1;optname2=val2...]}
			else:
				num_feats_in_this_row = 0
				# tabs are ignored here, now the {} becomes the delimiter
				tabs_stripped_out = cols[2].translate( None, '\t' )
				for fs_col, m in enumerate( cls.channel_col_finder.finditer( tabs_stripped_out ) ):

					# Start with a clean sample options template for each column
					base_sample_opts = deepcopy( global_sampling_options )
					base_sample_opts.name = cols[0]
					base_sample_opts.samplegroupid = ReturnSampleGroupID( cols[0] )
					base_sample_opts.label = cols[1]
					if not discrete:
						base_sample_opts.ground_truth = float(cols[1])
					base_sample_opts.fs_col = fs_col

					# Pull sampling options out of parentheses and load into template
					col_finder_dict = m.groupdict()
					opts_str = col_finder_dict['opts']
					if opts_str:
						col_opts = dict( [ mm.groups() for opt in opts_str.split(';') \
						                    for mm in cls.channel_opt_finder.finditer( opt ) ] )
						if None in col_opts:
							# value given without key is taken to be default channel
							col_opts['channel'] = col_opts[None]
							del col_opts[None]
						base_sample_opts.Update( **col_opts )

					# Now we deal with the path:
					if col_finder_dict['path'] == None:
						# If no path given in this channel column, implies string in cols[0]
						# is the path upon with the options refer to.
						# Doesn't make sense to have a sig file be in column 1
						# then have channel options listed, so cols[0] must be a path to a tiff
						path = cols[0]
					else:
						path = col_finder_dict['path']
					if not isfile( path ):
						# Try prepending current working directory
						path_to_sample = join( cwd, path )
						if not isfile( path_to_sample ):
							# Try prepending path to this FOF
							path_to_sample = join( dir_containing_fof, path )
							if not isfile( path_to_sample ):
								# Don't know what else to tell ya, pal.
								raise ValueError( "Can't find sample \"{0}\"".format( path_to_sample ) )
						path = path_to_sample

					if path.endswith('sig'):
						# Not possible to tile over a feature vector, so it's just one and done here
						base_sample_opts.LoadSigFile( path )
						samples.append( base_sample_opts )
					else:
						base_sample_opts.source_filepath = path
						for col_index in xrange( tile_num_cols ):
							for row_index in xrange( tile_num_rows ):
								fv = deepcopy( base_sample_opts )
								fv.Update( tile_row_index=row_index, tile_col_index=col_index )
								samples.append( fv )
					# By now (after perhaps needing to load sig file) we know how many features in this sample
					num_feats_in_this_row += base_sample_opts.num_features

				# END loop over {} columns
				# FIXME: Fix the number of channel column options to be uniform for now,
				# may be allowable in the future to have different number of channel opt columns
				# say for when the user wants to pull a ROI for feature calculation for a certain
				# sample, but use the whole image for all the others.
				if num_fs_columns == None:
					num_fs_columns = fs_col
				elif fs_col != num_fs_columns:
					err_smg = "File {0}, line {1} has {2} channel cols, when the rest has {3}"
					raise ValueError( err_msg.format( pathname, line_num, fs_col, num_fs_columns) )

				# FIXME: This is kinda kludgy since it doesn't evaluate major versions across columns
				if feature_set_version is None:
					# More than one col = must be nonstandard FS
					# Take major version sample opts from most recent col and call that this
					# FeatureSapce's major version
					if fs_col > 0:
						feature_set_version = base_sample_opts.feature_set_version.split('.',1)[0] + '.0'
					else:
						base_sample_opts.feature_set_version

			# END if len( cols ) < 3

			# Num features across entire row, no matter the columns, must be uniform
			if num_features == None:
				num_features = num_feats_in_this_row
			elif num_features != num_feats_in_this_row:
				errmsg = 'FOF "{0}" row {1} calls for {2} features where previous rows ' + \
					'call for {3} features.'
				raise ValueError(
					errmsg.format( pathname, line_num, num_feats_in_this_row, num_features ) )

		# END iterating over lines in FOF
		fof.close()

		# FIXME: Here's where the parallization magic can (will!) happen.
		[ fv.GenerateFeatures( write_sig_files_to_disk, quiet ) for fv in samples ]

		assert num_features > 0

		return cls._NewFromListOfFeatureVectors( samples, name=file_name, source_filepath=pathname,
		       num_samples=len(samp_name_to_samp_group_id_dict)*num_samples_per_group,
		       num_samples_per_group=num_samples_per_group, num_features=num_features,
		       feature_set_version=feature_set_version, discrete=discrete, quiet=quiet )

	#==============================================================
	@classmethod
	def _NewFromListOfFeatureVectors( cls, feature_vectors_list, name, source_filepath=None,
	  num_samples=None, num_samples_per_group=1, num_features=None, feature_set_version=None,
	  discrete=True, quiet=True ):
		"""Input is list of FeatureVectors WHOSE FEATURES HAVE ALREADY BEEN CALCULATED.

		The arguments to this method call are intentionally very similar to the args of
		FeatureSpace.__init__()"""

		new_fs = cls( name, source_filepath, num_samples, num_samples_per_group,
		             num_features, discrete, feature_set_version )

		new_fs.Print()

		# For compound samples, e.g., multichannel, need to know the column offsets.
		# key: col index, value: index in data_matrix demarking rightmost feature for this column
		feature_set_col_offset = {}

		# column enumeration starts at 0, and need to know where the 0th column has its
		# LEFT boundary, i.e. what is column 0 - 1's right most boundary (exclusive):
		feature_set_col_offset[-1] = 0

		# Sort list of FeatureVectors by column so we can fill in data_matrix from left to right.
		for fv in sorted( feature_vectors_list, key=lambda fv: fv.fs_col ):
			if fv.values is None:
				raise ValueError( "Calls to this method require features to have already been calculated." )

			col_left_boundary_index = feature_set_col_offset[ fv.fs_col - 1 ]
			col_right_boundary_index = col_left_boundary_index + fv.num_features
			row_index = (fv.samplegroupid * num_samples_per_group) + fv.samplesequenceid

			# Fill in column metadata if we've not seen a feature vector for this col before
			if fv.fs_col not in feature_set_col_offset:
				feature_set_col_offset[ fv.fs_col ] = col_right_boundary_index
				new_fs.featurenames_list[ col_left_boundary_index : col_right_boundary_index ] = \
				  [ name.replace( '()', '({0})'.format( fv.fs_col ) ) for name in fv.featurenames_list ]

			# Fill in row metadata with FeatureVector data from column 0 only
			if fv.fs_col == 0: # (fs_col member must be > 0 and cannot be None)
				new_fs._contiguous_samplenames_list[ row_index ] = fv.name
				new_fs._contiguous_samplegroupid_list[ row_index ] = fv.samplegroupid
				new_fs._contiguous_samplesequenceid_list[ row_index ] = fv.samplesequenceid
				if discrete:
					new_fs._contiguous_ground_truths[ row_index ] = fv.label
				else:
					new_fs._contiguous_ground_truths[ row_index ] = fv.ground_truth

			new_fs.data_matrix[ row_index, col_left_boundary_index : col_right_boundary_index ] = \
			  fv.values

		# alias:
		ground_truths = new_fs._contiguous_ground_truths

		if discrete:
			# Uniquify the sample group list, maintaining order of input sample group list.
			seen = set()
			seen_add = seen.add
			new_fs.classnames_list = [ x for x in ground_truths if not (x in seen or seen_add(x) ) ]
			new_fs.interpolation_coefficients = \
				CheckIfClassNamesAreInterpolatable( new_fs.classnames_list )
			new_fs.num_classes = len( new_fs.classnames_list )
			new_fs.classsizes_list = \
			  [ ground_truths.count( label ) for label in new_fs.classnames_list ]

		if not discrete or new_fs.interpolation_coefficients is not None:
			new_fs._contiguous_ground_truths = [ float( val ) for val in ground_truths ]

		new_fs._RebuildViews()

		if not quiet:
			new_fs.Print()

		return new_fs

	#==============================================================
	def ContiguousDataMatrix( self ):
		"""This method should be called to access the class data_matrix field.
		In the case where there are both views into this matrix (e.g. data_list) as well as the matrix itself,
		this method ensures that the data_matrix is a vstack of all data_lists
		After this call, the self.data_matrix field will be contiguous, and all of the views in data_list
		are consistent with it.
		"""
		if not self.discrete:
			self._contiguous_samplenames_list = self.samplenames_list
			self._contiguous_samplegroupid_list = self.samplegroupid_list
			self._contiguous_samplesequenceid_list = self.samplesequenceid_list
			self._contiguous_ground_truths = self.ground_truths
			return self.data_matrix

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

		if copy_class == len( self.data_list ):
			raise RuntimeError( "Internal error: ContiguousDataMatrix: none of the class views had their own data. Is the data_matrix contiguous already?" )

		# resize the matrix
		if self.data_matrix is not None:
			self.data_matrix.resize (self.num_samples, self.num_features)
		else:
			#print "called with empty data_matrix"
			self.data_matrix = np.empty ([ self.num_samples, self.num_features ], dtype='double')
			copy_class = 0
			copy_row = 0

		# In addition, keep a list of sample names corresponding to the 
		# rows in the contiguous feature matrix
		self._contiguous_samplenames_list = [ None ] * self.num_samples
		self._contiguous_samplegroupid_list = [ None ] * self.num_samples
		self._contiguous_samplesequenceid_list = [ None ] * self.num_samples
		self._contiguous_ground_truths = [ None ] * self.num_samples

		# We need to start copying at the first non-view class mat to the end.
		for class_index in range (copy_class, len (self.data_list)):
			#print "copy class"+str(class_index)
			nrows = self.data_list[class_index].shape[0]
			self.data_matrix[copy_row : copy_row + nrows] = np.copy (self.data_list[class_index])
			self.data_list[class_index] = self.data_matrix[copy_row : copy_row + nrows]
			self._contiguous_samplenames_list[copy_row : copy_row + nrows] = \
			                                             self.samplenames_list[class_index]
			self._contiguous_samplegroupid_list[copy_row : copy_row + nrows] = \
			                                             self.samplegroupid_list[class_index]
			self._contiguous_samplesequenceid_list[copy_row : copy_row + nrows] = \
			                                             self.samplesequenceid_list[class_index]
			self._contiguous_ground_truths[copy_row : copy_row + nrows] = \
			                                             self.ground_truths[class_index]
			copy_row += nrows

		self.data_matrix_is_contiguous = True
		return self.data_matrix

	#==============================================================
	def Normalize( self, input_feat_container=None, inplace=True, quiet=False ):
		"""By convention, the range of feature values in the WND-CHARM algorithm are
		normalized on the interval [0,100]. Normalizing is useful in making the variation 
		of features human readable. Normalized samples are only comprable if they've been 
		normalized against the same feature maxima/minima."""

		if self.normalized_against:
			# I've already been normalized, and you want to normalize me again?
			raise ValueError( "{0} \"{1}\" has already been normalized against {2}.".format (
				self.__class__.__name__, self.name, self.normalized_against ) )

		newdata = {}

		if not input_feat_container:
			# Recalculate my feature space using my own maxima/minima
			if not quiet:
				print 'Normalizing {0} "{1}" ({2} images) against self'.format(
					self.__class__.__name__, self.name, self.num_samples )
			mins = None
			maxs = None
			newdata['normalized_against'] = 'self'
		else:
			# Recalculate my feature space according to maxima/minima in input_feat_container
			if input_feat_container.featurenames_list != self.featurenames_list:
				err_str = "Can't normalize {0} \"{1}\" against {2} \"{3}\": Features don't match.".format(
				  self.__class__.__name__, self.name,
					input_feat_container.__class__.__name__, input_feat_container.name )
				raise ValueError( err_str )
			if not self.CompatibleFeatureSetVersion( input_feat_container ):
				err_str = 'Incompatible feature versions: "{0}" ({1}) and "{2}" ({3})'
				raise ValueError( err_str.format( self.name, self.feature_set_version,
				    input_feat_container.name, input_feat_container.feature_set_version ) )
			
			if not quiet:
				print 'Normalizing "{0}" ({1} samples) against "{2}" ({3} samples)'.format(
					self.name, self.num_samples, input_feat_container.name, input_feat_container.num_samples )

			# Need to make sure there are feature minima/maxima to normalize against:
			if not input_feat_container.normalized_against:
				input_feat_container.Normalize( quiet=quiet )

			mins = input_feat_container.feature_minima
			maxs = input_feat_container.feature_maxima
			newdata['normalized_against'] = input_feat_container

		newdata['data_matrix'] = np.copy( self.ContiguousDataMatrix() )
		newdata['feature_minima'], newdata['feature_maxima'] = \
		    normalize_by_columns( newdata['data_matrix'], mins, maxs )

		if inplace:
			self.Update( **newdata )
			self._RebuildViews()
			return self
		else:
			return self.Derive( **newdata )

	#==============================================================
	def FeatureReduce( self, requested_features, inplace=False ):
		"""Returns a new FeatureSpace that contains a subset of the data by dropping
		features (columns), and/or rearranging columns.

		requested_features := an object with a "featurenames_list" member
		    (FeatureVector/FeatureSpace/FeatureWeights) or an iterable containing 
		    strings that are feature names.

		Implementation detail: compares input "requested_features" to self.featurenames_list,
		and "requested_features" becomes the self.featurenames_list of the returned FeatureSpace."""

		try:
			requested_features = requested_features.featurenames_list
		except AttributeError:
			# assume it's already a list then
			pass

		# Check that self's faturelist contains all the features in requested_features
		selfs_features = set( self.featurenames_list )
		their_features = set( requested_features )
		if not their_features <= selfs_features:
			missing_features_from_req = their_features - selfs_features
			err_str = "Feature Reduction error:\n"
			err_str += '{0} "{1}" is missing '.format( self.__class__.__name__, self.name )
			err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
					len( missing_features_from_req ), len( requested_features ) )
			err_str += "\nDid you forget to convert the feature names into their modern counterparts?"
			raise ValueError( err_str )

		num_features = len( requested_features )
		shape = ( self.num_samples, num_features )

		newdata = {}
		newdata[ 'shape' ] = shape
		newdata[ 'source_filepath' ] = self.source_filepath + "(feature reduced)"
		newdata[ 'name' ] = self.name + "(feature reduced)"
		newdata[ 'featurenames_list' ] = requested_features
		newdata[ 'num_features' ] = num_features
		data_matrix = np.empty( shape , dtype='double' )

		# Columnwise operations in Numpy are a pig:
		# %timeit thing = shuffle_my_cols[:,desired_cols]
		# 1 loops, best of 3: 3.5 s per loop
		#
		# Yikes. Slightly faster to do the copying yourself:
		# %%timeit
		# thing = np.empty( (numrows, numcols/2))
		# for new_index, old_index in enumerate( desired_cols ):
		#    thing[ :, new_index ] = shuffle_my_cols[ :, old_index ]
		# 1 loops, best of 3: 2.25 s per loop

		new_order = [ self.featurenames_list.index( name ) for name in requested_features ]
		for new_index, old_index in enumerate( new_order ):
			data_matrix[ :, new_index ] = self.data_matrix[ :, old_index ]
		newdata[ 'data_matrix' ] = data_matrix

		if self.feature_maxima is not None:
			newdata[ 'feature_maxima' ] = self.feature_maxima[ new_order ]
		if self.feature_minima is not None:
			newdata[ 'feature_minima' ] = self.feature_minima[ new_order ]

		# If the feature vectors sizes changed then they are no longer standard feature vectors.
		if self.feature_set_version is not None and num_features != self.num_features:
			newdata[ 'feature_set_version' ] = \
					"{0}.0".format( self.feature_set_version.split('.',1)[0] )

		if inplace:
			return self.Update( **newdata )._RebuildViews()
		return self.Derive( **newdata )

	#==============================================================
	def SampleReduce( self, leave_in_samplegroupid_list=None, leave_out_samplegroupid_list=None ):
		"""Returns a new FeatureSpace that contains a subset of the data by dropping
		samples (rows), and/or rearranging rows.

		leave_in_sample_group_list := indicate the composition of the FeatureSpace to be returned.
			For discrete/classification FeatureSpaces:
		    an iterable of iterables of sample group indices indicating desired sample groups;
			For continuous/regression FeatureSpaces:
		    a iterable of desired sample group indices.

		leave_out_samplegroupid_list := a list containing sample group ids
		    that should be left out
		Returns a near-deep copy of self including only the sample groups specified in the list.
		If no tiles, sample group reduces to just sample index."""

		if leave_in_samplegroupid_list is None and leave_out_samplegroupid_list is None:
			raise ValueError( 'Invalid input, both leave_in_samplegroupid_list and leave_out_samplegroupid_list were None')

		if self.normalized_against:
			errmsg = 'Cannot perform SampleReduce on FeatureSpace "{0}" '.format( self.name ) + \
			    "because it's features have already been normalized and is therefore immutable."
			raise ValueError( errmsg )

		# Helper nested functions:
		#==================================
		def CheckForValidListOfInts( the_list ):
			"""Items in list must be ints that are valid sample group ids for this FeatureSpace"""
			for item in the_list:
				if type( item ) is not int:
					raise TypeError( "Input must be an int or a flat iterable containing only ints.")

			if not set( the_list ) < set( self._contiguous_samplegroupid_list ):
				msg = "Input contains sample group ids that aren't " + \
							'contained in FeatureSpace "' + self.name + '", specifically: ' + \
				      str( sorted( list( set( the_list ) - set( self._contiguous_samplegroupid_list ) ) ) )
				raise ValueError( msg )

		def CheckForValidLISTOFLISTSOfInts( the_list ):
			try:
				for item in the_list:
					if type( item ) is bool:
						continue
					elif type( item ) is list:
						CheckForValidListOfInts( item )
			except TypeError:
				raise TypeError( "Input must be an iterable containing either booleans or iterables containing only ints.")

		def UniquifySansLeaveOutList( sg_list, leave_out ):
			seen = set()
			seen_add = seen.add
			uniq_sgids = [ x for x in sg_list  if not( x in seen or seen_add(x) ) and ( x not in leave_out) ]
			return uniq_sgids
		#==================================

		if leave_out_samplegroupid_list is not None:
			if type( leave_out_samplegroupid_list ) is int:
				leave_out_samplegroupid_list = [ leave_out_samplegroupid_list ]
			CheckForValidListOfInts( leave_out_samplegroupid_list )

			# build up a leave IN list, excluding the SGids that the user indicated
			if self.discrete:
				leave_in_samplegroupid_list = []
				for class_sgid_list in self.samplegroupid_list:
					class_leave_in_sg_list = UniquifySansLeaveOutList( class_sgid_list, leave_out_samplegroupid_list )
					leave_in_samplegroupid_list.append( class_leave_in_sg_list  )
			else:
				leave_in_samplegroupid_list = \
				  UniquifySansLeaveOutList( self.samplegroupid_list, leave_out_samplegroupid_list )
		else: # user provided leave in list
			if self.discrete:
				CheckForValidLISTOFLISTSOfInts( leave_in_samplegroupid_list )
			else: # if continuous
				if type( leave_in_samplegroupid_list ) is int:
					leave_in_samplegroupid_list = [ leave_in_samplegroupid_list ]
				CheckForValidListOfInts( leave_in_samplegroupid_list )

		# Dummyproofing over.
		# Now we can count on the fact that leave_in_samplegroupid_list is defined,
		# either by the user or by the above code.

		# How many total training groups are requested?
		if self.discrete:
			try:
				total_num_sample_groups = \
					sum( len( class_list ) for class_list in leave_in_samplegroupid_list if class_list )
			except TypeError:
				errmsg = 'Leave in list for discrete FeatureSpaces has to be a list (of length ' + \
				         'num_classes) of lists of ' + \
				         'desired sample group ids. Did you mean to pass it in as the leave OUT list?'
				raise TypeError( errmsg )
		else:
			total_num_sample_groups = len( leave_in_samplegroupid_list )

		total_num_samples = total_num_sample_groups * self.num_samples_per_group
		shape =  (total_num_samples, self.num_features)

		newdata = {}
		newdata[ 'shape' ] = shape
		newdata[ 'source_filepath' ] = self.source_filepath + " (subset)"
		newdata[ 'name' ] = self.name + " (subset)"
		newdata[ 'num_samples' ] = total_num_samples
		data_matrix = np.empty( shape, dtype='double' )
		_contiguous_samplegroupid_list = [None] * total_num_samples
		_contiguous_samplenames_list = [None] * total_num_samples
		_contiguous_samplesequenceid_list = [None] * total_num_samples
		_contiguous_ground_truths = [None] * total_num_samples

		j = 0
		if self.discrete:
			# If there's a False in the list of lists instead of a list, skip the class whose
			# index is in the same position as the False's index.
			newdata['classsizes_list' ] = classsizes_list = \
			    [ self.num_samples_per_group * len(class_group_list) \
			        for class_group_list in leave_in_samplegroupid_list if class_group_list ]
			newdata[ 'num_classes' ] = num_classes = len( classsizes_list )

			# If user requests more classes than exists in self, that's ok, but you have to makeup
			# classnames. Throw a letter on the end of Class, and if they want more than
			# 26 classes, well they can inherit from this class and reimplement this function
			if num_classes <= self.num_classes:
				newdata[ 'classnames_list' ] = [ self.classnames_list[i] \
			      for i, num_groups in enumerate( leave_in_samplegroupid_list ) if num_groups ]
				if self.interpolation_coefficients:
					newdata[ 'interpolation_coefficients' ] = [ self.interpolation_coefficients[i] \
			      for i, num_groups in enumerate( leave_in_samplegroupid_list ) if num_groups ]
			else:
				newdata[ 'classnames_list' ] = [ "Class" + letter for i, letter in \
						zip( leave_in_samplegroupid_list, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' ) if i ]
				newdata[ 'interpolation_coefficients' ] = None
			for class_group_list in leave_in_samplegroupid_list:
				if not class_group_list:
					continue
				for samp_group_id in class_group_list:
					_contiguous_samplegroupid_list[ j : j + self.num_samples_per_group ] = \
					    [samp_group_id] * self.num_samples_per_group
					j += self.num_samples_per_group
			  
		else:
			for samp_group_id in leave_in_samplegroupid_list:
				_contiguous_samplegroupid_list[ j : j + self.num_samples_per_group ] = \
				    [samp_group_id] * self.num_samples_per_group
				j += self.num_samples_per_group

		assert( len( _contiguous_samplegroupid_list ) == total_num_samples )

		for i in xrange( 0, total_num_samples, self.num_samples_per_group ):
			groupid = _contiguous_samplegroupid_list[i]
			original_index = self._contiguous_samplegroupid_list.index( groupid )
			np.copyto( data_matrix[ i : i + self.num_samples_per_group ],
							self.data_matrix[ original_index : original_index + self.num_samples_per_group ] )
			_contiguous_samplenames_list[ i : i + self.num_samples_per_group ] = \
			   self._contiguous_samplenames_list[ original_index : original_index +  self.num_samples_per_group]
			_contiguous_samplesequenceid_list[ i : i + self.num_samples_per_group ] = \
			   self._contiguous_samplesequenceid_list[ original_index : original_index +  self.num_samples_per_group]
			_contiguous_ground_truths[ i : i + self.num_samples_per_group ] = \
			   self._contiguous_ground_truths[ original_index : original_index +  self.num_samples_per_group ]

		newdata[ 'data_matrix' ] = data_matrix
		newdata[ '_contiguous_samplenames_list' ] = _contiguous_samplenames_list
		newdata[ '_contiguous_samplesequenceid_list' ] = _contiguous_samplesequenceid_list
		newdata[ '_contiguous_ground_truths' ] = _contiguous_ground_truths
		newdata[ '_contiguous_samplegroupid_list' ] = _contiguous_samplegroupid_list
		return self.Derive( **newdata )

	#==============================================================
	def RemoveClass( self, class_index ):
		"""Virtual method."""
		raise NotImplementedError()

	#==============================================================
	def ScrambleGroundTruths( self ):
		"""Virtual method. Produce an instant negative control training set"""
		raise NotImplementedError()
	#==============================================================
	def NewFromSQLite( self ):
		"""Virtual method."""
		raise NotImplementedError()
	#==============================================================
	def NewFromHDF5( self ):
		"""Virtual method."""
		raise NotImplementedError()
	#==============================================================
	def Split( self, train_size=None, test_size=None, random_state=True,
					balanced_classes=True, quiet=False ):
		"""Used for dividing the current FeatureSpace into two subsets used for classifier
		cross-validation (i.e., training set and test set).

		Analogous to Scikit-learn's cross_validation.train_test_split().

		Parameters (stolen directly from scikit-learn's documentation):

		test_size : float, int, or None (default is None)
		            If float, should be between 0.0 and 1.0 and represent the proportion
		            of the dataset to include in the test split (rounded up). If int, 
		            represents the absolute number of test samples. If None, the value is
		            automatically set to the complement of the train size. If train size 
		            is also None, test size is set to 0.25.

		train_size : float, int, or None (default is None)
                If float, should be between 0.0 and 1.0 and represent the proportion
		            of the dataset to include in the train split (rounded down). If int,
		            represents the absolute number of train samples. If None, the value
		            is automatically set to the complement of the test size.

		random_state : int or RandomState
                If true, generate a new random split. If int or Pseudo-random number
		            generator state used for random sampling. If value evaluates to false,
		            then do not randomize, but take the first samples in the order
		            the occur in the FeatureSpace/class."""

		# Step 1: Determine composition of split classes, i.e.,
		# figure out how many images/samples goes into the train and test sets respectively.

		if random_state:
			from numpy.random import RandomState

			if random_state is True:
				from numpy.random import shuffle
			elif type( random_state ) is RandomState:
				shuffle = random_state.shuffle
			elif type( random_state ) is int:
				shuffle = RandomState( random_state ).shuffle
			else:
				raise ValueError( 'Arg random_state must be an instance of numpy.random.RandomState, an int, or the value True')

		training_set_only = None

		if test_size == None:
			test_size = 0.25

		# ----- BEGIN NESTED FUNCTION -----------------------------------------------
		def CalcTrainTestSampleGroupMembership( samplegroupid_list, _max=None ):
			"""Takes a list of sample group ids, and returns a two subset lists of ids whose lengths
			are determined by train_size and test_size.

			Raises ValueError if train_set/test_set params are not on interval 0<val<1 for float.

			train_size/test_size indicate number of sample groups desired, unless self.num_samples_per_group == 1,
			in which case they refer to samples.

			group_list = list containing sample group ids to be divided among training/test sets
			according to proportions prescribed in train_size & test_size (pulled from from outer scope)

			_max = Pare the list of sample group ids down to this length. Used to enforce balanced classes.
			
			Returns list of sample group ids to pass to SampleReduce()"""

			from math import ceil, floor

			# Uniquify the sample group list, maintaining order of input sample group list.
			seen = set()
			seen_add = seen.add
			unique_samplegroup_ids = [ x for x in samplegroupid_list if not (x in seen or seen_add(x) ) ]
			if _max and _max > len( unique_samplegroup_ids ):
				# Chop off the end
				unique_samplegroup_ids = unique_samplegroup_ids[ : _max ]

			num_samplegroups = len( unique_samplegroup_ids )

			# TEST SET SIZE:
			if type( test_size ) is int:
				if test_size < 0 or test_size >= num_samplegroups: # test sets of size 0 are allowed
					errmsg = 'Arg test_size ({0}) must be 0 <= test_size < {1} (# images/sample groups).'
					raise ValueError( errmsg.format( test_size, num_samplegroups ) )
				n_samplegroups_in_test_set = test_size
			elif type( test_size ) is float:
				if test_size < 0 or test_size >= 1: # test sets of size 0 are allowed
					errmsg = 'Arg test_size fraction {0} must be 0 <= 1.0'
					raise ValueError( errmsg.format( test_size ) )
				n_samplegroups_in_test_set = int( ceil( num_samplegroups * test_size ) )
			else:
				raise ValueError( 'Invalid input: test_size={0}'.format( test_size ) )

			# TRAIN SET SIZE:
			if type( train_size ) is int:
				if train_size < 0 or train_size >= num_samplegroups: # train sets of size 1 are allowed
					errmsg = 'Arg train_size ({0}) must be 0 <= train_size < {1} (# images/sample groups).'
					raise ValueError( errmsg.format( train_size, num_samplegroups ) )
				n_samplegroups_in_training_set = train_size
			elif type( train_size ) is float:
				if train_size < 0 or train_size >= 1: # train sets of size 1 are allowed
					errmsg = 'Arg train_size fraction {0} must be 0 <= 1.0'
					raise ValueError( errmsg.format( train_size ) )
				n_samplegroups_in_training_set = int( floor( num_samplegroups * train_size ) )
				if n_samplegroups_in_training_set <= 0:
					raise ValueError( 'Please input train_size fraction such that there aren\'t 0 images (or sample groups) in training set' )
			elif train_size == None:
				n_samplegroups_in_training_set = num_samplegroups - n_samplegroups_in_test_set
			else:
				raise ValueError( 'Invalid input: train_size={0}'.format( test_size ) )

			if( n_samplegroups_in_training_set + n_samplegroups_in_test_set ) > num_samplegroups:
					raise ValueError( 'User input specified train/test feature set membership contain more samples than are available.' )

			if random_state:
				shuffle( unique_samplegroup_ids )

			train_groups = sorted( unique_samplegroup_ids[ : n_samplegroups_in_training_set ] )

			if n_samplegroups_in_test_set:
				test_groups = sorted( unique_samplegroup_ids[ n_samplegroups_in_training_set : \
					        n_samplegroups_in_training_set + n_samplegroups_in_test_set ] )
			else:
				test_groups = False

			return train_groups, test_groups
		# ----- END NESTED FUNCTION -----------------------------------------------

		if not self.discrete: # If classless data:
			train_groups, test_groups = CalcTrainTestSampleGroupMembership( self.samplegroupid_list )

		else: # Discrete classes
			train_groups = []
			test_groups = []

			if balanced_classes:
				if self.num_samples_per_group > 1:
					num_groups_per_class = [ num / self.num_samples_per_group for num in self.classsizes_list ]
				else:
					num_groups_per_class = self.classsizes_list
				smallest_class_size = min( num_groups_per_class )
			else:
				smallest_class_size=None

			for class_index in xrange( self.num_classes ):
				try:
					class_train_groups, class_test_groups = \
					  CalcTrainTestSampleGroupMembership( self.samplegroupid_list[class_index], _max=smallest_class_size  )
				except ValueError as e:
					addl_msg = "Error with class index " + str(class_index) + \
					           '. For discrete FeatureSpaces (with classes), train_size and test_size' + \
					           ' are evaluated per-class. '
					raise ValueError( addl_msg + e.message )
				else:
					train_groups.append( class_train_groups )
					test_groups.append( class_test_groups )					

			if not any( test_groups ):
				training_set_only = True

		training_set = self.SampleReduce( train_groups )
		if not quiet:
			training_set.Print()
		if training_set_only:
			return training_set

		test_set = self.SampleReduce( test_groups )
		if not quiet:
			test_set.Print()
		return training_set, test_set

<<<<<<< HEAD

<<<<<<< HEAD
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

	#==============================================================
	def __init__( self, *args, **kwargs ):
		"""constructor"""

		super( FeatureSet_Discrete, self ).__init__( *args, **kwargs )

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

		if copy_class == len( self.data_list ):
			raise RuntimeError( "Internal error: ContiguousDataMatrix: none of the class views had their own data. Is the data_matrix contiguous already?" )

		# resize the matrix
		if self.data_matrix is not None:
			self.data_matrix.resize (self.num_images, self.num_features)
		else:
			#print "called with empty data_matrix"
			self.data_matrix = np.empty ([ self.num_images, self.num_features ], dtype='double')
			copy_class = 0
			copy_row = 0

		# In addition, keep a list of sample names corresponding to the 
		# rows in the contiguous feature matrix
		self._contiguous_samplenames_list = [ None ] * self.num_images
		self._contiguous_samplegroupid_list = [ None ] * self.num_images
		self._contiguous_samplesequenceid_list = [ None ] * self.num_images
		self._contiguous_ground_truths = [ None ] * self.num_images

		# We need to start copying at the first non-view class mat to the end.
		for class_index in range (copy_class, len (self.data_list)):
			#print "copy class"+str(class_index)
			nrows = self.data_list[class_index].shape[0]
			self.data_matrix[copy_row : copy_row + nrows] = np.copy (self.data_list[class_index])
			self.data_list[class_index] = self.data_matrix[copy_row : copy_row + nrows]
			self._contiguous_samplenames_list[copy_row : copy_row + nrows] = \
			                                             self.imagenames_list[class_index]
			self._contiguous_samplegroupid_list[copy_row : copy_row + nrows] = \
			                                             self.samplegroupid_list[class_index]
			self._contiguous_samplesequenceid_list[copy_row : copy_row + nrows] = \
			                                             self.samplesequenceid_list[class_index]
			self._contiguous_ground_truths[copy_row : copy_row + nrows] = \
			                                             self.ground_truths[class_index]
			copy_row += nrows

		self.data_matrix_is_contiguous = True
		return self.data_matrix

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
		if self.interpolation_coefficients != None and len( self.interpolation_coefficients ) != 0:
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

<<<<<<< HEAD
	#==============================================================
	def Split( self, randomize=True, balanced_classes=True, training_set_fraction=None,\
	           i=None, j=None, training_set_only=False, quiet=False ):
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

		if i and( i <= 0 or i > smallest_class_size ):
			raise ValueError( 'i must be greater than zero and less than total number of images'\
			    + ' in smallest class ({0})'.format( smallest_class_size ) )

		if j and( j <= 0 or j > smallest_class_size ):
			raise ValueError( 'j must be greater than zero and less than total number of images'\
			    + ' in smallest class ({0})'.format( smallest_class_size ) )

		if ( i and j ) and ( ( i + j ) > smallest_class_size ):
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

		if j and training_set_only and ( j > 0 ):
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
				num_samples_in_test_set = [ a - b for a, b in zip( self.classsizes_list, num_samples_in_training_set ) ]

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
		if self.feature_vector_version:
			training_set.feature_vector_version = self.feature_vector_version
		else:
			training_set.feature_vector_version = '2.0'
		if self.interpolation_coefficients != None and len( self.interpolation_coefficients ) != 0:
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
			if self.feature_vector_version:
				test_set.feature_vector_version = self.feature_vector_version
			else:
				test_set.feature_vector_version = '3.0'
			if self.interpolation_coefficients != None and len( self.interpolation_coefficients ) != 0:
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

		training_set.ContiguousDataMatrix()

		if training_set_only:
			return training_set

		test_set.ContiguousDataMatrix()
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

=======
>>>>>>> Major refactor of FeatureSet, ArtificialFeatureSets, unittests
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


	#==============================================================
	def __init__( self, data_dict = None ):

		# call parent constructor
		super( FeatureSet_Continuous, self ).__init__()

	#==============================================================
	def ScrambleGroundTruths( self ):
		"""Produce an instant negative control training set"""

		import random
		random.shuffle( self.ground_truths )
		self.source_path += " (scrambled)"

# END FeatureSet_Continuous class definition
=======
=======
>>>>>>> Adding return self statements at the end of methods to allow chaining of multiple function calls
# END FeatureSpace class definition
>>>>>>> API changes. FeatureSet->FeatureSpace. member names more consistent across objects

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

	def __init__( self ):
		self.name = None
		self.source_filepath = None
		self.ground_truth_value = None
		self.predicted_value = None
		self.batch_number = None

		#: For the future: a member to indicate the position of the ROI within an image
		self.tile_index = None

	@output_railroad_switch
	def Print():
		raise NotImplementedError


#=================================================================================
class DiscreteImageClassificationResult( ImageClassificationResult ):
	"""Concrete class that contains the result of a classification for a single image/ROI.
	which includes predicted class, marginal probabilities, etc."""

	def __init__( self ):
		"""Constructor"""
		super( DiscreteImageClassificationResult, self ).__init__()

		self.marginal_probabilities = []
		self.normalization_factor = None
		self.marginal_probabilities = None

		#: predicted_class_name will always be a string
		#: the interpolated value, if applicable, gets stored in self.predicted_vlaue
		self.predicted_class_name = 'UNKNOWN'
		self.ground_truth_class_name = 'UNKNOWN'

	#==============================================================
	def __eq__( self, other ):

		for a, b in zip( self.marginal_probabilities, other.marginal_probabilities ):
			# FIXME: For now, define equal as marg_probs agreeing within 5%
			if abs( a-b ) > 0.005:
				return False
		# The following is pretty useless since normalization factors are often in the 10e-20
		# range.
		#if abs( self.normalization_factor - other.normalization_factor ) > 0.001:
		#	return False
		return True

	#==============================================================
	def __ne__( self, other ):
		return not self == other

	#==============================================================
	@output_railroad_switch
	def Print( self, line_item=False ):
		"""Output classification line item data, including predicted class and marginal
		probabilities"""
		
		if line_item:
			print str(self)
		else:
			print "Image:             \t{0}".format( self.source_filepath )
			print "Normalization Factor:\t{0}".format( self.normalization_factor )
			print "Marg. Probabilities:\t" + "\t".join(\
					[ "{val:0.3f}".format( val=prob ) for prob in self.marginal_probabilities ] )
			print "Ground truth class:\t {0}".format( self.ground_truth_class_name ) 
			print "Predicted class:\t {0}".format( self.predicted_class_name ) 
			if self.predicted_value is not None:
				print "Interpolated Value:\t{0}".format( self.predicted_value )
				#print "Interpolated Value:\t{val=0.3f}".format( val=self.predicted_value )

	#==============================================================
	def __str__( self ):
		# img name:
		output_str = self.source_filepath if self.source_filepath else ""
		if self.tile_index is not None:
			if self.tile_index == 'AVG':
				output_str += " (AVG)"
			else:
				output_str += " ({0}/{1})".format( self.tile_index + 1, self.num_samples_in_group )
		# normalization factor:
		if self.normalization_factor is None:
			# no normalization factor means this is a non-call
			print output_str + "\t--COLLISION--"
			return
		output_str += "\t{0:0.3g}\t".format( self.normalization_factor )

		# marginal probabilities:
		output_str += "\t".join(\
				 [ "{0:0.3f}".format( prob ) for prob in self.marginal_probabilities ] )
		output_str += "\t"
		# actual class:
		if self.ground_truth_class_name:
			output_str += "{0}\t".format( self.ground_truth_class_name )
		else:
			output_str += "*\t"
		# predicted class:
		output_str += self.predicted_class_name + "\t"
		# interpolated value, if applicable
		if self.predicted_value is not None:
			output_str += "{0:0.3f}".format( self.predicted_value )
		return output_str

	#==============================================================
	def __repr__( self ):
		return str( self )

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

			denom = num_tiles - num_collisions
			if denom == 0:
				# This sample collided with every sample in the test set
				# return a non-call
				return cls()
			class_similarities[ class_index ] /= denom
#			print "class_similarities: "+str(class_similarities)

		result = cls()
		norm_factor = sum( class_similarities )
		result.normalization_factor = norm_factor 
		result.marginal_probabilities = [ x / norm_factor for x in class_similarities ]
		return result

	#=================================================================================
	@classmethod
	def NewWND5( cls, training_set, feature_weights, test_samp, quiet = False ):
		"""@brief: A wrapper function for _ClassifyOneImageWND5 that does dummyproofing
		@return: An instance of a DiscreteBatchClassificationResult"""

		if not isinstance( training_set, FeatureSpace ):
			raise ValueError( 'First argument to NewWND5 must be of type "FeatureSpace", you gave a {0}'.format( type( training_set ).__name__ ) )
		
		if not isinstance( feature_weights, FeatureWeights ):
			raise ValueError( 'Second argument to NewWND5 must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )

		if not isinstance( test_samp, FeatureVector ):
			raise ValueError( 'Third argument to NewWND5 must be of type "FeatureVector", you gave a {0}'.format( type( test_samp ).__name__ ) )

		train_set_len = len( training_set.featurenames_list )
		test_set_len = len( test_samp.featurenames_list )
		feature_weights_len = len( feature_weights.featurenames_list )

		if test_samp.featurenames_list != feature_weights.featurenames_list:
			raise ValueError("Can't classify, features in signature don't match features in weights." )

		if test_samp.featurenames_list != training_set.featurenames_list:
			raise ValueError("Can't classify, features in signature don't match features in training_set." )

		if not quiet:
			print "Classifying image '{0}' ({1} features) against test set '{2}' ({3} features)".\
			 format( test_samp.name, train_set_len, training_set.name, test_set_len )

		result = cls._WND5( training_set, test_samp.values, feature_weights.values )

		result.source_filepath = test_samp.source_filepath
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
	def __init__( self ):
		"""Constructor"""
		super( ContinuousImageClassificationResult, self ).__init__()

	#==============================================================
	@output_railroad_switch
	def Print( self, line_item = False ):
		"""Output results."""

		if line_item:
			# img name:
			output_str = str( self.source_filepath )
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
			print "Image:             \t{0}".format( self.source_filepath )
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

	#==============================================================
	@classmethod
	def _LeastSquaresRegression( cls, one_image_features, training_set ):
		"""Produce a predicted value for a single image based on numpy.linalg.lstsq().
		
		Don't call this function directly, but instead use the member function
		NewLeastSquaresRegression() on the ContinuousBatchClassificationResult class
		which has dummyproofing and performs the necessary matrix augmentation."""

		from numpy.linalg import lstsq
		from numpy import dot

		A = lstsq( training_set.data_matrix, np.array( training_set.ground_truths ) )[0]
		
		result = cls()
		result.predicted_value = dot( one_image_features, A )
		return result

#=================================================================================
class BatchClassificationResult( ClassificationResult ):
	"""An abstract base class which serves as container for individual 
	ImageClassificationResult instances.
	
	The equivalent concept to a BatchClassificationResult in the C++ implementation
	of WND-CHARM is the split, i.e., command line arg of -n10 would generate 10 
	train/test splits."""

	#==============================================================
	def __init__( self, training_set=None, test_set=None, feature_weights=None, name=None,
		batch_number=None ):
		"""BatchClassificationResult constructor"""

		self.training_set = training_set
		self.test_set = test_set
		self.feature_weights = feature_weights
		self.name = name
		if self.name is None and training_set and training_set.name:
			self.name = training_set.name
		self.individual_results = []
		self.tiled_results = None
		self.batch_number = batch_number

		self.num_classifications = 0

		self.figure_of_merit = None
		self.predicted_values = None
		self.ground_truth_values = None

		# This stuff is for correllation analysis
		self.pearson_coeff = None
		self.pearson_p_value = None
		self.pearson_std_err = None
		self.spearman_coeff = None
		self.spearman_p_value = None

	#==============================================================	
	def __str__( self ):
		outstr = self.__class__.__name__
		if self.name is None and self.batch_number is None:
			return outstr + " id({})".format( id( self ))
		if self.batch_number is not None:
			outstr += ' batch # {}'.format( self.batch_number )
		if self.name is not None:
			outstr += ' "{}"'.format( self.name )
		return outstr
	#==============================================================
	def __repr__( self ):
		return str(self)

	#==============================================================
	def __len__( self ):
		return len( self.individual_results )

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

			from math import sqrt
			from scipy.stats import linregress, spearmanr

			self.figure_of_merit = sqrt( err_sum / self.num_classifications )

			# no point in doing regression stuff if there's only 1 individual result:
			if len( self.individual_results ) == 1:
				# For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
				np.seterr (under='ignore')
				slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
			             linregress( self.ground_truth_values, self.predicted_values )

				try:
					self.spearman_coeff, self.spearman_p_value =\
			       spearmanr( self.ground_truth_values, self.predicted_values )
				except FloatingPointError:
					# to avoid: "FloatingPointError: invalid value encountered in true_divide"
					self.spearman_coeff, self.spearman_p_value = ( 0, 1 )

				np.seterr (all='raise')

	#==============================================================
	def NewNFold( self, num_folds=5, *args, **kwargs ):
		"""Base method that's called by daughter classes.
		
		If this BatchClassificationResult contains ImageClassificationResults that
		has known ground truth as well as predicted values, this method will calculate
		statistics describing how well predicted values correlate with ground truth.

		Requires scipy.stats package to be installed"""
		raise NotImplementedError


	#==============================================================
	def RankOrderSort( self ):
		"""Rank-order sorts ground truth/predicted value data points for purposes of 
		being graphed."""

		if not self.ground_truth_values or not self.predicted_values:
			self.GenerateStats()

		# Check again:
		if not self.ground_truth_values:
			raise ValueError( "Can't rank-order sort: no ground truth sample values" )
		if not self.predicted_values:
			# FIXME: this might be wrong, since the member predicted_values may contain a
			# non-graphable label string
			raise ValueError( "Can't rank-order sort: no sample predicted values" )

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

	#: batch_count class attribute
	batch_count = 0

	#==============================================================
	def __init__( self, *args, **kwargs ):
		"""Possible kwargs, with defaults:
		training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

		super( DiscreteBatchClassificationResult, self ).__init__( *args, **kwargs )

		self.num_correct_classifications = None
		self.classification_accuracy = None

		self.num_classifications_per_class = None
		self.num_correct_classifications_per_class = None

		self.confusion_matrix = None
		self.similarity_matrix = None
		self.average_class_probability_matrix = None

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
		from collections import defaultdict # introduced Python 2.5

		self.confusion_matrix = {}
		self.average_class_probability_matrix = {}

		# initialize the rows
		# Does the test set have ground truth?
		if self.test_set.classnames_list:
			for test_set_class_name in self.test_set.classnames_list:
				self.confusion_matrix[ test_set_class_name ] = defaultdict(int) 
				self.average_class_probability_matrix[ test_set_class_name ] = defaultdict(float)
		else:
			self.confusion_matrix[ "UNKNOWN" ] = defaultdict( int )
			self.average_class_probability_matrix[ "UNKNOWN" ] = defaultdict( float )

		self.num_correct_classifications = 0

		self.num_classifications_per_class = defaultdict( int )
		self.num_correct_classifications_per_class = defaultdict( int )

		if self.tiled_results:
			classification_results = self.tiled_results
		else:
			classification_results = self.individual_results

		for indiv_result in classification_results:
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
				if indiv_result.marginal_probabilities:
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
					if self.similarity_matrix[ row ][ row ]:
						self.similarity_matrix[ row ][ col ] /= self.similarity_matrix[ row ][ row ]
						self.similarity_matrix[ col ][ row ] /= self.similarity_matrix[ row ][ row ]

		self.classification_accuracy = float( self.num_correct_classifications) / float( self.num_classifications )
		return self

	#==============================================================
	@output_railroad_switch
	def Print( self ):
		"""Prints out statistics from this batch's image classifications, which include 
		classification accuracy, confusion matrix, similarity matrix, and average class 
		probability matrix."""

		if self.classification_accuracy == None:
			self.GenerateStats()

		if self.tiled_results:
			classification_results = self.tiled_results
		else:
			classification_results = self.individual_results

		print "==========================================="
		print "Batch summary:"
		if self.name:
			print "Name: ", self.name
		print "Total number of classifications: {0}".format( len( classification_results ) )
		print "Total number of CORRECT classifications: {0}".format( self.num_correct_classifications )
		print "Total classification accuracy: {0:0.4f}".format( self.classification_accuracy )
		if self.figure_of_merit is not None:
			print "Standard Error: {0:0.4f}".format( self.figure_of_merit )
		if self.pearson_coeff is not None:
			print "Pearson Coefficient: {0:0.4f}".format( self.pearson_coeff )
		if self.spearman_coeff is not None:
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
	def PrintDistanceMatrix( self, method='mps' ):
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
	def New( cls, training_set, test_set, feature_weights, batch_number=None,
					batch_name=None, quiet=False, norm_factor_threshold=None):
		"""The equivalent of the "wndcharm classify" command in the command line implementation
		of WND-CHARM. Input a training set, a test set, and feature weights, and returns a
		new instance of a DiscreteBatchClassificationResult, with self.individual_results
		filled with a new instances of DiscreteImageClassificationResult.

		FIXME: What happens when the ground truth is not known? Currently they would all be shoved
		into class 1, might not be a big deal since class name should be something
		like "UNKNOWN"
		"""

		# type checking
		if not isinstance( training_set, FeatureSpace ):
			raise ValueError( 'First argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
		if not isinstance( test_set, FeatureSpace ):
			raise ValueError( 'Second argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
		if not isinstance( feature_weights, FeatureWeights ):
			raise ValueError( 'Third argument to New must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )
	
		# feature comparison
		if test_set.featurenames_list != feature_weights.featurenames_list:
			raise ValueError( "Can't classify, features in test set don't match features in weights. Try translating feature names from old style to new, or performing a FeatureReduce()" )
		if test_set.featurenames_list != training_set.featurenames_list:
			raise ValueError( "Can't classify, features in test set don't match features in training set. Try translating feature names from old style to new, or performing a FeatureReduce()" )

		np.seterr (under='ignore')

		train_set_len = len( training_set.featurenames_list )
		test_set_len = len( test_set.featurenames_list )
		feature_weights_len = len( feature_weights.featurenames_list )

		# instantiate myself
		batch_result = cls( training_set, test_set, feature_weights, batch_name )

		# Are the samples to be classified tiled?
		if test_set.num_samples_per_group > 1:
			batch_result.tiled_results = []

		# Say what we're going to do
		if not quiet:
			print "Classifying test set '{0}' ({1} features) against training set '{2}' ({3} features)".\
					format( test_set.name, test_set_len, training_set.name, train_set_len )
			if batch_result.tiled_results:
				print "Performing tiled classification."
			column_header = "image\tnorm. fact.\t"
			column_header +=\
				"".join( [ "p(" + class_name + ")\t" for class_name in training_set.classnames_list ] )
			column_header += "act. class\tpred. class\tpred. val."
			print column_header

		# Give myself a number so that it looks good when I print out results
		if not batch_number:
			batch_number = DiscreteBatchClassificationResult.batch_count
			DiscreteBatchClassificationResult.batch_count += 1
		batch_result.batch_number = batch_number

		# Will there be a numeric predicted value associated with this classification?
		train_set_interp_coeffs = None
		if training_set.interpolation_coefficients != None and len( training_set.interpolation_coefficients) != 0:
			train_set_interp_coeffs = np.array( training_set.interpolation_coefficients )
			batch_result.predicted_values = []

		# Are there numeric ground truth values associated with the input samples?
		test_set_interp_coeffs = None
		if test_set.interpolation_coefficients != None and len( test_set.interpolation_coefficients ) != 0:
			test_set_interp_coeffs = np.array( test_set.interpolation_coefficients )
			batch_result.ground_truth_values = []

		for test_class_index in range( test_set.num_classes ):
			num_class_imgs, num_class_features = test_set.data_list[ test_class_index ].shape

			# Get tiling ready if need be
			if test_set.num_samples_per_group > 1:
				tile_results_in_this_sample_group = []

			for test_image_index in range( num_class_imgs ):
				one_image_features = test_set.data_list[ test_class_index ][ test_image_index,: ]
				result = DiscreteImageClassificationResult._WND5( training_set, one_image_features, feature_weights.values )
				
				if norm_factor_threshold and (result.normalization_factor > norm_factor_threshold):
					continue
				result.source_filepath = test_set.samplenames_list[ test_class_index ][ test_image_index ]
				result.ground_truth_class_name = test_set.classnames_list[ test_class_index ]
				result.batch_number = batch_number
				result.name = batch_name
				if result.marginal_probabilities:
					# Sometimes the result comes back with a non-call, like when the sample image
					# collides with every test image
					marg_probs = np.array( result.marginal_probabilities )
					result.predicted_class_name = training_set.classnames_list[ marg_probs.argmax() ]
				else:
					marg_probs = np.array( [np.nan] * training_set.num_classes )
					result.predicted_class_name = "*Collided with every training sample*"

				# interpolated value, if applicable
				if train_set_interp_coeffs is not None:
					interp_val = np.sum( marg_probs * train_set_interp_coeffs )
					result.predicted_value = interp_val
					batch_result.predicted_values.append( interp_val )

				if test_set_interp_coeffs is not None:
					result.ground_truth_value = test_set_interp_coeffs[ test_class_index ]
					batch_result.ground_truth_values.append( test_set_interp_coeffs[ test_class_index ] )

				if test_set.num_samples_per_group > 1:
					result.tile_index = len( tile_results_in_this_sample_group )
					result.num_samples_in_group = test_set.num_samples_per_group

				if not quiet:
					result.Print( line_item=True )
				batch_result.individual_results.append( result )

				# TILING SECTION:
				# Create a whole image classification result that
				# is the average of all the calls from all the tiles
				if test_set.num_samples_per_group > 1:
					tile_results_in_this_sample_group.append( result )
					if len( tile_results_in_this_sample_group ) >= test_set.num_samples_per_group:
						aggregated_result = DiscreteImageClassificationResult()
						# Use the last result from above
						aggregated_result.source_filepath = result.source_filepath
						aggregated_result.tile_index = 'AVG'
						aggregated_result.ground_truth_class_name = result.ground_truth_class_name
						marg_prob_lists = [ [] for i in xrange( training_set.num_classes ) ]
						norm_factor_list = []
						for tile_result in tile_results_in_this_sample_group:

							if tile_result.marginal_probabilities:
								# Sometimes the result comes back with a non-call, like when the sample image
								# collides with every test image
								for class_index, val in enumerate( tile_result.marginal_probabilities ):
									marg_prob_lists[ class_index ].append( val )
								norm_factor_list.append( tile_result.normalization_factor )

						if any( [ len( class_marg_prob_list ) <= 0 for class_marg_prob_list in marg_prob_lists ] ):
							aggregated_result.marginal_probabilities = np.array( [np.nan] * training_set.num_classes )
							aggregated_result.predicted_class_name = "*Collided with every training image*"
						else:
							marg_probs = [ sum( l ) / float( len( l ) ) for l in marg_prob_lists ]
							aggregated_result.marginal_probabilities = marg_probs
							np_marg_probs = np.array( marg_probs )
							aggregated_result.predicted_class_name = training_set.classnames_list[ np_marg_probs.argmax() ]
							aggregated_result.normalization_factor = sum(norm_factor_list)/float(len(norm_factor_list))

							# interpolated value, if applicable
							if train_set_interp_coeffs is not None:
								aggregated_result.predicted_value = np.sum( np_marg_probs * train_set_interp_coeffs )
							if test_set_interp_coeffs is not None:
								aggregated_result.ground_truth_value = test_set_interp_coeffs[ test_class_index ]
						# Save references to the tiled results on the aggregated result
						# To facilitate creation of graphs, heatmaps, etc
						aggregated_result.tiled_results = tile_results_in_this_sample_group
						batch_result.tiled_results.append( aggregated_result )
						if not quiet:
							aggregated_result.Print( line_item=True )

						# now reset for the next sample group
						tile_results_in_this_sample_group = []
						tile_count = 0

		np.seterr (all='raise')

		return batch_result

#=================================================================================
class ContinuousBatchClassificationResult( BatchClassificationResult ):
	"""Concrete class that is a container for ContinuousImageClassificationResult instances.
	Use this class to run classification experiments and to generate statistics on their results.
	For continuous (Pearson) classifiers, the figure_of_merit is the standard error between
	predicted and ground truth."""

	#: batch_count class attribute
	batch_count = 0

	def __init__( self, *args, **kwargs ):
		"""Possible kwargs, with defaults:
		training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

		super( ContinuousBatchClassificationResult, self ).__init__( *args, **kwargs )
		self.predicted_values = []

	#=====================================================================
	def GenerateStats( self ):
		"""Calls base class method to run a standard error analysis if ground_truth/
		predicted vals exist."""
		super( ContinuousBatchClassificationResult, self ).GenerateStats()
		return self

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
	def New( cls, test_set, feature_weights, quiet=False, batch_number=None, batch_name=None):
		"""Uses Pearson-coefficient weighted Regression-Voting classifier."""

		# type checking
		if not isinstance( test_set, FeatureSpace ):
			raise ValueError( 'First argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )	
		if not isinstance( feature_weights, ContinuousFeatureWeights ):
			raise ValueError( 'Second argument to New must be of type "ContinuousFeatureWeights", you gave a {0}'.format( type( feature_weights ).__name__ ) )

		# feature comparison
		if test_set.featurenames_list != feature_weights.featurenames_list:
			raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )

		# say what we're gonna do
		if not quiet:
			out_str = 'Classifying test set "{0}" ({1} images, {2} features)\n\tagainst training set "{3}" ({4} images)'
			print out_str.format( test_set.name, \
			                      test_set.num_samples, \
			                      len( test_set.featurenames_list ), \
			                      feature_weights.associated_training_set.name, \
			                      feature_weights.associated_training_set.num_samples )

		if not quiet:
			column_header = "image\tground truth\tpred. val."
			print column_header

		batch_result = cls( feature_weights.associated_training_set, test_set, feature_weights )
		batch_result.name = batch_name

		if test_set.ground_truths is not None and len( test_set.ground_truths ) != 0:
			batch_result.ground_truth_values = test_set.ground_truths

		for test_image_index in range( test_set.num_samples ):
			one_image_features = test_set.data_matrix[ test_image_index,: ]
			result = ContinuousImageClassificationResult._LinearRegression( one_image_features, feature_weights )
			result.batch_number = batch_number
			result.name = batch_name
			result.source_filepath = test_set.samplenames_list[ test_image_index ]
			result.ground_truth_value = test_set.ground_truths[ test_image_index ]
			batch_result.predicted_values.append( result.predicted_value )

			if not quiet:
				result.Print( line_item = True )
			batch_result.individual_results.append( result )

		batch_result.GenerateStats()
		return batch_result

	#=====================================================================
	@classmethod
	def NewLeastSquaresRegression( cls, training_set, test_set, feature_weights,
					leave_one_out=False, quiet=False, batch_number=None, batch_name=None):
		"""Uses Linear Least Squares Regression classifier in a feature space filtered/weighed
		by Pearson coefficients.
		
		Use cases:
		1. if training_set != test_set and both not none: straight classification
		2. if training_set is not None and test_set is None = Leave one out cross validation
		3. if training_set == test_set == not None: Shuffle/Split cross validation
		"""

		# Type checking
		# First arg must be a FeatureSpace
		if not isinstance( training_set, FeatureSpace ):
			raise ValueError( 'First argument to NewLeastSquaresRegression must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
		# Second arg must be a FeatureSpace or a None
		if (not isinstance( test_set, FeatureSpace )) and (test_set is not None):
			raise ValueError( 'Second argument to NewLeastSquaresRegression must be of type "FeatureSpace" or None, you gave a {0}'.format( type( test_set ).__name__ ) )
		if not isinstance( feature_weights, ContinuousFeatureWeights ):
			raise ValueError( 'Third argument to New must be of type "ContinuousFeatureWeights", you gave a {0}'.format( type( feature_weights ).__name__ ) )

		# If there's both a training_set and a test_set, they both have to have the same features
		if training_set and test_set:
			if training_set.featurenames_list != training_set.featurenames_list:
				raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )
		# Check feature_weights
		if training_set.featurenames_list != feature_weights.featurenames_list:
			raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )

		# Least squares regression requires a featres matrix augmented with ones
		# signifying the constant, i.e., the y-intercept
		# Make a copy of the feature sets to augment while leaving original matrices unchanged.
		from copy import deepcopy
		augmented_train_set = deepcopy( training_set )

		# figure out what we're gonna do
		if training_set and not test_set:
			leave_one_out = True
			cross_validation = True
			augmented_test_set = augmented_train_set 
		elif training_set is test_set:
			cross_validation = True
			augmented_test_set = augmented_train_set
		else:
			augmented_test_set = deepcopy( test_set )
			cross_validation = False

		# say what we're gonna do
		if not quiet:
			if cross_validation:
				out_str = 'Cross validation of training set "{0}" ({1} images, {2} features)'.format(
			          training_set.name, training_set.num_samples, len( training_set.featurenames_list ) )

				out_str += "\nWITH"
				if not leave_one_out:
					out_str += "OUT"
				out_str += " using leave-one-out analysis.\n"
			
			else:
				out_str = 'Classifying test set "{0}" ({1} images, {2} features)'.format(
			          test_set.name, test_set.num_samples, len( test_set.featurenames_list ) )
				out_str += '\n\tagainst training set "{0}" ({1} images)'.format(
				            training_set.name, training_set.num_samples )
			print out_str

		# Now, build the augmented feature matrices, which includes multiplying the feature
		# space by the weights, and augmenting the matrices with 1's

		oldsettings = np.seterr(all='ignore')

		augmented_train_set.data_matrix *= feature_weights.values
		augmented_train_set.data_matrix = np.hstack( 
		      [ augmented_train_set.data_matrix, np.ones( ( augmented_train_set.num_samples, 1 ) ) ] )
		augmented_train_set.num_features += 1 # Tell the object it has a new feature column
		if not cross_validation:
			augmented_test_set.data_matrix *= feature_weights.values
			augmented_test_set.data_matrix = np.hstack( 
		      [ augmented_test_set.data_matrix, np.ones( ( augmented_test_set.num_samples, 1 ) ) ] )
			augmented_test_set.num_features += 1 # Tell the object it has a new feature column

		if not quiet:
			print "image\tground truth\tpred. val."

		batch_result = cls( training_set, test_set, feature_weights )
		batch_result.name = batch_name

		if augmented_test_set.ground_truths is not None and len( augmented_test_set.ground_truths ) != 0:
			batch_result.ground_truth_values = augmented_test_set.ground_truths

		intermediate_train_set = augmented_train_set
		for test_image_index in range( augmented_test_set.num_samples ):
			if leave_one_out:
				sample_group_id = augmented_train_set.samplegroupid_list[ test_image_index ]
				intermediate_train_set = augmented_train_set.SampleReduce(\
						leave_out_samplegroupid_list = sample_group_id,)
			one_image_features = augmented_test_set.data_matrix[ test_image_index,: ]
			result = ContinuousImageClassificationResult._LeastSquaresRegression( \
			               one_image_features, intermediate_train_set )
			result.batch_number = batch_number
			result.name = batch_name
			result.source_filepath = augmented_test_set.samplenames_list[ test_image_index ]
			result.ground_truth_value = augmented_test_set.ground_truths[ test_image_index ]
			batch_result.predicted_values.append( result.predicted_value )

			if not quiet:
				result.Print( line_item = True )
			batch_result.individual_results.append( result )

		# return settings to original
		np.seterr(**oldsettings)

		batch_result.GenerateStats()
		return batch_result

#============================================================================
class ClassificationExperimentResult( BatchClassificationResult ):
	"""Abstract base class which serves as a container for BatchClassificationResult instances
	and their associated statistics.
	
	N.B. This is a container for batches, not for individual results, i.e.,
	the list self.individual_results contains instances of type BatchClassificationResult."""

	def __init__( self, *args, **kwargs ):
		"""Possible kwargs, with defaults:
		training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

		super( ClassificationExperimentResult, self ).__init__( *args, **kwargs )

		#: A dictionary where the name is the key, and the value is a list of individual results
		self.accumulated_individual_results = None

		self.feature_weight_statistics = None

		#: keep track of stats related to predicted values for reporting purposes
		self.individual_stats = None

	#=====================================================================
	def GenerateStats( self ):
		"""Only partially implemented.
		
		Completed functionality:
		1. Simple aggregation of ground truth->predicted value pairs across splits.
		   Use the function PerSampleStatistics() to average results for specific
			 images across splits.
		2. Calculation of aggregated feature weight statistics

		Considerations for future implementation.
		1. The test set may or may not have ground truth (discrete and continuous)
		2. The results may not have a predicted value (discrete only)
		3. Continuous classifications do not have marginal probabilities
		4. Hybrid test sets (discrete test sets loaded into a continuous test set)
		   have "pseudo-classes," i.e., easily binnable ground truth values."""

		# Simple result aggregation.
		# Calls GenerateStats() on the individual batches if the 
		# ground truth->predicted value pairs haven't been scraped
		# from the batch's list of individual ImageClassificationResults.
		from itertools import chain

		lists_of_ground_truths = []
		lists_of_predicted_values = []

		self.figure_of_merit = 0
		for batch_result in self.individual_results:
			self.num_classifications += len( batch_result.individual_results )
			if batch_result.figure_of_merit == None:
				batch_result.GenerateStats()
			if batch_result.ground_truth_values:
				lists_of_ground_truths.append( batch_result.ground_truth_values )
			if batch_result.predicted_values:
				lists_of_predicted_values.append( batch_result.predicted_values )

		if lists_of_ground_truths:
			self.ground_truth_values = list( chain( *lists_of_ground_truths ) )
		if lists_of_predicted_values:
			self.predicted_values = list( chain( *lists_of_predicted_values ) )

		# Aggregate feature weight statistics across splits, if any:
		feature_weight_lists = {}
		for batch_result in self.individual_results:
			if not batch_result.feature_weights:
				continue
			weight_names_and_values = zip( batch_result.feature_weights.featurenames_list, 
					                                        batch_result.feature_weights.values)
			for name, weight in weight_names_and_values:
				if not name in feature_weight_lists:
					feature_weight_lists[ name ] = []
				feature_weight_lists[ name ].append( weight )

		if feature_weight_lists is not {}:
			for feature_name in feature_weight_lists:
				feature_weight_lists[ feature_name ] = np.array( feature_weight_lists[ feature_name ] )

			feature_weight_stats = []
			for fname in feature_weight_lists:
				fwl = feature_weight_lists[ fname ]
				count = len( fwl )
				fwl_w_zeros = np.zeros( len( self.individual_results ) )
				fwl_w_zeros[0:count] = fwl
				feature_weight_stats.append( ( np.mean( fwl_w_zeros ),
				                             count, np.std(fwl), np.min(fwl), np.max(fwl), fname ) )

			# Sort on mean values, i.e. index 0 of tuple created above
			self.feature_weight_statistics = sorted( feature_weight_stats, key=lambda a: a[0], reverse=True )

		# Remember, there's no such thing as a confusion matrix for a continuous class
		return self

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
				format( self.training_set.source_filepath )

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

	#=====================================================================
	@classmethod
	def NewShuffleSplit( cls, feature_set, n_iter=5, name=None, features_size=0.15,
		                   train_size=None, test_size=None, random_state=True, classifier=None,
		                   quiet=False):
		"""args train_size, test_size, and random_state are all passed through to Split()
		feature_size if a float is feature usage fraction, if in is top n features."""

		experiment = cls( training_set=feature_set, test_set=feature_set, name=name )
		if isinstance( features_size, float ):
			if features_size < 0 or features_size > 1.0:
				raise ValueError('Arg "features_size" must be on interval [0,1] if a float.')
			num_features = int( round( features_size * feature_set.num_features ) )
		elif isinstance( features_size, int ):
			if features_size < 0 or features_size > feature_set.num_features:
				raise ValueError( 'must specify num_features or feature_usage_fraction in kwargs')
			num_features = features_size
		else:
			raise ValueError( 'Arg "features_size" must be valid float or int.' )

		if not quiet:
			print "using top "+str (num_features)+" features"

		# If you passed the same random_state into Split, you'd get the same exact split for
		# all n_iter. Therefore use the seed passed in here to predictably generate seeds
		# for the Split() iterations.
		if random_state:
			from numpy.random import RandomState
			from functools import partial
			maxint = 2 ** 32 - 1
			if random_state is True:
				from numpy.random import randint as np_randint
				randint = partial( np_randint, low=0, high=maxint )
			elif type( random_state ) is int:
				randint = partial( RandomState( random_state ).randint, low=0, high=maxint )
			elif type( random_state ) is RandomState:
				randint = partial( random_state.randint, low=0, high=maxint )
			else:
				raise ValueError( "arg random_state must be an int, instance of numpy.random.RandomState, or True")
		else:
			randint = lambda: None # samples split the same way all iterations - not recommended!

		for split_index in range( n_iter ):
			train_set, test_set = feature_set.Split(
			    train_size, test_size, random_state=randint(), quiet=quiet )
			train_set.Normalize( quiet=quiet )
			
			if feature_set.discrete:
				weights = \
				  FisherFeatureWeights.NewFromFeatureSpace( train_set ).Threshold( num_features )
			else:	
				weights = \
				  ContinuousFeatureWeights.NewFromFeatureSpace( train_set ).Threshold( num_features )

			weights.Print()
			reduced_train_set = train_set.FeatureReduce( weights )
			reduced_test_set = test_set.FeatureReduce( weights )
			reduced_test_set.Normalize( reduced_train_set, quiet=quiet )

			if feature_set.discrete:
				batch_result = DiscreteBatchClassificationResult.New( reduced_train_set, \
		         reduced_test_set, weights, batch_number=split_index, quiet=quiet )
			else:
				if classifier == 'voting':
					batch_result = ContinuousBatchClassificationResult.New(
							reduced_train_set, weights, batch_number=split_index, quiet=quiet )
				else: # default classifier == 'lstsq':
					batch_result = ContinuousBatchClassificationResult.NewLeastSquaresRegression(
						reduced_train_set, reduced_test_set, weights, batch_number=split_index, quiet=quiet )

			experiment.individual_results.append( batch_result )

		return experiment

#============================================================================
class DiscreteClassificationExperimentResult( ClassificationExperimentResult ):
	"""Concrete class which serves as a container for DiscreteBatchClassificationResults
	and their associated statistics. The information contained here comprises everything
	that would appear in a HTML file generated by the C++ implementation of WND-CHARM.

	In this subclass, the figure of merit is self.classification_accuracy"""

	def __init__( self, *args, **kwargs ):
		"""Possible kwargs, with defaults:
		training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

		super( DiscreteClassificationExperimentResult, self ).__init__( *args, **kwargs )

		self.num_correct_classifications = None

		self.confusion_matrix = None
		self.average_similarity_matrix = None
		self.average_class_probability_matrix = None

	#=====================================================================
	def GenerateStats( self ):
		"""Not fully implemented yet. Need to implement generation of confusion, similarity, and
		average class probability matrices from constituent batches."""

		# Base class does feature weight analysis, ground truth-pred. value aggregation
		super( DiscreteClassificationExperimentResult, self ).GenerateStats()
		
		self.num_classifications = 0
		self.num_correct_classifications = 0
		for batch_result in self.individual_results:
			for indiv_result in batch_result.individual_results:
				self.num_classifications += 1
				if indiv_result.ground_truth_class_name == indiv_result.predicted_class_name:
					self.num_correct_classifications += 1

		self.figure_of_merit = float( self.num_correct_classifications) / float( self.num_classifications )
		self.classification_accuracy = self.figure_of_merit

		#FIXME: Create confusion, similarity, and class probability matrices
		return self

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
		print "Total classification accuracy: {0:0.4f}\n".format( self.classification_accuracy )

		print "Batch Accuracies:"
		print "#\tAcc.\tName"
		print "------------------------------------"

		for batch_result in sorted( self.individual_results, key=lambda X: X.classification_accuracy, reverse=True ):
			if batch_result.figure_of_merit == None:
				batch_result.GenerateStats()
			print "{0}\t{1:0.4f}\t{2}".format( batch_result.batch_number,
			                  batch_result.classification_accuracy, batch_result.name )

		outstr = "{0}\t{1:0.3f}\t{2:>3}\t{3:0.3f}\t{4:0.3f}\t{5:0.3f}\t{6}"
		print "Feature Weight Analysis (top 50 features):"
		print "Rank\tmean\tcount\tStdDev\tMin\tMax\tName"
		print "----\t----\t-----\t------\t---\t---\t----"
		for count, fw_stat in enumerate( self.feature_weight_statistics, 1 ):
			print outstr.format( count, *fw_stat )
			if count >= 50:
					break

	#=====================================================================
	@classmethod
	@output_railroad_switch
	def NewFromHTMLReport( cls, path_to_html ):
		"""Takes as input an HTML report generated by C++ implementation wndchrm
		Parses the report and builds up a Pychrm representation of the results,
		which facilitates analysis, graphing, etc."""

		import re
		row_re = re.compile( r'<tr>(.+?)</tr>' )
		name_re = re.compile( r'"(.+?)"' )
		num_re = re.compile( r'(\d*\.?\d+)' )

		# FIXME: This should fail if there isn't some part of the class names that are interpretable
		# as a number, specifically when it tries to calculate an "interpolated" (predicted) value
		# for the sample based on marginal probabilities.

		def ParseClassSummaryHTML( the_html ):
			rows = row_re.findall( the_html )
			ts = FeatureSpace()
			ts.num_classes = 0
			ts.interpolation_coefficients = []
			ts.classnames_list = []
			for rownum, row in enumerate( rows ):
				if rownum == 0:
					continue # skip column header
				ts.num_classes += 1
				classname = re.search( r'<th>(.+?)</th>', row ).group(1)
				ts.classnames_list.append( classname )
			ts.interpolation_coefficients = CheckIfClassNamesAreInterpolatable( ts.classnames_list )
			return ts

		# The following will be determined once the number of classes has been ascertained
		normalization_col = 1
		mp_col = 2
		ground_truth_col = None
		predicted_col = None
		interp_val_col = None # We don't even use this
		name_col = None

		_training_set = None
		_test_set = None
		exp = cls()
		exp.name = path_to_html

		trainingset_definition = False
		trainingset_html = ""
		testset_definition = False
		testset_html = ""

		insidesplit = False
		split = None
		splitcount = 0
		split_linecount = None
		with open( path_to_html ) as file:
			for line in file:
				if 'trainset_summary' in line:
					trainingset_definition = True
				elif trainingset_definition == True:
					trainingset_html += line.strip()
					if '</table>' in line:
						trainingset_definition = False
						ts = _training_set = ParseClassSummaryHTML( trainingset_html )
						ground_truth_col = ts.num_classes + 3
						predicted_col = ts.num_classes + 4
						interp_val_col = ts.num_classes + 6
						name_col = ts.num_classes + 7

				elif 'testset_summary' in line:
					testset_definition = True
				elif testset_definition == True:
					testset_html += line.strip()
					if '</table>' in line:
						testset_definition = False
						_test_set = ParseClassSummaryHTML( testset_html )

				elif line.startswith( '<TABLE ID="IndividualImages_split' ):
					# If we haven't seen a test set definition by now, we ain't gonna see one period.
					if not _test_set:
						_test_set = _training_set
					insidesplit = True
					split = DiscreteBatchClassificationResult( 
					                                   training_set=_training_set, test_set=_test_set )
					split.predicted_values = []
					split.ground_truth_values = []
					splitcount += 1
					split_linecount = 0
				elif insidesplit and line.startswith( '</table><br><br>' ):
					insidesplit = False
					exp.individual_results.append( split )
				elif insidesplit:
					split_linecount += 1
					if split_linecount == 1:
						# First line in individual results is column headers
						# Pure clasification without interpolation
						if 'Interpolated Value' not in line:
							interp_val_col = None
							name_col = ts.num_classes + 6 # one less than set above -- that column won't exist
						continue
					noends = line.strip( '<trd/>\n' ) # take the tr and td tags off front end
					values = noends.split( '</td><td>' )
					result = DiscreteImageClassificationResult()

					result.normalization_factor = float( values[ normalization_col ] )
					result.marginal_probabilities = \
							[ float( val.strip( '</b>' ) ) for val in values[ mp_col : mp_col + _training_set.num_classes ] ]
					result.predicted_class_name = values[ predicted_col ]
					# Sometimes c-chrm labels classes with a * to say it's not part of the training set
					result.ground_truth_class_name = values[ ground_truth_col ].strip('*')
					result.name = name_re.search( values[ name_col ] ).groups()[0]
					result.source_filepath = result.name
					if ts.interpolation_coefficients is not None:
						result.ground_truth_value = \
					    ts.interpolation_coefficients[ ts.classnames_list.index(result.ground_truth_class_name ) ]
						#result.predicted_value = float( values[ interp_val_col ] )
						result.predicted_value = \
					    sum( [ x*y for x,y in zip( result.marginal_probabilities, _training_set.interpolation_coefficients ) ] )
						split.predicted_values.append( result.predicted_value )
						split.ground_truth_values.append( result.ground_truth_value )
					#result.Print( line_item = True )
					result.batch_number = splitcount
					split.individual_results.append(result)

		exp.training_set = _training_set
		exp.test_set = _test_set

		exp.GenerateStats()
		return exp

	#=====================================================================
	@output_railroad_switch
	def PerSampleStatistics( self ):
		"""This function is meant to elucidate the amount of variability of classifications 
		across batches. ImageClassificationResult imformation is aggregated for each individual
		image/ROI encountered, and statistics are generated for each image/ROI and printed out."""

		if self.individual_results == 0:
			raise ValueError( 'No batch results to analyze' )

		#self.predicted_values = []
		self.ground_truth_values = []

		self.accumulated_individual_results = {}
		self.individual_stats = {}

		for batch in self.individual_results:
			for result in batch.individual_results:
				if not result.source_filepath in self.accumulated_individual_results:
					# initialize list of individual results for this file
					self.accumulated_individual_results[ result.source_filepath ] = []
				self.accumulated_individual_results[ result.source_filepath ].append( result )

		for filename in self.accumulated_individual_results:

			# Get marginal probability averages
			mp_totals = None
			for result in self.accumulated_individual_results[filename]:
				if not mp_totals:
					mp_totals = result.marginal_probabilities[:]
				else:
					new_total = []
					for class_total, new_mp in zip( mp_totals, result.marginal_probabilities ):
						new_total.append( class_total + new_mp )
					mp_totals = new_total

			mp_avgs = [ float(mp_totals[i]) / len( self.accumulated_individual_results[filename] ) for i in range( len( mp_totals ) ) ]
			#vals = np.array ([result.predicted_value for result in self.accumulated_individual_results[filename] ])
			vals = [result.predicted_class_name for result in self.accumulated_individual_results[filename] ]
			#self.ground_truth_values.append( self.accumulated_individual_results[filename][0].ground_truth_value )
			gt_class = self.accumulated_individual_results[filename][0].ground_truth_class_name
			self.ground_truth_values.append( gt_class )
			#self.predicted_values.append( np.mean(vals) )
			self.individual_stats[filename] = ( len(vals), float( vals.count( gt_class ) ) / len(vals), mp_avgs, gt_class )

		print "==========================================="
		print 'Experiment name: "{0}"'.format( self.name ) + ' Individual results\n'

		mp_delim = "  "
		discrlineoutstr = "\tsplit {split_num:02d}: pred: {pred_class}\tact: {actual_class}\tnorm factor: {norm_factor:0.3g},\tmarg probs: ( {norm_dists} )"
		outstr = "\t---> Tested {0} times, avg correct: {1:0.3f}, avg marg probs ( {2} )"

		#create view
		res_dict = self.accumulated_individual_results

		# sort by ground truth, then alphanum
		sort_func = lambda A, B: cmp( A, B ) if res_dict[A][0].ground_truth_class_name == res_dict[B][0].ground_truth_class_name else cmp( res_dict[A][0].source_filepath, res_dict[B][0].source_filepath  ) 
		sorted_images = sorted( self.accumulated_individual_results.iterkeys(), sort_func )

		for samplename in sorted_images:
			print 'File "' + samplename + '"'
			for result in self.accumulated_individual_results[ samplename ]:
				marg_probs = [ "{0:0.3f}".format( num ) for num in result.marginal_probabilities ]
				print discrlineoutstr.format( split_num = result.batch_number, \
				                         pred_class = result.predicted_class_name, \
				                         actual_class = result.ground_truth_class_name, \
				                         norm_factor = result.normalization_factor, \
				                         norm_dists = mp_delim.join( marg_probs ) )

			marg_probs = [ "{0:0.3f}".format( num ) for num in self.individual_stats[ samplename ][2] ]
			print outstr.format( self.individual_stats[ samplename ][0], self.individual_stats[ samplename ][1], mp_delim.join( marg_probs ) )


		# If 2 or 3 class problem, plot individuals in marginal probability space



# END class definition for DiscreteClassificationExperimentResult

#============================================================================
class ContinuousClassificationExperimentResult( ClassificationExperimentResult ):
	"""Concrete class which serves as a container for ContinuousBatchClassificationResult instances
	and their associated statistics. The information contained here comprises everything
	that would appear in a HTML file generated by the C++ implementation of WND-CHARM.

	In this subclass, the figure of merit is the average standard error arcoss batches."""

	def __init__( self, *args, **kwargs):
		super( ContinuousClassificationExperimentResult, self ).__init__( *args, **kwargs )

	#=====================================================================
	def GenerateStats( self ):
		"""Calculates statistics describing how well predicted values
		correlate with ground truth across all batches.

		Requires scipy.stats package to be installed"""

		# Base class does feature weight analysis, ground truth-pred. value aggregation
		super( ContinuousClassificationExperimentResult, self ).GenerateStats()
	
		gt = np.array( self.ground_truth_values )
		pv = np.array( self.predicted_values )
		
		diffs = gt - pv
		diffs = np.square( diffs )
		err_sum = np.sum( diffs )

		from math import sqrt
		from scipy.stats import linregress, spearmanr

		self.figure_of_merit = sqrt( err_sum / self.num_classifications )

		# For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
		np.seterr (under='ignore')
		slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
								 linregress( self.ground_truth_values, self.predicted_values )

		try:
			self.spearman_coeff, self.spearman_p_value =\
					 spearmanr( self.ground_truth_values, self.predicted_values )
		except FloatingPointError:
			self.spearman_coeff, self.spearman_p_value = ( 0, 1 )

		np.seterr (all='raise')
		return self

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

		outstr = "{0}\t{1:0.3f}\t{2:>3}\t{3:0.3f}\t{4:0.3f}\t{5:0.3f}\t{6}"
		print "Feature Weight Analysis (top 50 features):"
		print "Rank\tmean\tcount\tStdDev\tMin\tMax\tName"
		print "----\t----\t-----\t------\t---\t---\t----"
		for count, fw_stat in enumerate( self.feature_weight_statistics, 1 ):
			print outstr.format( count, *fw_stat )
			if count >= 50:
					break

	#=====================================================================
	@output_railroad_switch
	def PerSampleStatistics( self ):
		"""This function is meant to elucidate the amount of variability of classifications
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
				if not result.source_filepath in self.accumulated_individual_results:
					# initialize list of individual results for this file
					self.accumulated_individual_results[ result.source_filepath ] = []
				self.accumulated_individual_results[ result.source_filepath ].append( result )

		for filename in self.accumulated_individual_results:
			vals = np.array( [result.predicted_value for result in self.accumulated_individual_results[filename] ])
			self.ground_truth_values.append( self.accumulated_individual_results[filename][0].ground_truth_value )
			self.predicted_values.append( np.mean(vals) )
			self.individual_stats[filename] = ( len(vals), np.min(vals), np.mean(vals), \
																							np.max(vals), np.std(vals) ) 

		print "==========================================="
		print 'Experiment name: "{0}"'.format( self.name ) + ' Individual results\n'
		mp = "  "
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
				print contlineoutstr.format( split_num = result.batch_number, \
				                         batch_name = result.name, \
				                         actual_class = result.ground_truth_value, \
				                         pred_val = result.predicted_value )
			print outstr.format( *self.individual_stats[ samplename ] )

#============================================================================
class BaseGraph( object ):
	"""An abstract base class that is supposed to hold onto objects on which to call
	matplotlib.pyplot API methods."""

	def __init__( self ):

		# general stuff:
		self.chart_title = None
		self.file_name = None
		self.batch_result = None

		# pyplot-specific stuff
		self.figure = None
		self.main_axes = None

	def SaveToFile( self, filepath ):
	
		if self.figure == None:
			raise ValueError( 'No figure to save!' )
		self.figure.savefig( filepath )
		print 'Wrote chart "{0}" to file "{1}"'.format( self.chart_title, filepath )
			
#============================================================================
class PredictedValuesGraph( BaseGraph ):
	"""This is a concrete class that can produce two types of graphs that are produced
	from ImageClassificationResult data stored in a BatchClassificationResult."""

	#=================================================================
	def __init__( self, result, name=None ):
		"""Constructor sorts ground truth values contained in BatchClassificationResult
		and loads them into self.grouped_coords"""

		self.batch_result = result
		self.grouped_coords = {}

		# Right now these are only dealing
		self.classnames_list = result.test_set.classnames_list
		self.class_values = result.test_set.interpolation_coefficients
		self.num_classes = result.test_set.num_classes

		if name:
			self.chart_title = name
		elif result and result.name:
			self.chart_title = result.name
		else:
			self.chart_title = None

		#FIXME: implement user-definable bin edges

		result.RankOrderSort()
		whole_list = zip( result.ground_truth_values, result.predicted_values )

		for class_index, class_name in enumerate( self.classnames_list ):
			self.grouped_coords[ class_name ] = []
			for coords in whole_list:
				if coords[0] == self.class_values[ class_index ]:
					self.grouped_coords[ class_name ].append( coords )

	#=====================================================================
	@classmethod
	@output_railroad_switch
	def NewFromHTMLFile( cls, filepath ):
		"""Helper function to facilitate the fast generation of graphs from C++-generated
		HTML Report files."""

		exp = DiscreteClassificationExperimentResult.NewFromHTMLReport( filepath )
		exp.GenerateStats()
		exp.PredictedValueAnalysis()
		newgraphobj = cls( exp )
		return newgraphobj

	#=====================================================================
	def SaveToFile( self, filepath ):
		"""Calls base class function"""
		super( PredictedValuesGraph, self ).SaveToFile( filepath )

	#=====================================================================
	def RankOrderedPredictedValuesGraph( self, chart_title=None ):
		"""This graph visualizes the distribution of predicted values generated by classification.
		For each individual ImageClassificationResult with ground truth value (i.e., class id) and
		predicted value, all results are grouped within their class, sorted by predicted value
		in ascending order, then ploted side-by-side.
		
		Required the package matplotlib to be installed."""

		print "Rendering rank-ordered predicted values graph"
		import matplotlib
		# Need following line to generate images on servers, see
		# http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		from itertools import cycle
		color = cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

		self.figure = plt.figure()
		self.main_axes = self.figure.add_subplot(111)

		if chart_title:
			self.chart_title = chart_title
		self.main_axes.set_title( self.chart_title )
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
	def KernelSmoothedDensityGraph( self, chart_title=None ):
		"""This graph visualizes the distribution of predicted values generated by classification.
		A kernel-smoothed probability density function is plotted for each image class on
		the same chart allowing comparison of distribution of predicted values amoung image class.
		
		Requires the packages matplotlib and scipy. Uses scipy.stats.gaussian_kde to
		generate kernel-smoothed probability density functions."""

		print "Rendering kernel-smoothed probability density estimate graph"
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		from itertools import cycle
		color = cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

		self.figure = plt.figure()
		self.main_axes = self.figure.add_subplot(111)
		if chart_title:
			self.chart_title = chart_title

		self.main_axes.set_title( self.chart_title )
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

	def __init__( self, training_set, feature_weights, test_image_path,
	    chart_title=None, max_num_features=300 ):

		self.timing_axes = None
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
				sig = FeatureVector( source_filepath=test_image_path, featurenames_list=reduced_fw.featurenames_list ).GenerateFeatures()
				reduced_ts = training_set.FeatureReduce( reduced_fw )
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

		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		x_vals = list( range( 1, max_num_features + 1 ) )

		self.figure = plt.figure()
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
class AccuracyVersusNumFeaturesGraph( BaseGraph ):
	"""Graphing the figure of merit a a function of number of features"""

	# FIXME: roll this class into FeatureTimingVersusAccuracyGraph, allowing
	# both Discrete and continuous data

	def __init__( self, training_set, feature_weights, chart_title=None, min_num_features=1, max_num_features=None, step=5, y_min=None, y_max=None, quiet=False):

		ls_experiment = ContinuousClassificationExperimentResult( training_set, training_set, feature_weights, name="Least Squares Regression Method")
		voting_experiment = ContinuousClassificationExperimentResult( training_set, training_set, feature_weights, name="Voting Method")
		if max_num_features is None:
			max_num_features = len( feature_weights )

		x_vals = range( min_num_features, max_num_features + 1, step )

		for number_of_features_to_use in x_vals:
			reduced_fw = feature_weights.Threshold( number_of_features_to_use )
			reduced_ts = training_set.FeatureReduce( reduced_fw )
			if not quiet:
				reduced_fw.Print()

			ls_batch_result = ContinuousBatchClassificationResult.NewLeastSquaresRegression( reduced_ts, None, reduced_fw, batch_number=number_of_features_to_use, quiet=my_quiet )
			if not quiet:
				ls_batch_result.Print()
			ls_experiment.individual_results.append( ls_batch_result )

			voting_batch_result = ContinuousBatchClassificationResult.New( reduced_ts, reduced_fw, batch_number=number_of_features_to_use )
                        if not quiet:
			    voting_batch_result.Print()
			voting_experiment.individual_results.append( voting_batch_result )

		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		self.figure = plt.figure( figsize=(12, 8) )
		self.main_axes = self.figure.add_subplot(111)
		if chart_title == None:
			self.chart_title = "R vs. num features, two methods"
		else:
			self.chart_title = chart_title

		# need to make axes have same range
		ls_yvals = [ batch_result.figure_of_merit for batch_result in ls_experiment.individual_results ]
		voting_yvals = [ batch_result.figure_of_merit for batch_result in voting_experiment.individual_results ]

		min_ls_yval = min(ls_yvals)
		optimal_num_feats_ls = ls_yvals.index( min_ls_yval ) + 1 # count from 1, not 0
		min_voting_yval = min(voting_yvals)
		optimal_num_feats_voting = voting_yvals.index( min_voting_yval ) + 1 # count from 1, not 0

		all_vals = ls_yvals + voting_yvals

		if y_min is not None:
			try:
				y_min = float(y_min)
			except:
				raise ValueError( "Can't convert {0} to float".format(y_min))
			_min = y_min
		else:
			_min = min( all_vals )
		if y_max is not None:
			try:
				y_max = float(y_max)
			except:
				raise ValueError( "Can't convert {0} to float".format(y_max))
			_max = y_max
		else:
			_max = max( all_vals )

		# Plot least Squares Data
		self.main_axes.set_title( self.chart_title )
		self.main_axes.set_xlabel( 'Number of features' )
		self.main_axes.set_ylabel( 'RMS Least Squares Method', color='b' )
		self.main_axes.set_ylim( [_min, _max ] )
		self.main_axes.plot( x_vals, ls_yvals, color='b', marker='o', linestyle='--' )
		for tl in self.main_axes.get_yticklabels():
			tl.set_color('b')

		self.main_axes.annotate( 'min R={0:.3f} @ {1}'.format(min_ls_yval, optimal_num_feats_ls),
		color='b',
		xy=( optimal_num_feats_ls, min_ls_yval ),
		xytext=( optimal_num_feats_ls, 0.8 * _max ),
		arrowprops=dict(facecolor='black', shrink=0.05),
		horizontalalignment='right' )

		# Plot Voting method data
		self.timing_axes = self.main_axes.twinx()
		self.timing_axes.set_ylabel( 'RMS Voting Method', color='r' )
		self.timing_axes.set_ylim( [_min, _max ] )
		self.timing_axes.plot( x_vals, voting_yvals, color='r', marker='o', linestyle='--' )
		for tl in self.timing_axes.get_yticklabels():
			tl.set_color('r')

		self.timing_axes.annotate( 'min R={0:.3f} @ {1}'.format(min_voting_yval, optimal_num_feats_voting),
		color='r',
		xy=( optimal_num_feats_voting, min_voting_yval ),
		xytext=( optimal_num_feats_voting, 0.6 * _max ),
		arrowprops=dict(facecolor='black', shrink=0.05),
		horizontalalignment='right' )

#============================================================================
class Dendrogram( object ):
	"""Not implemented. In the future might use scipy.cluster (no unrooted dendrograms though!)
	or Biopython.Phylo to visualize. Perhaps could continue C++ implementation's use of PHYLIP's
	Fitch-Margoliash program "fitch" to generate Newick phylogeny, and visualize using
	native python tools."""
	pass

#================================================================

initialize_module()

