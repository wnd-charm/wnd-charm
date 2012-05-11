#!/usr/bin/env python
import stringformat
import pychrm
import numpy as np
from StringIO import StringIO

#

def main():
	
	filename = '/Users/chris/projects/josiah_worms/terminal_bulb.fit'
	fitfile = open( filename )


	# FIXME: put all this in the TrainingSet constructor

	print "loading file {}".format( fitfile )
	num_classes = -1
	num_features = -1
	num_images = -1
	classnames_list = []
	featurenames_list = []
	tilenames_list = []

	data_lists = []

	name_line = False
	line_num = 0
	for line in fitfile:
		if line_num is 0:
			num_classes = int( line )
			# we sort samples according to their ground truth class on the fly as we read them in
			# Initialize the lists where the sorted samples will end up
			for i in range( num_classes ):
				data_lists.append( [] )
		elif line_num is 1:
			num_features = int( line )
		elif line_num is 2:
			num_images = int( line )
		elif line_num <= ( num_features + 2 ):
			featurenames_list.append( line.strip() )
		elif line_num == ( num_features + 3 ):
			pass # skip a line
		elif line_num <= ( num_features + 3 + num_classes ):
			classnames_list.append( line.strip() )
		else:
			# read in features
			if name_line:
				tilenames_list.append( line.strip() )
			if not name_line:
				# strip off the class identity value, which is the last in the array
				split_line = line.strip().rsplit( " ", 1)
				# The features are in split_line[0]
				# and the class_id is in split_line[1]
				#print "class {}".format( split_line[1] )
				data_lists[ int( split_line[1] ) - 1 ].append( split_line[0] )
			name_line = not name_line
		line_num += 1

	fitfile.close()

	string_data = "\n"
	
	sig_matrix = []
	for i in range( num_classes ):
		print "generating matrix for class {}".format( i )
		sig_matrix.append( np.genfromtxt( StringIO( string_data.join( data_lists[i] ) ) ) )

	# Return a fully constructed TrainingSet class

	# normalize the features at some point

	# read in weights from a wndchrm html file

	# now, give me only those columns that correspond with feature weights

	one_tile = sig_matrix[0][0]
	# do a classify operation

	fake_weights = [1] * num_features
	norm_factor, marg_probs = ClassifyWND5( sig_matrix, one_tile, fake_weights )

	print "norm factor {}, marg probs {}".format( norm_factor, marg_probs )


def ClassifyWND5( trainingset, testimg, feature_weights ):
	# If you're using this function, your training set data is not continuous
	# for N images and M features:
	#   trainingset is list of length L of N x M numpy matrices
	#   testtile is a 1 x M list of feature values
	# NOTE: the trainingset and test image must have the same number of features!!!
  # returns a list of length L comprised of marginal probabilities
	# FIXME: what about tiling??

	print "classifying..."
	num_test_img_features = len( testimg ) 
	
	epsilon = np.finfo( np.float ).eps
	#EpsTest = np.vectorize( lambda x: 0 if x < epsilon else x )

	weights_squared = np.square( feature_weights )

	#class_distances = [0]* len( trainingset )
	class_similarities = [0] * len( trainingset )

	for class_index in range( len( trainingset ) ):
		print "Calculating distances to class {}".format( class_index )
		num_tiles, num_features = trainingset[ class_index ].shape
		assert num_test_img_features == num_features,\
		"num features {}, num features in test img {}".format( num_features, num_test_img_features )

		# create a view
		sig_matrix = trainingset[ class_index ]
		wnd_sum = 0
		num_collisions = 0

		#print "num tiles: {}, num_test_img_features {}".format( num_tiles, num_test_img_features )
		for tile_index in range( (num_tiles) ):
			#print "{} ".format( tile_index )
			# dists = EpsTest( np.absolute( sig_matrix[ tile_index ] - testimg ) )
			# epsilon checking for each feature is too expensive
			# do this quick and dirty check until we can figure something else out
			dists = np.absolute( sig_matrix[ tile_index ] - testimg )
			w_dist = np.sum( dists )
			if w_dist < epsilon:
				num_collisions += 1
				continue
			dists = np.multiply( weights_squared, np.square( dists ) )
			w_dist = np.sum( dists )
			#class_distances[ class_index ] += w_dist
			class_similarities[ class_index ] += w_dist ** -5
		#print "\n"

		#class_distances[ class_index ] /= num_tiles - num_collisions
		class_similarities[ class_index ] /= num_tiles - num_collisions

	normalization_factor = sum( class_similarities )

	return ( normalization_factor, [ x / normalization_factor for x in class_similarities ] ) 

	

#================================================================
if __name__=="__main__":
	main()


