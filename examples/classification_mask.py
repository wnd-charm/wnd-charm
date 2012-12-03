#!/usr/bin/env python

# import pychrm
from pychrm.TrainingSet import *
import sys
import os
import re 

# Give a nice error message in case PIL isn't installed
try:
	from PIL import Image
except:
	print "This script depends on the Python Image Library (PIL) to write mask images."
	print "Please download and install PIL from <http://www.pythonware.com/products/pil/>"
	print "Or, use your system's software package manager to install PIL for Python"
	print "   Ubuntu, Debian and variants: sudo apt-get install python-imaging"
	print "   CentOS, RedHat and variants: sudo yum install python-imaging"
	sys.exit(0)

# We'll be using numpy to work with image masks
import numpy as np

# We're doing manual parameter processing, which is probably not a great idea...
if ( len(sys.argv) < 4 ):
	print "Generate mask images by running a classifier over an input image."
	print "Specify a classifier (.fit file and number of features or pickled features and weights),"
	print "the size of the window to scan, and an input tiff file"
	print "The output will be a set of masks stored as tiff files (one per class), where the pixel values are"
	print "the marginal probabilities of that class in the corresponding window location"
	print "Usage:"
	print "\t"+sys.argv[0]+" (classifier.fit [num features] | train_features.pickled feature_weights.pickled) size_xXsize_y input.tif"
	sys.exit(0)

from_scratch = False
pickled_features = None
pickled_weights = None
num_features = 200

# Get the classifier parameter(s)
input_filename = sys.argv[1]
next_arg=2
if ( input_filename.endswith (".fit") ):
	from_scratch = True
	if ( len (sys.argv) > 2 and sys.argv[next_arg].isdigit() ):
		num_features = int (sys.argv[next_arg])
		next_arg = next_arg + 1
	print "using top "+str (num_features)+" features"
elif ( input_filename.endswith (".fit.pickled") ):
	from_scratch = False
	pickled_features = input_filename
	try:
		if (sys.argv[next_arg].endswith (".weights.pickled")):
			pickled_weights = sys.argv[next_arg]
			next_arg = next_arg + 1
		else:
			raise Exception("expecting second argument to be pickled weights")
	except:
		print "expecting second argument to be pickled weights"
		sys.exit(0)
else:
	print "first argument is either a .fit file or pickled features and weights"

# Get the size of the window
try:
	re_size = re.search(r"^(\d+)[xX](\d+)$", sys.argv[next_arg])
	if ( re_size ):
		scan_x = int (re_size.groups()[0])
		scan_y = int (re_size.groups()[1])
		print "window size = ("+str(scan_x)+","+str(scan_y)+")"
		next_arg = next_arg + 1
	else:
		raise Exception("no scanning window arg")
except:
		print "expecting an argument specifying the width X height of scanning window (e.g. 20x20)"
		sys.exit(0)

# Get the input image
if re.search(r"\.tiff?$", sys.argv[next_arg], re.IGNORECASE):
	image_path = sys.argv[next_arg]
	if not os.path.exists( image_path ):
		raise ValueError( "The file '{0}' doesn't exist, maybe you need to specify the full path?".format( image_path ) )
	input_image = pychrm.ImageMatrix()
	if 1 != input_image.OpenImage( image_path, 0, None, 0, 0 ):
		raise ValueError( 'Could not build an ImageMatrix from {0}, check the file.'.format( image_path ) )
	next_arg = next_arg + 1
else:
	print "expecting an input image file (.tif, .tiff, .TIF, .TIFF)"
	sys.exit(0)





if from_scratch:
	# I preprocessed your training set and feature weights and pickled them for speed.
	# Pickle files are binary files that are super fast to load.
	# You don't need to use a pickle file though, you can make one from scratch
	# Here's how:

	# 1. Load the raw c-charm fit file
	full_training_set = DiscreteTrainingSet.NewFromFitFile( input_filename )

	# 2. C-charm uses "Lior-style" feature names. Translate them into the new "Ilya-style"
	full_training_set.featurenames_list = FeatureNameMap.TranslateToNewStyle( full_training_set.featurenames_list )

	# 3. Normalize the features:
	full_training_set.Normalize()

	# 4. Make Fisher scores based on the normalized training set
	full_fisher_weights = FisherFeatureWeights.NewFromTrainingSet( full_training_set )

	# 5. Take only the top 200 features
	reduced_fisher_weights = full_fisher_weights.Threshold( num_features )

	# 6. Reduce the training set feature space to contain only those top 200 features
	reduced_training_set = full_training_set.FeatureReduce( reduced_fisher_weights.names )

	# 7. Save your work:
	reduced_training_set.PickleMe( os.path.splitext(input_filename)[0] + ".fit.pickled" )
	reduced_fisher_weights.PickleMe( os.path.splitext(input_filename)[0] + "_w"+str(num_features)+".weights.pickled" )
else:
	# I've already done all that, just proceed from here:
	reduced_training_set = DiscreteTrainingSet.NewFromPickleFile( pickled_features )
	reduced_fisher_weights = FisherFeatureWeights.NewFromPickleFile( pickled_weights )


# create the tile image iterator
image_iter = SampleImageTiles (input_image, scan_x, scan_y, True)
print "Number of samples = "+str (image_iter.samples)

# Create a list of zero'd out, image-sized, 2-D byte numpys
masks = [] 
for i in range( reduced_training_set.num_classes ):
	masks.append ( np.zeros (shape=(image_iter.image.height,image_iter.image.width), dtype='uint8') )


# iterate over the image, classifying each tile
for sample in image_iter.sample():
	try:
		test_image_signatures = Signatures.NewFromFeatureNameList( sample, reduced_fisher_weights.names )
		test_image_signatures.Normalize( reduced_training_set )
		result = DiscreteImageClassificationResult.NewWND5( reduced_training_set, reduced_fisher_weights, test_image_signatures )
		for i in range( reduced_training_set.num_classes ):
			mask_val = int (result.marginal_probabilities[i] * 255.0)
			# Write the mask value into the numpy
			masks[i][image_iter.current_y:image_iter.current_y + image_iter.tile_height,
				image_iter.current_x:image_iter.current_x + image_iter.tile_width] = mask_val
			print "{0} ({1},{2}) {3}: {4}".format (
				i, image_iter.current_x, image_iter.current_y, reduced_training_set.classnames_list[i], mask_val)
	except:
		x_y_str = "{0}_{1}".format (image_iter.current_x, image_iter.current_y)
		tif_path = os.path.join (os.path.abspath(os.path.dirname(input_filename)), os.path.splitext(os.path.basename(image_path))[0] + "_" + x_y_str + ".tiff")
		print "Could not classify sample at ({0},{1}), saving tiff file {2}".format (
			image_iter.current_x, image_iter.current_y, tif_path)
		sample.SaveTiff (tif_path)
		for i in  range ( len(test_image_signatures.values) ):
			val = test_image_signatures.values[i]
			if np.isinf(val):
				print test_image_signatures.names[i]+" is INF"
			elif np.isnan(val):
				print test_image_signatures.names[i]+" is NAN"
		sys.exit(0)

# 
# Create tiff files from the numpys
mask_dir = os.path.abspath(os.path.dirname(input_filename))
print "mask tiff files will be saved in '{0}{1}'".format(mask_dir, os.sep)
for i in range( reduced_training_set.num_classes ):
	class_name = reduced_training_set.classnames_list[i]
	mask_path = os.path.join (mask_dir, os.path.splitext(os.path.basename(image_path))[0] + "_" + class_name + ".tiff")
	print 'creating tiff file mask for class {0}: {1}'.format(class_name, os.path.basename(mask_path) )
	# Make a PIL image out of the numpy and save it as a tiff.
	Image.fromarray(masks[i]).save(mask_path)
