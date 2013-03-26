#!/usr/bin/env python
""" Real-time image classification script using Pychrm
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

This is an example script that illustrates how to use (or create) a canned image classifier
to classify a single image.

"""

import sys
import os
import re 

# import pychrm
from pychrm.FeatureSet import *
from pychrm import __version__ as pychrm_version
print "pychrm "+pychrm_version

# We're doing manual parameter processing, which is probably not a great idea...
if ( len(sys.argv) < 3 ):
	print "Generate mask images by running a classifier over an input image."
	print "Specify a classifier (.fit file and number of features or pickled features and weights),"
	print "the size of the window to scan, and an input tiff file"
	print "The output will be a set of masks stored as tiff files (one per class), where the pixel values are"
	print "the marginal probabilities of that class in the corresponding window location"
	print "Usage:"
	print "\t"+sys.argv[0]+" (classifier.fit [num features] | train_features.pickled feature_weights.pickled) input.tif"
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

# Get the input image
if re.search(r"\.tiff?$", sys.argv[next_arg], re.IGNORECASE):
	image_path = sys.argv[next_arg]
	if not os.path.exists( image_path ):
		raise ValueError( "The file '{0}' doesn't exist, maybe you need to specify the full path?".format( image_path ) )
	next_arg = next_arg + 1
else:
	print "expecting an input image file (.tif, .tiff, .TIF, .TIFF)"
	sys.exit(0)

# For real time classification, it is best practice to preprocess your FeaturesSet and
# feature weights and pickle them for speed.
# Pickle files are binary files that are super fast to load.
# You don't need to use a pickle file though, you can make one from scratch
# Here's how:

if from_scratch:
	 
	# 1a. Instantiate a FeaturesSet from a file, perhaps a ".fit" file from the
	#    legacy C++ WND-CHARM implementation (a.k.a. "C-charm")
	full_training_set = FeatureSet_Discrete.NewFromFitFile( input_filename )

	# 2. Normalize the features:
	full_training_set.Normalize()

	# 3. Make Fisher scores based on the normalized training set
	full_fisher_weights = FisherFeatureWeights.NewFromFeatureSet( full_training_set )

	# 4. Take only the top 200 features
	reduced_fisher_weights = full_fisher_weights.Threshold( num_features )

	# 5. Reduce the training set feature space to contain only those top 200 features
	reduced_training_set = full_training_set.FeatureReduce( reduced_fisher_weights.names )

	# 6. Save your work:
	reduced_training_set.PickleMe( os.path.splitext(input_filename)[0] + "_w"+str(num_features) + ".fit.pickled" )
	reduced_fisher_weights.PickleMe( os.path.splitext(input_filename)[0] + "_w"+str(num_features)+".weights.pickled" )

else:
	# If you've already done all that, just proceed from here:
	reduced_training_set = FeatureSet_Discrete.NewFromPickleFile( pickled_features )
	reduced_fisher_weights = FisherFeatureWeights.NewFromPickleFile( pickled_weights )


# Calculate features for the test image, but only those features we need 
test_image_signatures = Signatures.NewFromFeatureNameList( image_path, reduced_fisher_weights.names )

# It might be useful to hold onto the sigs, so write them out to a file
test_image_signatures.WriteFeaturesToASCIISigFile()
# Note that if you call the function without specifying a filename, it will generate a
# default name and path for the sigs file based on the original image
# In the future, if you want to load the sig file, use this function call:
future = False
if future:
	test_image_signatures = Signatures.NewFromSigFile( "whatever_the_path_is.pysig" )

# Normalize the newly calculated signatures against the training set
test_image_signatures.Normalize( reduced_training_set )

# Classify away! Return all the pertinent results including marginal probabilities,
# normalization factor, and interpolated value inside the variable "result"
result = DiscreteImageClassificationResult.NewWND5( reduced_training_set, reduced_fisher_weights, test_image_signatures )

# See what we got... Print out the results to STDOUT
result.Print()

