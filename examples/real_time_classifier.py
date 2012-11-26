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

# Pull the input image's filename from the command line
import sys
input_filename = sys.argv[1]

# import pychrm
from pychrm import *

# For real time classification, it is best practice to preprocess your FeaturesSet and
# feature weights and pickle them for speed.
# Pickle files are binary files that are super fast to load.
# You don't need to use a pickle file though, you can make one from scratch
# Here's how:

from_scratch = False
if from_scratch:
	 
	# 1a. Instantiate a FeaturesSet from a file, perhaps a ".fit" file from the
	#    legacy C++ WND-CHARM implementation (a.k.a. "C-charm")
	full_training_set = FeaturesSet_Discrete.NewFromFitFile( "OfficialState.fit" )

	# 1b. Translate feature names from C-chrm style, to Pychrm style
	FeatureNameMap.TranslateToNewStyle( full_training_set )

	# 2. Normalize the features:
	full_training_set.Normalize()

	# 3. Make Fisher scores based on the normalized training set
	full_fisher_weights = FisherFeatureWeights.NewFromTrainingSet( full_training_set )

	# 4. Take only the top 200 features
	fisher_weights_subset = full_fisher_weights.Threshold( 200 )

	# 5. Reduce the training set feature space to contain only those top 200 features
	reduced_training_set = fisher_weights_subset.FeatureReduce( full_training_set )

	# 6. Save your work:
	reduced_training_set.PickleMe( "OfficialState_normalized_200_features.fit.pickled" )
	fisher_weights_subset.PickleMe( "feature_weights_len_200.weights.pickled" )


# If you've already done all that, just proceed from here:
reduced_training_set = FeaturesSet_Discrete.NewFromPickleFile( "OfficialState_normalized_200_features.fit.pickled" )
fisher_weights_subset = FisherFeatureWeights.NewFromPickleFile( "feature_weights_len_200.weights.pickled" )


# Calculate features for the test image, but only those features we need 
test_image_signatures = Signatures.NewFromFeatureWeights( input_filename, fisher_weights_subset )

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
result = DiscreteImageClassificationResult.NewWND5( reduced_training_set, fisher_weights_subset, test_image_signatures )

# See what we got... Print out the results to STDOUT
result.Print()

