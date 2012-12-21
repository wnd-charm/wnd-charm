#!/usr/bin/env python
"""
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

Meant to exercize train/test/split functionality"""


# import pychrm
from pychrm.FeatureSet import *
from pychrm import __version__ as pychrm_version
print "pychrm "+pychrm_version

import argparse

parser = argparse.ArgumentParser( description="perform classifier cross validation" )
parser.add_argument( '-n', help='specify number of train/test splits',
                     type=int, metavar='<integer>', default=5 )
parser.add_argument( '-f', help='specify number feature usage fraction on interval [0.0,1.0]',
                     type=float, metavar='<float>',default = None)
parser.add_argument( '-F', help='specify number of features',
                     type=int, metavar='<integer>',default = 200)
parser.add_argument( 'classifier_file_path', help='path to Pychrm classifier file, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'output_filepath', help='Results are written to this file, otherwise to STDOUT',
                     nargs='?', default = None )
args = parser.parse_args()

num_splits = args.n
feature_usage_fraction = args.f
num_features = args.F
input_filename = args.classifier_file_path[0]
outpath = args.output_filepath

full_training_set = None

# Get the classifier parameter(s)
if ( input_filename.endswith (".fit") ):
	full_training_set = FeatureSet_Discrete.NewFromFitFile( input_filename )
	full_training_set.featurenames_list = FeatureNameMap.TranslateToNewStyle( full_training_set.featurenames_list )
elif ( input_filename.endswith (".fit.pickled") ):
	full_training_set = FeatureSet_Discrete.NewFromPickleFile( input_filename )
elif ( input_filename.endswith (".fof") ):
	full_training_set = FeatureSet_Discrete.NewFromFileOfFiles( input_filename )
else:
	raise Exception( 'The classifier must either end in .fit, .fit.pickled, or .fof' )

full_training_set.Print()

if feature_usage_fraction:
	if feature_usage_fraction < 0 or feature_usage_fraction > 1.0:
		raise Exception('Feature usage fraction must be on interval [0,1]')
	num_features = int( feature_usage_fraction * full_training_set.num_features )
print "using top "+str (num_features)+" features"

experiment = DiscreteClassificationExperimentResult( training_set = full_training_set )


for i in range( num_splits ):

	training_set, test_set = full_training_set.Split()
	training_set.Normalize()
	test_set.Normalize( training_set )

	fisher_weights = FisherFeatureWeights.NewFromFeatureSet( training_set )
	fisher_weights = fisher_weights.Threshold( num_features )

	reduced_test_set = test_set.FeatureReduce( fisher_weights.names )
	reduced_training_set = training_set.FeatureReduce( fisher_weights.names )

	batch_result = DiscreteBatchClassificationResult.New( reduced_training_set, \
	                                    reduced_test_set, fisher_weights, batch_number = i )

	experiment.individual_results.append( batch_result )

if outpath:
	experiment.Print( output_filepath=outpath, mode='w' )
	experiment.PredictedValueAnalysis( output_filepath=outpath, mode= 'a' )
else:
	experiment.Print()
	experiment.PredictedValueAnalysis()

