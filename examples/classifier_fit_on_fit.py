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

def get_featureset (input_filename):
	# Get the classifier parameter(s)
	if ( input_filename.endswith (".fit") ):
		featureset = FeatureSet_Discrete.NewFromFitFile( input_filename )
	elif ( input_filename.endswith (".fit.pickled") ):
		featureset = FeatureSet_Discrete.NewFromPickleFile( input_filename )
	elif ( input_filename.endswith (".fof") ):
		featureset = FeatureSet_Discrete.NewFromFileOfFiles( input_filename )
	else:
		raise Exception( 'The classifier must either end in .fit, .fit.pickled, or .fof' )
	featureset.Print()
	return (featureset)


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
parser.add_argument( 'training_set_file_path', help='path to Pychrm classifier file for training, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'test_set_file_path', help='path to Pychrm classifier file for testing, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'output_filepath', help='Results are written to this file, otherwise to STDOUT',
                     nargs='?', default = None )
args = parser.parse_args()

num_splits = args.n
feature_usage_fraction = args.f
num_features = args.F
training_filename = args.training_set_file_path[0]
testing_filename = args.test_set_file_path[0]
outpath = args.output_filepath

training_featureset = get_featureset (training_filename)
testing_featureset = get_featureset (testing_filename)



if feature_usage_fraction:
	if feature_usage_fraction < 0 or feature_usage_fraction > 1.0:
		raise Exception('Feature usage fraction must be on interval [0,1]')
	num_features = int( feature_usage_fraction * training_featureset.num_features )
print "using top "+str (num_features)+" features"

experiment = DiscreteClassificationExperimentResult( training_set = training_featureset )

training_featureset.Normalize()
testing_featureset.Normalize( training_featureset )

fisher_weights = FisherFeatureWeights.NewFromFeatureSet( training_featureset )
fisher_weights = fisher_weights.Threshold( num_features )

reduced_test_set = testing_featureset.FeatureReduce( fisher_weights.names )
reduced_training_set = training_featureset.FeatureReduce( fisher_weights.names )

batch_result = DiscreteBatchClassificationResult.New( reduced_training_set,
	reduced_test_set, fisher_weights, batch_number = 0 )


experiment.individual_results.append( batch_result )

if outpath:
	experiment.Print( output_filepath=outpath, mode='w' )
	experiment.PredictedValueAnalysis( output_filepath=outpath, mode= 'a' )
else:
	experiment.Print()
	experiment.PredictedValueAnalysis()

