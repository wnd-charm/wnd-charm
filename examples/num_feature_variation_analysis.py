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


This script performs discrete classifier cross-validation while varying the number
of features used. For each range of features, the script produces Rank Ordered Predicted
Values graphs as well as Kernel-Smoothed Probability Density graphs showing how predicted 
values change as a function of number of features used. Finally the graphs are coallated into
a single PDF using ImageMagick's convert command."""


# import pychrm
from pychrm.FeatureSet import *
from pychrm import __version__ as pychrm_version
print "pychrm "+pychrm_version

import argparse

parser = argparse.ArgumentParser( description="Cross-validate a continuous classifier across a range of features, and produce graphs showing how predicted values change as a function of feature range used." )
parser.add_argument( '-n', help='specify number of train/test splits per feature range',
                     type=int, metavar='<integer>', default=5 )
parser.add_argument( '--start', help='Defining the low end of the feature number range',
                     type=int, metavar='<integer>', default = 15)
parser.add_argument( '--end', help='Defining the hi end of the feature number range',
                     type=int, metavar='<integer>', default = 431)
parser.add_argument( '--step', help='Defining the step inside feature number range',
                     type=int, metavar='<integer>', default = 25)
parser.add_argument( '-b', help='number of bins that the feature space is divided into for the purposes of creating a range of features',
                     type=int, metavar='<integer>', default = 10)
parser.add_argument( 'classifier_file_path', help='path to Pychrm classifier file, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'output_filepath', help='Results are written to this file, otherwise to STDOUT',
                     nargs='?')
parser.add_argument( '-D', help='Write the training set to a pickle file', metavar='<optional path>',
                     default='unset', nargs='?')

args = parser.parse_args()


num_splits = args.n
num_bins = args.b
input_filename = args.classifier_file_path[0]
outpath = args.output_filepath
dump_pickle = args.D
range_step = args.step
range_start = args.start
range_end = args.end

num_bins = (range_end - range_start) / range_step


if ( input_filename.endswith (".fit") ):
	full_set = FeatureSet_Discrete.NewFromFitFile( input_filename )
	full_set.featurenames_list = FeatureNameMap.TranslateToNewStyle( full_set.featurenames_list )
elif ( input_filename.endswith (".fit.pickled") ):
	full_set = FeatureSet_Discrete.NewFromPickleFile( input_filename )
elif ( input_filename.endswith (".fof") ):
	full_set = FeatureSet_Discrete.NewFromFileOfFiles( input_filename )
else:
	raise Exception( 'The classifier must either end in .fit, .fit.pickled, or .fof' )


if not dump_pickle == 'unset':
	if dump_pickle:
		# user used -D to specify a name for their training set pickle
		full_set.PickleMe( dump_pickle )
	else:
		# user used -D as a flag, use default pickle name pattern
		full_set.PickleMe()


rank_ordered_graphs = []
ks_density_graphs = []

for bin_index in range( num_bins + 1 ):
	num_features = range_start + bin_index * range_step
	name = "{0}\nusing {1} features".format( full_set.source_path, num_features )
	print name
	experiment = DiscreteClassificationExperimentResult( name )
	experiment.test_set = full_set

	for i in range( num_splits ):

		full_training_set, full_test_set = full_set.Split( quiet = True, balanced_classes = True )
		full_training_set.Normalize( quiet = True )
		full_test_set.Normalize( full_training_set, quiet = True )

		full_weights = FisherFeatureWeights.NewFromFeatureSet( full_training_set )
		# This just does a reorder
		weights_subset = full_weights.Threshold( num_features )

		reduced_training_set = full_training_set.FeatureReduce( weights_subset.names )
		reduced_test_set = full_test_set.FeatureReduce( weights_subset.names )

		batch_result = DiscreteBatchClassificationResult.New( reduced_training_set, reduced_test_set, \
		                                 weights_subset, quiet = True, batch_number = i )

		experiment.individual_results.append( batch_result )

	experiment.Print()

	grapher = PredictedValuesGraph( experiment )
	grapher.RankOrderedPredictedValuesGraph( name )
	grapher.SaveToFile( "rank_ordered_features_{0:03d}".format( num_features ) )

	grapher.KernelSmoothedDensityGraph( name )
	grapher.SaveToFile( "ks_density_features_{0:03d}".format( num_features ) )

	bin_offset += num_features_per_bin 


import subprocess
subprocess.call( [ "convert", "rank_ordered_features*", "rank_ordered.pdf"] )
subprocess.call( [ "convert", "ks_density_features*", "ks_density.pdf"] )

