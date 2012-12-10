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

Test of various continuous classification functionality.

This script performs continuous classifier cross-validation while varying the range
of features used. For each range of features, the script produces Rank Ordered Predicted
Values graphs as well as Kernel-Smoothed Probability Density graphs showing how predicted 
values change as a function of feature range used. Finally the graphs are coallated into
a single PDF using ImageMagick's convert command."""


# import pychrm
from pychrm.TrainingSet import *
from pychrm import __version__ as pychrm_version
print "pychrm "+pychrm_version

import argparse

parser = argparse.ArgumentParser( description="Cross-validate a continuous classifier across a range of features, and produce graphs showing how predicted values change as a function of feature range used." )
parser.add_argument( '-n', help='specify number of train/test splits per feature range',
                     type=int, metavar='<integer>', default=5 )
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

if ( input_filename.endswith (".fit") ):
	full_set = FeatureSet_Continuous.NewFromFitFile( input_filename )
	full_set.featurenames_list = FeatureNameMap.TranslateToNewStyle( full_set.featurenames_list )
elif ( input_filename.endswith (".fit.pickled") ):
	full_set = FeatureSet_Continuous.NewFromPickleFile( input_filename )
elif ( input_filename.endswith (".fof") ):
	full_set = FeatureSet_Continuous.NewFromFileOfFiles( input_filename )
else:
	raise Exception( 'The classifier must either end in .fit, .fit.pickled, or .fof' )


if not dump_pickle == 'unset':
	if dump_pickle:
		# user used -D to specify a name for their training set pickle
		full_set.PickleMe( dump_pickle )
	else:
		# user used -D as a flag, use default pickle name pattern
		full_set.PickleMe()

num_features_per_bin = int( float( len( full_set.featurenames_list ) ) / float( num_bins) )
bin_offset = 0

rank_ordered_graphs = []
ks_density_graphs = []

for bin_index in range( num_bins ):
	name = "{0}\nFeatures ranked {1}-{2}".format( full_set.source_path,
	                                     bin_offset + 1, bin_offset + num_features_per_bin )
	print name
	experiment = ContinuousClassificationExperimentResult( name )
	experiment.test_set = full_set

	for i in range( num_splits ):

		full_training_set, full_test_set = full_set.Split( quiet = True )
		full_training_set.Normalize( quiet = True )
		full_test_set.Normalize( full_training_set, quiet = True )

		full_weights = ContinuousFeatureWeights.NewFromFeatureSet( full_training_set )
		# This just does a reorder
		full_weights = full_weights.Threshold( len( full_weights.names ) )
		weights_subset = full_weights.Slice( bin_offset, bin_offset + num_features_per_bin -1 )
		#weights_subset.Print( print_legend = False )

		reduced_training_set = full_training_set.FeatureReduce( weights_subset.names )
		reduced_test_set = full_test_set.FeatureReduce( weights_subset.names )

		batch_result = ContinuousBatchClassificationResult.New( reduced_test_set, \
		                                 weights_subset, quiet = True, batch_number = i )

		experiment.individual_results.append( batch_result )

	experiment.Print()

	grapher = PredictedValuesGraph( experiment )
	grapher.RankOrderedPredictedValuesGraph( name )
	grapher.SaveToFile( "rank_ordered_features_{0:03d}-{1:03d}".format( bin_offset + 1, bin_offset + num_features_per_bin  ) )

	grapher.KernelSmoothedDensityGraph( name )
	grapher.SaveToFile( "ks_density_features_{0:03d}-{1:03d}".format( bin_offset + 1, bin_offset + num_features_per_bin ) )

	bin_offset += num_features_per_bin 


import subprocess
subprocess.call( [ "convert", "rank_ordered_features*", "rank_ordered.pdf"] )
subprocess.call( [ "convert", "ks_density_features*", "ks_density.pdf"] )

