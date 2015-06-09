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
 Written by:  Christopher Coletta (github.com/colettace)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Meant to exercize train/test/split functionality"""

def get_featureset( path, quiet=True ):
    if path.endswith( ".fit" ):
        fs = FeatureSpace.NewFromFitFile( path )
    elif path.endswith( ".fit.pickled" ):
        fs = FeatureSpace.NewFromPickleFile( path )
    elif path.endswith( ".fof" ):
        fs = FeatureSpace.NewFromFileOfFiles( path )
    else:
        raise Exception( 'The classifier must either end in .fit, .fit.pickled, or .fof' )
    if not quiet:
        fs.Print()
    return fs

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights
from wndcharm.FeatureSpacePrediction import FeatureSpaceClassification
from wndcharm.FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment

from wndcharm import __version__ as wndcharm_version
print "wndcharm "+wndcharm_version

import argparse

parser = argparse.ArgumentParser( description="obtain classifier *training* accuracy" )
parser.add_argument( '-n', help='specify number of train/test splits',
                     type=int, metavar='<integer>', default=1 )
parser.add_argument( '-f', help='specify number feature usage fraction on interval [0.0,1.0]',
                     type=float, metavar='<float>',default=None)
parser.add_argument( '-F', help='specify number of features',
                     type=int, metavar='<integer>',default=None)
parser.add_argument( 'training_set_file_path', help='path to WND-CHARM classifier file for training, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'test_set_file_path', help='path to WND-CHARM classifier file for testing, could be WND-CHARM .fit file, Pychrm .fit.pickled file, or Pychrm/WND-CHRM file of files .fof',
                     nargs=1 )
parser.add_argument( 'output_filepath', help='Results are written to this file, otherwise to STDOUT',
                     nargs='?', default=None )
parser.add_argument( '-D', help='write all intermediates to disk, including .sig and .fit',
                     action='store_true', default=False)
args = parser.parse_args()

num_splits = args.n
feature_usage_fraction = args.f
num_features = args.F
training_filename = args.training_set_file_path[0]
testing_filename = args.test_set_file_path[0]
outpath = args.output_filepath
write_intermediates = args.D

train_set = get_featureset( training_filename )
if write_intermediates:
    train_set.ToFitFile()

if training_filename == testing_filename:
    test_set = train_set
else:
    test_set = get_featureset( testing_filename )
    if write_intermediates:
        train_set.ToFitFile()

if feature_usage_fraction:
    if feature_usage_fraction < 0 or feature_usage_fraction > 1.0:
        raise Exception('Feature usage fraction must be on interval [0,1]')
    num_features = int( feature_usage_fraction * train_set.num_features )

if num_features:
    print "Using top {0} Fisher-ranked features.".format( num_features )
else:
    print "Using top 15% Fisher-ranked features."

experiment = FeatureSpaceClassificationExperiment( training_set=train_set )

train_set.Normalize( inplace=True )
weights = FisherFeatureWeights.NewFromFeatureSpace( train_set ).Threshold( num_features )
train_set.FeatureReduce( weights, inplace=True )

if train_set != test_set:
    test_set.FeatureReduce( weights, inplace=True ).Normalize( train_set )

for i in range( num_splits ):
    split = FeatureSpaceClassification.NewWND5( train_set, test_set, weights, batch_number=i )
    experiment.individual_results.append( split )

if outpath:
    experiment.Print( output_filepath=outpath, mode='w' )
    #experiment.PerSampleStatistics( output_filepath=outpath, mode= 'a' )
else:
    experiment.Print()
    #experiment.PerSampleStatistics()

