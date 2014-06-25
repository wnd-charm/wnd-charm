#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2013 University of Dundee & Open Microscopy Environment.
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#
#

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

import re
import numpy as np

from os.path import sep, dirname, realpath
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

from pychrm.FeatureSet import FeatureSet, FeatureSet_Discrete, FisherFeatureWeights,\
        DiscreteBatchClassificationResult, FeatureSet_Continuous, ContinuousFeatureWeights,\
        ContinuousBatchClassificationResult

class TestFeatureSet( unittest.TestCase ):
    """
    Test the pychrm module to check it still works with OmeroPychrm
    """

    def setUp(self):
        pass

    def test_ContinuousFitOnFit( self ):
        from pychrm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete

        fs_discrete = CreateArtificialFeatureSet_Discrete( n_samples=1000, n_classes=10,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True)

        tempdir = mkdtemp()
        path_to_fit = tempdir + sep + 'Artificial.fit'
 
        try:
          fs_discrete.ToFitFile( path_to_fit )
          fs_continuous = FeatureSet_Continuous.NewFromFitFile(
                  path_to_fit, discrete=False )

          fs_continuous.Normalize()
          fw_reduced = ContinuousFeatureWeights.NewFromFeatureSet( fs_continuous ).Threshold()
          fs_reduced = fs_continuous.FeatureReduce( fw_reduced.names )
          batch_result = ContinuousBatchClassificationResult.New( fs_reduced, fw_reduced )
          #batch_result.Print()

        finally:
          rmtree( tempdir )

    def test_DiscreteTrainTestSplitNoTiling( self ):
        """Uses binucleate test set"""

        fitfilepath = wndchrm_test_dir + sep + 'test-l.fit'
        fs = FeatureSet_Discrete.NewFromFitFile( fitfilepath )

        from numpy.random import RandomState
        prng = RandomState(42)
        full_train, full_test = fs.Split( random_state=prng )
        full_train.Normalize()
        reduced_fw = FisherFeatureWeights.NewFromFeatureSet( full_train ).Threshold()
        reduced_train = full_train.FeatureReduce( reduced_fw.names )

        reduced_test = full_test.FeatureReduce( reduced_fw.names )
        reduced_test.Normalize( reduced_train )

        batch_result = DiscreteBatchClassificationResult.New( reduced_train,
            reduced_test, reduced_fw )

    def test_DiscreteTrainTestSplitWithTiling( self ):
        """Uses a curated subset of the IICBU 2008 Lymphoma dataset, preprocessed as follows:
        auto-deconvolved, eosin channel only, tiled 5x6, 3 classes, 10 imgs per class,
        300 samples per class.
        """

        # Inflate the zipped test fit into a temp file
        import zipfile
        zipped_file_path = pychrm_test_dir + sep + 'lymphoma_t5x6_10imgseach.fit.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf.extractall( tempdir )

        try:
            fitfilepath = tempdir + sep + zf.namelist()[0]
            #fs = FeatureSet.NewFromFitFile( fitfilepath  )
            fs = FeatureSet_Discrete.NewFromFitFile( fitfilepath, tile_options=(5,6) )
            from numpy.random import RandomState
            prng = RandomState(42)
            #fs.Print( verbose=True )
            #print "\n\n\n********************\n\n\n"
            full_train, full_test = fs.Split( random_state=prng )
            #full_train.Print( verbose=True )
            #full_test.Print( verbose=True )
            full_train.Normalize()
            fw = FisherFeatureWeights.NewFromFeatureSet( full_train ).Threshold()
            reduced_train = full_train.FeatureReduce( fw.names )
            reduced_test = full_test.FeatureReduce( fw.names )
            reduced_test.Normalize( reduced_train )

        finally:
            rmtree( tempdir )

    def test_ContinuousTrainTestSplitWithTiling( self ):
        """Uses a synthetic preprocessed as follows: 500 total samples, 25 tiles per group
        240 total features"""

        from pychrm.ArtificialFeatureSets import CreateArtificialFeatureSet_Continuous

        fs = CreateArtificialFeatureSet_Continuous( n_samples=500,
                num_features_per_signal_type=20, n_samples_per_group=25 )

        from numpy.random import RandomState
        prng = RandomState(42)
        #fs.Print( verbose=True )
        #print "\n\n\n********************\n\n\n"
        full_train, full_test = fs.Split( random_state=prng )
        #full_train.Print( verbose=True )
        #full_test.Print( verbose=True )

        full_train.Normalize()
        fw = ContinuousFeatureWeights.NewFromFeatureSet( full_train ).Threshold()
        reduced_train = full_train.FeatureReduce( fw.names )
        reduced_test = full_test.FeatureReduce( fw.names )
        reduced_test.Normalize( reduced_train )

    def test_FitOnFitClassification( self ):

        fitfile_path = wndchrm_test_dir + sep + 'test-l.fit'
        #fs = FeatureSet.NewFromFitFile( fitfile_path )
        fs = FeatureSet_Discrete.NewFromFitFile( fitfile_path )
        fs.Normalize()
        reduced_fw = FisherFeatureWeights.NewFromFeatureSet( fs ).Threshold()
        reduced_fs = fs.FeatureReduce( reduced_fw.names )
        batch_result = DiscreteBatchClassificationResult.New( 
                                       reduced_fs, reduced_fs, reduced_fw )
        batch_result.Print()

    @unittest.skip( "test tile options after it's implemented" )
    def test_TileOptions( self ):

        tile_options = None
        fs = FeatureSet.NewFromFitFile( wndchrm_test_dir + sep + 'test-l.fit', tile_options  )

        # pass an int
        # pass a tuple
        # pass a None
        # Wht if tile options passed doesn't match fit file?

    @unittest.skip( "test split options after it's implemented" )
    def test_SplitOptions( self ):

        tile_options = None
        fs = FeatureSet.NewFromFitFile( wndchrm_test_dir + sep + 'test-l.fit', tile_options  )

        # What if the feature set number of groups within a class are less than called for
        # when specifying by integer?
      

if __name__ == '__main__':
    unittest.main()
