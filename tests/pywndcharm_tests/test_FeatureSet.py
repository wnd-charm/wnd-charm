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

from wndcharm.FeatureSet import FeatureSet, FeatureSet_Discrete, FisherFeatureWeights,\
        DiscreteBatchClassificationResult, FeatureSet_Continuous, ContinuousFeatureWeights,\
        ContinuousBatchClassificationResult

class TestFeatureSet( unittest.TestCase ):
    """
    Test the wndcharm module to check it still works with OmeroPychrm
    """

    def setUp(self):
        pass

    def test_ContinuousFitOnFit( self ):
        from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete

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

        from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Continuous

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

    def test_SplitOptions( self ):
        from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete

        fs_discrete = CreateArtificialFeatureSet_Discrete( n_samples=1000, n_classes=10,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True)


        from numpy.random import RandomState
        prng = RandomState(42)
        # default
        train_set, test_set = fs_discrete.Split( random_state=prng, quiet=True )
        self.assertEqual( train_set.shape, (750, 360) )
        self.assertEqual( test_set.shape, (250, 360) )

        # dummyproofing

        self.assertRaises( ValueError, fs_discrete.Split, train_size='trash' )
        self.assertRaises( ValueError, fs_discrete.Split, train_size=1.1 )
        self.assertRaises( ValueError, fs_discrete.Split, test_size='trash' )
        self.assertRaises( ValueError, fs_discrete.Split, test_size=1.1 )

        # What if the feature set number of groups within a class are less than called for
        # when specifying by integer?
        self.assertRaises( ValueError, test_set.Split, test_size=25 )


    def test_SampleReduce( self ):
        from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete

        n_classes = 10
        #========================================================
        # Section 1: LEAVE IN, Untiled Discrete (w/ classes) FeatureSets
        fs_discrete = CreateArtificialFeatureSet_Discrete( n_samples=1000, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True)

        # For discrete feature sets, passing in a flat list of ints ain't gonna cut it,
        # they have to be arranged into lists of lists indicating class composition
        desired = range(15)
        self.assertRaises( TypeError, fs_discrete.SampleReduce, desired )

        # Reduce to 9 classes from 10, one sample per class
        desired = [ [val] for val in xrange(50, 950, 100)]
        desired.insert( 2, False )
        A = fs_discrete.SampleReduce( desired )

        correct_samplenames = ['FakeClass-100_050', 'FakeClass-77.78_050', 'FakeClass-55.56_050', 'FakeClass-33.33_050', 'FakeClass-11.11_050', 'FakeClass11.11_050', 'FakeClass33.33_050', 'FakeClass55.56_050', 'FakeClass77.78_050']
        self.assertEqual( correct_samplenames, A._contiguous_samplenames_list )

        # Note the absense of the classname at index 2, "FakeClass-55.56"
        correct_classnames = ['FakeClass-100', 'FakeClass-77.78', 'FakeClass-33.33', 'FakeClass-11.11', 'FakeClass11.11', 'FakeClass33.33', 'FakeClass55.56', 'FakeClass77.78', 'FakeClass100']
        self.assertEqual( correct_samplenames, A._contiguous_samplenames_list )
        del A

        # Request more classes than exists in the original
        desired = [ [val] for val in xrange(50, 1000, 50)]
        B = fs_discrete.SampleReduce( desired )
        self.assertEqual( B.num_classes, len( desired ) )
        del B

        #========================================================
        # Section 2: LEAVE OUT, UNTiled Feature sets, Discrete FeatureSets

        # check that function barfs when passed an iterable containing anythin other than ints
        UNdesired = [ [val] for val in xrange(50, 950, 100)]
        self.assertRaises( TypeError, fs_discrete.SampleReduce,
                leave_out_samplegroupid_list=UNdesired )

        UNdesired = range(50, 950, 100)
        C = fs_discrete.SampleReduce( leave_out_samplegroupid_list=UNdesired )
        self.assertEqual( C.num_images, fs_discrete.num_images - len( UNdesired ) )

        # Single integers for leave_out_list is ok
        UNdesired = 50
        C = fs_discrete.SampleReduce( leave_out_samplegroupid_list=UNdesired )
        self.assertEqual( C.num_images, fs_discrete.num_images - 1 )
        del C

        #========================================================
        # Section 3: LEAVE IN, Tiled Feature sets, Discrete FeatureSets
        num_tiles = 4
        fs_discrete = CreateArtificialFeatureSet_Discrete( n_samples=1000, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=num_tiles, interpolatable=True)

        desired = [ [val] for val in xrange(5, 95, 10) ] # Rearrange into 9 classes
        D = fs_discrete.SampleReduce( desired )
        # Total num samples should be 9 classes, 1 sample group per class, 4 tiles per SG = 36
        self.assertEqual( num_tiles * len( desired ), D.num_images )
        del D

        #========================================================
        # Section 4: LEAVE OUT, WITH Tiled Feature sets, Discrete FeatureSets

        # check that function barfs when passed an iterable containing anythin other than ints
        UNdesired = [ [val] for val in xrange(50, 950, 100)]
        self.assertRaises( TypeError, fs_discrete.SampleReduce,
                leave_out_samplegroupid_list=UNdesired )

        # You can't leave out a sample group that doesn't exist
        UNdesired = range(50000, 50010)
        self.assertRaises( ValueError, fs_discrete.SampleReduce,
                leave_out_samplegroupid_list=UNdesired )

        # Can't leave out trash
        UNdesired = ['foo', 'bar']
        self.assertRaises( TypeError, fs_discrete.SampleReduce,
                leave_out_samplegroupid_list=UNdesired )

        # This input is ok:
        UNdesired = range(5, 95, 10)
        E = fs_discrete.SampleReduce( leave_out_samplegroupid_list=UNdesired )
        self.assertEqual( E.num_images, fs_discrete.num_images - len( UNdesired ) * num_tiles )
        del E

        #========================================================
        # Section 5: LEAVE IN, Untiled Continuous FeatureSets
        from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Continuous

        fs_cont = CreateArtificialFeatureSet_Continuous( n_samples=1000,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1)

        # dummyproof
        desired = ['foo', 'bar']
        self.assertRaises( TypeError, fs_cont.SampleReduce, desired )

        # flat lists only please
        desired = [ [val] for val in xrange(50, 950, 100)]
        self.assertRaises( TypeError, fs_cont.SampleReduce, desired )

        desired = range(50, 950)
        F = fs_cont.SampleReduce( desired )
        self.assertEqual( F.num_images, len(desired) )
        del F

        #========================================================
        # Section 6: LEAVE OUT, Untiled Continuous FeatureSets

        UNdesired = range(50, 950)
        G = fs_cont.SampleReduce( leave_out_samplegroupid_list=UNdesired )
        self.assertEqual( G.num_images, fs_cont.num_images - len(UNdesired) )
        del G

        # single int is ok
        H = fs_cont.SampleReduce( leave_out_samplegroupid_list=998 )
        self.assertEqual( H.num_images, fs_cont.num_images - 1 )
        del H

        #========================================================
        # Section 7: LEAVE IN, TILED Continuous FeatureSets

        fs_cont = CreateArtificialFeatureSet_Continuous( n_samples=1000,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=num_tiles)


        desired = range(50, 95)
        I = fs_cont.SampleReduce( desired )
        self.assertEqual( I.num_images, len(desired) * num_tiles )
        del I

        # single int is ok, ALTHOUGH NOT SURE WHY YOU'D EVER WANT A FS WITH A SINGLE SAMPLE
        J = fs_cont.SampleReduce( 98 )
        self.assertEqual( J.num_images, num_tiles )
        del J


        #========================================================
        # Section 8: LEAVE OUT, TILED Continuous FeatureSets

        UNdesired = range(50, 95)
        K = fs_cont.SampleReduce( leave_out_samplegroupid_list=UNdesired )
        self.assertEqual( K.num_images, fs_cont.num_images - len(UNdesired) * num_tiles )
        del K

        # single int is ok
        L = fs_cont.SampleReduce( leave_out_samplegroupid_list=98 )
        self.assertEqual( L.num_images, fs_cont.num_images - num_tiles  )
        del L

if __name__ == '__main__':
    unittest.main()
