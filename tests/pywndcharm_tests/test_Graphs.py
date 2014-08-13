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

from os.path import sep, dirname, realpath
import filecmp
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

from wndcharm.FeatureSet import FisherFeatureWeights, DiscreteBatchClassificationResult,\
        PredictedValuesGraph, DiscreteClassificationExperimentResult

from wndcharm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete

try:
    import matplotlib
    HasMatplotlib = True
except ImportError:
    HasMatplotlib = False

class TestGraphs( unittest.TestCase ):
    """Test WND-CHARM's graph-making functionality."""
    
    fs = CreateArtificialFeatureSet_Discrete( n_samples=1000, n_classes=10, 
            initial_noise_sigma=100, noise_gradient=10, random_state=43 )
    train_set, test_set = fs.Split( randomize=False, quiet=True )
    train_set.Normalize( quiet=True )
    fw = FisherFeatureWeights.NewFromFeatureSet( train_set ).Threshold()

    reduced_train_set = train_set.FeatureReduce( fw.names )
    reduced_test_set = test_set.FeatureReduce( fw.names )
    reduced_test_set.Normalize( reduced_train_set, quiet=True )

    batch_result = DiscreteBatchClassificationResult.New(
            reduced_train_set, reduced_test_set, fw, quiet=True )

    def CompareGraphs( self, graph, testfilename ):
        """Helper function to check output graphs"""
        tempdir = mkdtemp()
        tempfile = tempdir + sep + testfilename
        try:
            graph.SaveToFile( tempfile )
            self.assertTrue( filecmp.cmp( testfilename, tempfile ) ) 
        finally:
            rmtree( tempdir )

    @unittest.skipIf( HasMatplotlib, "Skipped if matplotlib IS installed" )
    def test_ErrMsgIfMatplotibNotInstalled( self ):
        """Fail gracefully with informative message if matplotlib"""

        graph = PredictedValuesGraph( self.batch_result )
        with self.assertRaises( ImportError ):
            graph.RankOrderedPredictedValuesGraph()
        with self.assertRaises( ImportError ):
            graph.KernelSmoothedDensityGraph()

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_RankOrderedFromBatchClassificationResult( self ):
        """Rank Ordered Predicted values graph from a single split"""

        testfilename = 'test_graph_rank_ordered.png'
        graph = PredictedValuesGraph( self.batch_result )
        graph.RankOrderedPredictedValuesGraph()
        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_KernelSmoothedFromBatchClassificationResult( self ):
        """Kernel Smoothed Probability density graph from a single split"""

        testfilename = 'test_graph_kernel_smoothed.png'
        graph = PredictedValuesGraph( self.batch_result )
        graph.KernelSmoothedDensityGraph()
        self.CompareGraphs( graph, testfilename )

    @unittest.skip( "Skip until ShuffleSplit has a RandomState param to pass in" )
    def test_FromDiscreteClassificationExperimentResults( self ):
        """Rank Ordered Predicted values graph from an experiment result (multiple splits)"""

        testfilename = 'test_graph_rank_ordered_experiment.png'

        small_fs = CreateArtificialFeatureSet_Discrete( n_samples=100, n_classes=5, 
                initial_noise_sigma=100, noise_gradient=10, random_state=42 )
        experiment = DiscreteClassificationExperimentResult.NewShuffleSplit( small_fs, quiet=True )
        graph = PredictedValuesGraph( self.batch_result )
        graph.RankOrderedPredictedValuesGraph()
        # graph.SaveToFile( testfilename ) # remove after RandomState for ShuffleSplit is implemented
        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_FromHTML( self ):
        """Rank Ordered Predicted values graph from an experiment result (multiple splits)"""

        testfilename = 'test_graph_fromHTML.png' 
         # Inflate the zipped html file into a temp file
        import zipfile
        zipped_file_path = pychrm_test_dir + sep + 'c_elegans_terminal_bulb.html.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf.extractall( tempdir )
        tempfile = tempdir + sep + testfilename
        try:
            htmlfilepath = tempdir + sep + zf.namelist()[0]
            _exp_result = DiscreteClassificationExperimentResult.NewFromHTMLReport( htmlfilepath )
            graph = PredictedValuesGraph( _exp_result )
            graph.RankOrderedPredictedValuesGraph()
            graph.SaveToFile( testfilename )
            graph.SaveToFile( tempfile )
            self.assertTrue( filecmp.cmp( testfilename, tempfile ) ) 
        finally:
            rmtree( tempdir )


    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOTinstalled" )
    def test_IfNotInterpolatable( self ):
        """You can't graph predicted values if the classes aren't interpolatable."""

        testfilename = "noop" 
        non_interp_fs = CreateArtificialFeatureSet_Discrete( n_samples=100, n_classes=5, 
                initial_noise_sigma=100, noise_gradient=10, random_state=42, interpolatable=False )
        _train_set, _test_set = non_interp_fs.Split( randomize=False, quiet=True )
        _train_set.Normalize( quiet=True )
        _fw = FisherFeatureWeights.NewFromFeatureSet( _train_set ).Threshold()

        _reduced_train_set = _train_set.FeatureReduce( _fw.names )
        _reduced_test_set = _test_set.FeatureReduce( _fw.names )
        _reduced_test_set.Normalize( _reduced_train_set, quiet=True )

        _batch_result = DiscreteBatchClassificationResult.New( _reduced_train_set, _reduced_test_set, _fw, quiet=True )
        graph = PredictedValuesGraph( _batch_result )
        with self.assertRaises( ValueError ):
            graph.RankOrderedPredictedValuesGraph()


if __name__ == '__main__':
    unittest.main()
