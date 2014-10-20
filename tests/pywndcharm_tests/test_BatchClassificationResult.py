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

from wndcharm.FeatureSet import FeatureSet_Discrete, FisherFeatureWeights,\
        DiscreteBatchClassificationResult, FeatureSet_Continuous, ContinuousFeatureWeights,\
        ContinuousBatchClassificationResult

class TestDiscreteBatchClassificationResult( unittest.TestCase ):
    """
    Test the classification functionality
    """

    def setUp( self ):
        pass
    

    
    def test_TiledFitOnFit( self ):
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
            #fs.Print( verbose=True )
            #print "\n\n\n********************\n\n\n"
            #full_train, full_test = fs.Split( random_state=42, quiet=True )
            #full_train.Print( verbose=True )
            #full_test.Print( verbose=True )
            fs.Normalize( quiet=True )
            fw = FisherFeatureWeights.NewFromFeatureSet( fs ).Threshold()
            reduced_fs = fs.FeatureReduce( fw.names )
            #reduced_test = full_test.FeatureReduce( fw.names )
            #reduced_test.Normalize( reduced_train, quiet=True )

            batch_result = DiscreteBatchClassificationResult.New( reduced_fs, reduced_fs, fw )
            batch_result.Print()

        finally:
            rmtree( tempdir )

            

if __name__ == '__main__':
    unittest.main()
