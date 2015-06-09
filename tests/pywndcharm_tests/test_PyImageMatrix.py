"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2015 National Institutes of Health 

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

from wndcharm.PyImageMatrix import PyImageMatrix

#from wndcharm.utils import compare

from os.path import dirname, sep, realpath, join, abspath, splitext, basename
from tempfile import mkdtemp
from shutil import rmtree
from numpy.testing import assert_equal

import numpy as np
import matplotlib.pyplot as plt

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )

test_dir = wndchrm_test_dir

class TestPyImageMatrix( unittest.TestCase ):
    """WND-CHARM object that handles pixels and feature calculation"""

    def test_OpenImage( self ):
        """Testing ROI open functionality"""

        tempdir = mkdtemp()
        # crop the test image (size=1388x1040) to be bottom right tile of the test image
        # tiled with 6 cols and 5 rows => ROI= 231x208+1155+832

        orig_big = join( pychrm_test_dir, 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif' )
        orig_cropped = join( pychrm_test_dir, 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E-t6x5_5_4-l.tiff' )
        test_cropped_path = join( tempdir, 'TEST_OpenImageROI.tif' )

        x = 1155
        y = 832
        w = 231
        h = 208
       
        from wndcharm import rect
        bb = rect()
        bb.x = x
        bb.y = y
        bb.w = w
        bb.h = h

        # setting mean = 0 is flag to not use mean in ImageMatrix::OpenImage()
        downsample = 0
        #bb = None
        mean = 0.0
        stddev = 0.0 
        try:
            cropped_on_load_im = PyImageMatrix()

            if 1 != cropped_on_load_im.OpenImage( orig_big, downsample, bb, mean, stddev ):
                self.fail( 'Could not build an ImageMatrix from ' + orig_big )

            cropped_on_load_im.SaveTiff( test_cropped_path )

            orig_pixels = plt.imread( orig_cropped )
            cropped_pixels = plt.imread( test_cropped_path )

            assert_equal( orig_pixels, cropped_pixels )

        finally:
            rmtree( tempdir )

    def test_SaveImage( self ):
        """SaveImage pixels to file."""

        tempdir = mkdtemp()
        orig = join( pychrm_test_dir, 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif' )

        # setting mean = 0 is flag to not use mean in ImageMatrix::OpenImage()
        downsample = 0
        bb = None
        mean = 0.0
        stddev = 0.0 
        try:
            origim = PyImageMatrix()
            if 1 != origim.OpenImage( orig, downsample, bb, mean, stddev ):
                self.fail( 'Could not build an ImageMatrix from ' + orig )

            copied_tiff = join( tempdir, 'TEST1.tif') 
            origim.SaveTiff( copied_tiff )
            
            orig_pixels = plt.imread( orig )
            copied_pixels = plt.imread( copied_tiff )

            assert_equal( orig_pixels, copied_pixels )

        finally:
            rmtree( tempdir )

    def test_allocate( self ):
        """make an empty pixel plane"""

        tempdir = mkdtemp()
        test_path = join( tempdir, "TEST_allocate.tif" )

        num_rows, num_cols = shape = (123,456)
        
        try:
            origim = PyImageMatrix()
	    # virtual void allocate (unsigned int w, unsigned int h);
            origim.allocate( num_cols, num_rows  )
            origim.SaveTiff( test_path )
            pixels = plt.imread( test_path )
            self.assertEqual( pixels.shape, shape )

        finally:
            rmtree( tempdir )

    def test_submatrix( self ):
        """crop"""

        tempdir = mkdtemp()

        # crop the test image (size=1388x1040) to be bottom right tile of the test image
        # tiled with 6 cols and 5 rows => ROI= 231x208+1155+832
        orig_big = join( pychrm_test_dir, 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif' )
        orig_cropped = join( pychrm_test_dir, 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E-t6x5_5_4-l.tiff' )
        test_cropped_path = join( tempdir, 'TEST_submatrix.tif' )
        x = 1155
        y = 832
        w = 231
        h = 208

        x1 = x
        y1 = y
        x2 = x1 + w - 1
        y2 = y1 + h - 1
        
        # setting mean = 0 is flag to not use mean in ImageMatrix::OpenImage()
        downsample = 0
        bb = None
        mean = 0.0
        stddev = 0.0 
        try:
            origim = PyImageMatrix()
            if 1 != origim.OpenImage( orig_big, downsample, bb, mean, stddev ):
                self.fail( 'Could not build an ImageMatrix from ' + orig_big )

            # API calls for copying desired pixels into empty ImageMatrix instance:
            # the_tiff is garbage collected on return
            cropped_im = PyImageMatrix()
            # void ImageMatrix::submatrix (const ImageMatrix &matrix, const unsigned int x1, const unsigned int y1, const unsigned int x2, const unsigned int y2)
            cropped_im.submatrix( origim, x1, y1, x2, y2 ) # no retval

            cropped_im.SaveTiff( test_cropped_path )

            orig_pixels = plt.imread( orig_cropped )
            cropped_pixels = plt.imread( test_cropped_path )

            assert_equal( orig_pixels, cropped_pixels )

        finally:
            rmtree( tempdir )


if __name__ == '__main__':
    unittest.main()
