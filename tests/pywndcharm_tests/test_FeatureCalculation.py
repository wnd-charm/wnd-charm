#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2014 National Institutes of Health
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

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

from wndcharm.FeatureSet import FeatureSpace, FisherFeatureWeights,\
        DiscreteImageClassificationResult, FeatureVector

from os.path import dirname, realpath, join

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )

test_dir = wndchrm_test_dir

class TestFeatureCalculation( unittest.TestCase ):
    """Feature Calculation"""

    sig_file_path = join( test_dir,'010067_301x300-l_precalculated.sig' )
    test_tif_path = join( test_dir,'010067_301x300.tif' )

    # An alternate target for sig calculation:
    # sig_file = os.path.join (test_dir,'t1_s01_c05_ij-l_precalculated.sig')
    # test_tif = os.path.join (test_dir,'t1_s01_c05_ij.tif')

    def compare( self, a_list, b_list, atol=1e-7 ):

        for count, (a_raw, b_raw) in enumerate( zip( a_list, b_list ) ):

            if a_raw == b_raw:
                continue
            if abs( float( a_raw ) - float( b_raw ) ) < atol:
                continue

            a_str = "{0:0.6g}".format( a_raw )
            b_str = "{0:0.6g}".format( b_raw )

            # These deal with the 1e-6 ~ 6.93e-7 comparison issue
            # a_addl_zero = ""
            # b_addl_zero = ""

            exp_digits = 0
            e_in_a_str = 'e' in a_str
            e_in_b_str = 'e' in b_str
            if e_in_a_str != e_in_b_str:
                errmsg = "Index {0}: \"{1}\" and \"{2}\" exponents don't match."
                self.fail( errmsg.format( count, a_str, b_str, ) )
            if e_in_a_str:
                a_coeff, a_exp = a_str.split( 'e' )
                b_coeff, b_exp = b_str.split( 'e' )
                if a_exp != b_exp:
                    # AssertionError: Index 623: "1e-06" and "6.93497e-07" exponents don't match.
                    a_exp = int( a_exp )
                    b_exp = int( b_exp )
#                    if a_exp > b_exp:
#                        a_addl_zero = '0'* abs( a_exp - b_exp )
#                    else:
#                        b_addl_zero = '0'* abs( a_exp - b_exp )
                    exp_digits = abs( a_exp - b_exp )

                    #errmsg = "Index {0}: \"{1}\" and \"{2}\" exponents don't match."
                    #self.fail( errmsg.format( count, a_raw, b_raw, ) )
                # FIXME: lstrip doesn't properly deal with negative numbers
                a_int_str = a_coeff.translate( None, '.' ).lstrip('0')
                b_int_str = b_coeff.translate( None, '.' ).lstrip('0')
            else:
                a_int_str = a_str.translate( None, '.' ).lstrip('0')
                b_int_str = b_str.translate( None, '.' ).lstrip('0')

            a_len = len( a_int_str )
            b_len = len( b_int_str )
            diff_digits = abs( a_len - b_len ) + exp_digits
            tail = '0' * diff_digits

            #msg = 'a_str "{0}" (len={1}), b_str "{2}" (len={3}), tail="{4}"'
            #print msg.format( a_int_str, a_len, b_int_str, b_len, tail )

            if a_len > b_len:
#                a = int( a_int_str + a_addl_zero )
#                b = int( b_int_str + tail + b_addl_zero )
#            else:
#                a = int( a_int_str + tail + a_addl_zero )
#                b = int( b_int_str + b_addl_zero )
                a = int( a_int_str )
                b = int( b_int_str + tail )
            elif b_len > a_len:
                a = int( a_int_str + tail )
                b = int( b_int_str )
            else:
                a = int( a_int_str )
                b = int( b_int_str )
            # Rounding is useless, since due to floating point's inexact representation
            # it's possible to have 2.65 round to 2.6
            #a = round( a, -1 * diff_digits )
            #b = round( b, -1 * diff_digits )

            diff = abs( a - b )

            #print "{0}->{1}=={2}<-{3} : {4} <= {5}".format( a_raw, a, b, b_raw, diff, 10 ** diff_digits )
            if diff > 10 ** diff_digits:      
                errstr = "Index {0}: {1} isn't enough like {2}"
                self.fail( errstr.format( count, a_raw, b_raw ) )

    # --------------------------------------------------------------------------
    @unittest.skip( "Not doing anything with this right now" )
    def test_ProfileLargeFeatureSet( self ):
        """Profiling for calculating sigs"""

        import cProfile
        import tempfile
        import pstats
        prof = tempfile.NamedTemporaryFile()
        cmd = 'FeatureVector( source_filepath="{0}", long=True ).GenerateFeatures()'.format( self.test_tif_path )
        cProfile.run( cmd, prof.name, 'time')
        p = pstats.Stats(prof.name)
        p.sort_stats('time').print_stats(5)
        prof.close()

#Loaded features from file /Users/chris/src/wnd-charm/tests/wndchrm_tests/010067_301x300-l_precalculated.sig
#.Fri Jan 23 15:15:35 2015    /var/folders/cr/vsd9_15x6xbc3np6rvx12mqm0000gp/T/tmpf7Oof0
#
#         17595 function calls in 14.956 seconds
#
#   Ordered by: internal time
#   List reduced from 55 to 5 due to restriction <5>
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#        1   14.926   14.926   14.926   14.926 {_wndcharm.FeatureComputationPlanExecutor_run}
#        1    0.005    0.005   14.956   14.956 <string>:1(<module>)
#     2920    0.004    0.000    0.004    0.000 {_wndcharm.SwigPyIterator_next}
#     2922    0.004    0.000    0.004    0.000 {method 'format' of 'str' objects}
#        1    0.003    0.003   14.951   14.951 /Users/chris/src/wnd-charm/build/lib.macosx-10.9-x86_64-2.7/wndcharm/FeatureSet.py:1126(GenerateFeatures)
        # FIXME: Actually do some checking of the profile results

    # --------------------------------------------------------------------------
    def test_LargeFeatureSetGrayscale( self ):
        """Large feature set, grayscale image"""
        reference_sample = FeatureVector.NewFromSigFile( self.sig_file_path,
            image_path=self.test_tif_path )

        target_sample = FeatureVector( source_filepath=self.test_tif_path,
            long=True).GenerateFeatures( write_sig_files_to_disk=False )

#        This doesn't work since the ranges of features are so wide
#        Tried using relative tolerance, but no dice:
#        from numpy.testing import assert_allclose
#        assert_allclose( reference_sample.values, target_sample.values, rtol=1e-3 )

        # Remember we're reading these values in from strings. and the ranges are so wide
        # you only have 6 sig figs. Better apples to apples comparison is to 
        # compare strings.

        self.compare( target_sample.values, reference_sample.values )


if __name__ == '__main__':
    unittest.main()
