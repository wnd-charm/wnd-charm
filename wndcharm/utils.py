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


# wndcharm.py has the definitions of all the SWIG-wrapped primitive C++ WND_CHARM objects.
import wndcharm
import numpy as np

# ============================================================
# BEGIN: Initialize module level globals
Algorithms = []
Transforms = []

def initialize_module(): 
    """If you're going to calculate any features, you need this stuff.
    FIXME: Rig this stuff to load only on demand."""
    
    global Algorithms
    global Transforms
    # The verbosity is set by the environment variable WNDCHRM_VERBOSITY, in wndchrm_error.cpp
    # wndcharm.cvar.verbosity = 7

    # These are the auto-registered ComputationTasks, which come in different flavors
    all_tasks = wndcharm.ComputationTaskInstances.getInstances()
    for task in all_tasks:
        if task.type == task.ImageTransformTask:
            Transforms.append (task)
#            print task.name + " added to Transforms"
        elif task.type == task.FeatureAlgorithmTask:
            Algorithms.append (task)
#            print task.name + " added to Algorithms"

    # The standard feature plans get loaded as needed from C++ statics, so they don't need to be initialized.
    # Standard sets (other plan "parts" are in Tasks.h under StdFeatureComputationPlans)
    #     getFeatureSet();
    #     getFeatureSetColor();
    #     getFeatureSetLong();
    #     getFeatureSetLongColor();

    # e. g.:
#     small_feature_plan = wndcharm.StdFeatureComputationPlans.getFeatureSet()
#     print "small feature set groups:"
#     last_feature_group = None;
#     for i in range( 0, small_feature_plan.n_features):
#         feature_group = small_feature_plan.getFeatureGroupByIndex(i)
#         feature_name = small_feature_plan.getFeatureNameByIndex(i)
#         if feature_group.name != last_feature_group:
#             print "feature_group "+feature_group.name
#         last_feature_group = feature_group.name
#         print "  feature_name "+feature_name

    # while we're debugging, raise exceptions for numerical weirdness, since it all has to be dealt with somehow
    # In cases where numerical weirdness is expected and dealt with explicitly, these exceptions are
    # temporarily turned off and then restored to their previous settings.
    np.seterr (all='raise')

# ============================================================
def output_railroad_switch( method_that_prints_output ):
    """This is a decorator that optionally lets the user specify a file to which to redirect
    STDOUT. To use, you must use the keyword argument "output_filepath" and optionally
    the keyword argument "mode" """

    def print_method_wrapper( *args, **kwargs ):
        
        retval = None
        if "output_filepath" in kwargs:
            output_filepath = kwargs[ "output_filepath" ]
            del kwargs[ "output_filepath" ]
            if "mode" in kwargs:
                mode = kwargs[ "mode" ]
                del kwargs[ "mode" ]
            else:
                mode = 'w'
            print 'Saving output of function "{0}()" to file "{1}", mode "{2}"'.format(\
                  method_that_prints_output.__name__, output_filepath, mode )
            import sys
            backup = sys.stdout
            sys.stdout = open( output_filepath, mode )
            retval = method_that_prints_output( *args, **kwargs )
            sys.stdout.close()
            sys.stdout = backup
        elif "output_stream" in kwargs:
            output_stream = kwargs[ "output_stream" ]
            del kwargs[ "output_stream" ]
            print 'Saving output of function "{0}()" to stream'.format(\
                  method_that_prints_output.__name__)
            import sys
            backup = sys.stdout
            try:
                sys.stdout = output_stream
                retval = method_that_prints_output( *args, **kwargs )
            finally:
                sys.stdout = backup
        else:
            retval = method_that_prints_output( *args, **kwargs )
        return retval

    return print_method_wrapper

# ============================================================
def ReplaceNonReal( feature_matrix, mins=None, maxs=None ):
    """Useful for cross-validation, so you don't have to check for NaNs/INFs every split.
    Will clip to mins maxs if passed.

    Args:
        feature_matrix - numpy.ndarray
            matrix shape = NxM, where N = num samples, M = num features
        mins - numpy.ndarray, shape 1xM
        maxs - numpy.ndarray, shape 1xM
            use these to clip/set in max if PINF/NINF exists
    Returns:
        False if no modifications, True if modified

    Will clip to mins/maxs if passed as aruguments.
    NANs in the columns will be set to 0."""

    # Edge cases to deal with:
    #   Range determination:
    #     1. features that are nan, inf, -inf
    #        max and min determination must ignore invalid numbers
    #        nan -> 0, inf -> max, -inf -> min
    #   Normalization:
    #     2. feature values outside of range
    #        values clipped to range (-inf to min -> min, max to inf -> max) - leaves nan as nan
    #     3. feature ranges that are 0 result in nan feature values
    #     4. all nan feature values set to 0

    feature_matrix_m = np.ma.masked_invalid( feature_matrix, copy=False )
    # Masked cells are True, unmasked are False
    if not np.any( feature_matrix_m.mask ):
        # Nothing to do, no masked cells.
        return False

    # First take care of INFs:
    # NANs and +/-INFs have been masked above to facilitate computation of min/max
    if mins is None:
        maxs = feature_matrix_m.max( axis=0 )
    if maxs is None:
        mins = feature_matrix_m.min( axis=0 )
    # clip the values to the min-max range (NANs are left, but +/- INFs are taken care of)
    feature_matrix.clip( mins, maxs, feature_matrix )

    # Finally take care of NaNs
    feature_matrix[ np.isnan( feature_matrix ) ] = 0
    return True

# ============================================================
def normalize_by_columns( feature_matrix, mins=None, maxs=None, means=None, stdevs=None,
        zscore=False, non_real_check=True ):
    """Performs feature scaling/normalization in-place on input arg feature_matrix.

    Args:
        feature_matrix - numpy.ndarray
            matrix shape = NxM, where N = num samples, M = num features
        mins - numpy.ndarray, shape 1xM
        maxs - numpy.ndarray, shape 1xM
        means - numpy.ndarray, shape 1xM
        stdevs - numpy.ndarray, shape 1xM
            Reference vales to transform this feature space if passed, Must not contain NaNs/INFs
        zscore - bool
            If False, min max scaling to interval 0-100, else Z-score standardization
        non_real_check - bool
            Comb through matrix and replace non-real numbers in accordance with prior
            WND-CHARM convention (calls ReplaceNonReal)
    Returns:
        (mins, maxs, means, stdevs) - 2-tuple of numpy.ndarrays w/ shape 1xM
            Vals may be None if it wasn't pased, derived as part of running this function.

    Notes:
        zero-range columns will be set to 0.
        The normalized output range is hard-coded to 0-100."""

    # Step 1: if requested, replace NaN's/INFs
    if non_real_check:
        replaced = ReplaceNonReal( feature_matrix, mins, maxs )
        if (mins is not None or maxs is not None) and replaced == False:
            # This means there was no NaNs/INFs to replace, but feature_matrix still
            # needs to be clipped.
            feature_matrix.clip( mins, maxs, feature_matrix )
    # Step 2: make sure feature space has been clipped if it wasn't done in previous step.
    elif mins is not None or maxs is not None:
        feature_matrix.clip( mins, maxs, feature_matrix )

    # Turn off numpy warnings, since we're taking care of invalid values explicitly
    oldsettings = np.seterr( all='ignore' )

    if not zscore:
        # Do min-max scaling to interval 0-100
        if mins is None:
            mins = feature_matrix.min( axis=0 )
        if maxs is None:
            maxs = feature_matrix.max( axis=0 )
        feature_matrix -= mins
        feature_matrix /= (maxs - mins)
        feature_matrix *= 100
        nan_cols = (maxs - mins) == 0
    else:
        # Perform z-score normalization on feature space
        if means is None:
            means = feature_matrix.mean( axis=0 )
        if stdevs is None:
            stdevs = feature_matrix.std( axis=0 )
        feature_matrix -= means
        feature_matrix /= stdevs
        nan_cols = stdevs == 0

    # Selectively fill NaN cols created when dividing by 0
    if len( feature_matrix.shape ) > 1:
        feature_matrix[:, nan_cols ] = 0
    else:
        feature_matrix[ nan_cols ] = 0

    # return settings to original
    np.seterr( **oldsettings )

    return( mins, maxs, means, stdevs )

# END: Initialize module level globals
#===============================================================

initialize_module()

# BEGIN: Helper functions
def compare( a_list, b_list, atol=1e-7 ):
    """Helps to compare floating point values to values stored
    in text files (ala .fit and .sig files) where the number of
    significant figures is sometimes orders of magnitude different"""

    result = True
    errcount = 0
    for count, (a_raw, b_raw) in enumerate( zip( a_list, b_list ) ):
        if errcount > 20:
            break

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
            print errmsg
            result = False
            errcount += 1
            continue
        if e_in_a_str:
            a_coeff, a_exp = a_str.split( 'e' )
            b_coeff, b_exp = b_str.split( 'e' )
            if a_exp != b_exp:
                # AssertionError: Index 623: "1e-06" and "6.93497e-07" exponents don't match.
                a_exp = int( a_exp )
                b_exp = int( b_exp )
                exp_digits = abs( a_exp - b_exp )

                #errmsg = "Index {0}: \"{1}\" and \"{2}\" exponents don't match."
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
            errstr = "****Index {0}: {1} isn't enough like {2}".format( count, a_raw, b_raw )
            print errstr
            result = False
            errcount += 1
            continue
            #self.fail( errstr.format( count, a_raw, b_raw ) )
    return result

# ============================================================

def RunInProcess( fv ):
    #print "**DEBUG fv {} row-{} col{}".format( fv, fv.tile_row_index, fv.tile_col_index )
    fv.GenerateFeatures( write_to_disk=True )
    return True

def parallel_compute( samples, n_jobs=True ):
    """WND-CHARM implementation of symmetric multiprocessing, see:
    https://en.wikipedia.org/wiki/Symmetric_multiprocessing"""

    from multiprocessing import cpu_count, Queue, Pool, log_to_stderr
    import logging
    if n_jobs == True:
        n_jobs = cpu_count()

    logger = log_to_stderr()
    logger.setLevel(logging.INFO)
    try:
        pool = Pool( processes=n_jobs )
        pool.map( RunInProcess, samples, chunksize=1 )
        pool.close()
        pool.join()
    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        raise
 
