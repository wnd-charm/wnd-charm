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

import numpy as np

lbound = -100
ubound = 100

# signal types (picked to have good range of output for inputs [-100,100])
# inputs are np arrays 1XD, where D is num_classes or num_samples

# No singularities, functions defined and differentiable for input of all Real numbers
well_behaved_signals = {
'positiveconstant'     : lambda x: 100 * np.ones( x.shape ),
'negativeconstant'     : lambda x: -100 * np.ones( x.shape ),
'positive_linear'      : lambda x: np.copy( x ),
'negative_linear'      : lambda x: -x,
'trough_first_sin'     : lambda x: 100 * np.sin( x * np.pi / 40 ),
'peak_first_sin'       : lambda x: -100 * np.sin( x * np.pi / 40 ),
'exponential_quadI'    : lambda x: np.exp( 0.05 * x ),
'exponential_quadII'   : lambda x: np.exp( -0.05 * x ),
'exponential_quadIII'  : lambda x: -np.exp( -0.05 * x ),
'exponential_quadIV'   : lambda x: -np.exp( 0.05 * x ),
'quadraticup'          : lambda x: 0.01 * np.square( x ) - 50,
'quadraticdown'        : lambda x: -0.01 * np.square( x ) + 50,
'cubicup'              : lambda x: np.power( x, 3 ) / 4500 - x + 50,
'cubicdown'            : lambda x: -np.power( x, 3 ) / 4500 + x - 50,
'arctan_left_pos'      : lambda x: 50 * np.arctan( 0.1 * x + 10 ),
'arctan_middle_pos'    : lambda x: 50 * np.arctan( 0.1 * x ),
'arctan_right_pos'     : lambda x: 50 * np.arctan( 0.1 * x - 10 ),
'arctan_left_neg'      : lambda x: -50 * np.arctan( 0.1 * x + 10 ),
'arctan_middle_neg'    : lambda x: -50 * np.arctan( 0.1 * x ),
'arctan_right_neg'     : lambda x: -50 * np.arctan( 0.1 * x - 10 ),
}


singularity_signals = {
'log_quadI'            : lambda x: 100 * np.log10( 0.05 * x ),
'log_quadII'           : lambda x: 100 * np.log10( -0.05 * x ),
'log_quadIII'          : lambda x: -100 * np.log10( -0.05 * x ),
'log_quadIV'           : lambda x: -100 * np.log10( 0.05 * x ),
'stepup'               : lambda x: 100 * np.sign( x ),
'stepdown'             : lambda x: -100 * np.sign( x ),
'triangleup'           : lambda x: np.fabs( x ) - 50,
'triangledown'         : lambda x: -np.fabs( x ) + 50,
'tangent1'             : lambda x: 5 * np.tan( np.pi * x / 50 ),
'tangent2'             : lambda x: -5 * np.tan( np.pi * (x / 50 - 0.5 ) ),
}

from .FeatureSpace import FeatureSpace

def CreateArtificialFeatureSpace_Continuous( name="ContinuousArtificialFS", n_samples=100,
    num_features_per_signal_type=25, noise_gradient=10, initial_noise_sigma=10,
    n_samples_per_group=1, random_state=None, singularity=None, clip=None ):
    """
    Analogous to sklearn.datasets.make_regression, but simplified.

    Returns an instance of a FeatureSet_Continuous.

    The number of features contained in the feature set will be a multiple of the number of
    signals contained in the dict "signals" (len(signals) == 12 at the time of writing, so
    len( new_fs.feature_names ) == 12, 24, 36, etc ).

    The more features the user asks for via the argument "num_features_per_signal_type",
    the more each successive feature generated using that signal will have a greater degree of
    gaussian noise added to it (controlled via args "noise_gradient" and "initial_noise_sigma").

    singularity = If evals to True, singularity signals are used
    clip = {False, True, (user_lbound, user_ubound)} signal clipping not used,
           used with default lbound and ubound, or used with user specifications
    """
    if n_samples_per_group < 1 or type(n_samples_per_group) is not int:
      raise ValueError( "n_samples_per_group has to be an integer, and at least 1." )

    if random_state:
        from numpy.random import RandomState
        if random_state is True:
            from numpy.random import normal
        elif type( random_state ) is RandomState:
            normal = random_state.normal
        elif type( random_state ) is int:
            normal = RandomState( random_state ).normal
        else:
            raise ValueError( 'Arg random_state must be an instance of np.random.RandomState, an int, or the value True')
    else: # no noise added to feature values
        normal = lambda mu, sigma, n: np.zeros( n )

    # Figure out what signals to use
    if singularity:
        # combine dicts of signals
        from itertools import chain
        signals = dict(chain(well_behaved_signals.iteritems(), singularity_signals.iteritems()))
    else:
        signals = well_behaved_signals

    # Use signal clipping?
    try:
        if len( clip ) == 2:
            hi, lo= float( clip[0] ), float( clip[1] )
            from functools import partial
            clip = partial( np.clip, a_min=lo, a_max=hi )
        else:
            raise ValueError( "Arg clip requires True, false, or a 2-tuple of numbers" )
    except TypeError:
        if clip:
            from functools import partial
            clip = partial( np.clip, a_min=lbound, a_max=ubound )
        else:
            clip = lambda x: x

    num_features = num_features_per_signal_type * len( signals )

    if n_samples_per_group > 1:
      # make n_samples evenly divisible
      n_samples = int( n_samples // n_samples_per_group ) * n_samples_per_group

    # Instantiate and assign basic data members
    new_fs = FeatureSpace( name=name, source_filepath=name, num_samples=n_samples,
      num_samples_per_group=n_samples_per_group, num_features=num_features, discrete=False,
      feature_set_version='-1.0')

    # The function call np.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_samples)
    new_fs._contiguous_ground_truth_values = list( np.mgrid[ lbound : ubound : step ] )

    # Generate artificial feature names
    # N.B. The feature generation signals are sorted in alphanum order!
    new_fs.feature_names = [ "{0}{1:04d}".format(fname, i)\
                                    for fname in sorted( signals.keys() ) \
                                    for i in xrange( num_features_per_signal_type ) ]
    # Creating sample metadata
    if n_samples_per_group == 1:
        new_fs._contiguous_sample_names = [ "FakeContinuousSample{0:03d}".format( i )\
                                             for i in xrange( n_samples ) ]
        new_fs._contiguous_sample_group_ids = range( n_samples ) # not xrange
        new_fs._contiguous_sample_sequence_ids =  [1] * n_samples
    else:
        # Format: FakeContinuousSample_i<sample index>_g<group>_t<tile#>
        temp1 = [ "FakeContinuousSample_i{0:03d}".format( i ) for i in xrange( n_samples ) ]
        n_samplegroups = n_samples / n_samples_per_group
        temp2 = [ "_g{0:03d}_t{1:02d}".format( samplegroup_index, tile_index ) \
            for samplegroup_index in xrange( n_samplegroups ) \
                for tile_index in xrange( n_samples_per_group ) ]
        new_fs._contiguous_sample_names = [ a + b for a, b in zip( temp1, temp2 ) ]
        new_fs._contiguous_sample_group_ids = \
            [ samplegroup_index for samplegroup_index in xrange( n_samplegroups ) \
                for tile_index in xrange( n_samples_per_group ) ]
        new_fs._contiguous_sample_sequence_ids = \
            [ tile_index for samplegroup_index in xrange( n_samplegroups ) \
                for tile_index in xrange( n_samples_per_group ) ]

    old_settings = np.seterr( all='ignore' ) # Bring on the NaN's!

    # Create features across all classes at the same time
    feat_count = 0
    # N.B. The features are in sort order!
    ground_truth_values = np.array( new_fs._contiguous_ground_truth_values )
    for func_name in sorted( signals.keys() ):
        f = signals[ func_name ]
        raw_feature_values = clip( f( ground_truth_values ) )
        for feat_index in xrange( num_features_per_signal_type ):
            # Add noise proportional to the feature index
            noise_vector = normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples )
            new_fs.data_matrix[:,feat_count] = np.add( noise_vector, raw_feature_values )
            feat_count += 1

    np.seterr( **old_settings )

    new_fs.samples_sorted_by_ground_truth = True
    new_fs._RebuildViews( recalculate_class_metadata=False )
    return new_fs

def CreateArtificialFeatureSpace_Discrete( name="DiscreteArtificialFS", n_samples=100,
    n_classes=2, num_features_per_signal_type=25, noise_gradient=10,
    initial_noise_sigma=10, n_samples_per_group=1, interpolatable=True, random_state=None,
    singularity=None, clip=None ):
    """
    Analogous to sklearn.datasets.make_classification, but simplified.

    Number of samples is reduced to be evenly divisible by number of classes.

    The number of features contained in the feature set will be a multiple of the number of
    signals contained in the dict "signals" (len(signals) == 12 at the time of writing, so
    len( new_fs.feature_names ) == 12, 24, 36, etc ).

    The more features the user asks for via the argument "num_features_per_signal_type",
    the more each successive feature generated using that signal will have a greater degree of
    gaussian noise added to it (controlled via args "noise_gradient" and "initial_noise_sigma").
    The features inside each individual feature set will come from a single signal function,
    and will get progressively noisier and noisier via the noise gradient, which is the
    multiplier for the sigma term in the gaussian noise generator.

    n_samples_per_group is related to tiling, i.e., how many ROIS/tiles within an image
    will have features calculated.

    interpolatiable sets the self.interpolation_coefficients, off of which interpolated
    value functionalitiy is keyed.
    """

    if random_state:
        from numpy.random import RandomState

        if random_state is True:
            from numpy.random import normal
        elif type( random_state ) is RandomState:
            normal = random_state.normal
        elif type( random_state ) is int:
            normal = RandomState( random_state ).normal
        else:
            raise ValueError( 'Arg random_state must be an instance of np.random.RandomState, an int, or the value True')
    else: # no noise added to feature values
        normal = lambda mu, sigma, n: np.zeros( n )

    if n_samples_per_group < 1 or type(n_samples_per_group) is not int:
        raise ValueError( "n_samples_per_group has to be an integer, and at least 1." )

    # Figure out what signals to use
    if singularity:
        # combine dicts of signals
        from itertools import chain
        signals = dict(chain(well_behaved_signals.iteritems(), singularity_signals.iteritems()))
    else:
        signals = well_behaved_signals

    # Use signal clipping?
    try:
        if len( clip ) == 2:
            hi, lo= float( clip[0] ), float( clip[1] )
            from functools import partial
            clip = partial( np.clip, a_min=lo, a_max=hi )
        else:
            raise ValueError( "Arg clip requires True, false, or a 2-tuple of numbers" )
    except TypeError:
        if clip:
            from functools import partial
            clip = partial( np.clip, a_min=lbound, a_max=ubound )
        else:
            clip = lambda x: x

    # num_samples must be evenly divisible by the number of samples per sample group
    # and the number of classes - therefore the final number of total samples may be
    # less than the user asked for.
    n_samples_per_class = int( n_samples // n_classes )
    n_samplegroups_per_class = int( n_samples_per_class // n_samples_per_group )
    n_samples_per_class = int( n_samplegroups_per_class * n_samples_per_group ) # change number inputted
    n_samples = int( n_samples_per_class * n_classes ) # changes number inputted!
    if n_samples <= 0:
        raise ValueError( "Specify n_samples to be a multiple of n_classes ({0}) * n_samples_per_group ({1}) >= {2}".format( n_classes, n_samples_per_group, n_classes* n_samples_per_group ) )
    # Initialize the basic data members

    num_features = num_features_per_signal_type * len( signals )

    # Instantiate and assign basic data members
    new_fs = FeatureSpace( name=name, source_filepath=name, num_samples=n_samples,
        num_samples_per_group=n_samples_per_group, num_features=num_features, discrete=True,
        feature_set_version='-1.0')

    new_fs.num_classes = n_classes
    new_fs.class_sizes = [ n_samples_per_class for i in xrange( n_classes ) ]
 
    # The function call np.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_classes)

    # The artificial class names are chosen to be numbers evenly spaced between 
    # lbound=-100 and ubound=100, which is the interval within which the functions
    # chosen to create the artificial features behave "interestingly."
    new_fs.interpolation_coefficients = list( np.mgrid[ lbound : ubound : step ] )

    new_fs.class_names = []
    for val in new_fs.interpolation_coefficients:
        # Try to use integers for each individual class name
        if float(int(val)) == val:
            new_fs.class_names.append( "FakeClass{0:02d}".format( int(val) ) )
        else:
            new_fs.class_names.append( "FakeClass{0:.02f}".format( val ) )

    # Generate artificial feature names
    # N.B. The feature generation signals are sorted in alphanum order!
    new_fs.feature_names = [ "{0}{1:04d}".format(fname, i)\
                                    for fname in sorted( signals.keys() ) \
                                    for i in xrange( num_features_per_signal_type ) ]
    if n_samples_per_group >= 1:
        group_index = 0

    # Creating sample metadata
    if n_samples_per_group == 1:
      # Format: <class name>_<sample index within class>
      new_fs._contiguous_sample_names = \
        [ "{0}_{1:03d}".format( class_name, i ) \
             for class_name in new_fs.class_names for i in xrange( n_samples_per_class ) ]
      new_fs._contiguous_sample_group_ids = range( n_samples ) # not xrange
      new_fs._contiguous_sample_sequence_ids =  [1] * n_samples
    else:
      # Format: <class name>_i<sample index within class>_g<group>_t<tile#>
      temp1 = [ "{0}_i{1:03d}".format( class_name, i ) \
        for class_name in new_fs.class_names for i in xrange( n_samples_per_class ) ]
      n_samplegroups = n_samples / n_samples_per_group
      temp2 = [ "_g{0:03d}_t{1:02d}".format( samplegroup_index, tile_index ) \
        for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_sample_names = [ a + b for a, b in zip( temp1, temp2 ) ]
      new_fs._contiguous_sample_group_ids = \
          [ samplegroup_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_sample_sequence_ids = \
          [ tile_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]

    old_settings = np.seterr( all='ignore' ) # Bring on the NaN's!

    # Create features across all classes at the same time
    feat_count = 0
    # N.B. The features are in sort order!
    ground_truth_values = np.array( new_fs.interpolation_coefficients )
    for func_name in sorted( signals.keys() ):
        f = signals[ func_name ]
        raw_class_feature_values = clip( f( ground_truth_values ) )
        raw_feature_values = np.empty( n_samples, )
        for i, val in enumerate( raw_class_feature_values ):
            raw_feature_values[ i * n_samples_per_class: (i+1) * n_samples_per_class].fill( val )

        for feat_index in xrange( num_features_per_signal_type ):
            # Add noise proportional to the feature index
            noise_vector = normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples )
            new_fs.data_matrix[:,feat_count] = np.add( noise_vector, raw_feature_values )
            feat_count += 1

    np.seterr( **old_settings )

    # discrete data always gets labels and sometimes gets values
    new_fs._contiguous_ground_truth_labels = [ new_fs.class_names[i] \
        for i in xrange( n_classes ) \
            for j in xrange( n_samples_per_class ) ]

    if not interpolatable:
        # delete the coefficients if user asks for a pure classification problem feat. set.
        new_fs.interpolation_coefficients = None
    else:
        new_fs._contiguous_ground_truth_values = [ new_fs.interpolation_coefficients[i] \
            for i in xrange( n_classes ) for j in range( n_samples_per_class ) ]

    new_fs.samples_sorted_by_ground_truth = True
    new_fs._RebuildViews( recalculate_class_metadata=False )
    return new_fs
