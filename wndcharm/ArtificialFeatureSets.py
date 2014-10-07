from math import sin, tan, e, atan, log
import numpy

lbound = -100
ubound = 100

# signal types (picked to have good range of output for inputs [-100,100])
signals = {
'constant'        : lambda x: 100,                       # y:100
'positive_linear' : lambda x: x,                         # y:x
'negative_linear' : lambda x: -1*x,                      # y:-x
'step'            : lambda x: -100 if x < 0 else 100,    # y:-100, x<0; y:100 x>0
'sinusoidal'      : lambda x: 100 * sin( 0.1 * x ),      # y:100*sin(0.1x)
'tangent'         : lambda x: tan( 0.0156 * x ),         # y:tan(0.0156x) # corners are distinct
'exponential'     : lambda x: e ** (0.05 * x ),          # y:e^(0.7x) => e^7 = 1096
'quadratic'       : lambda x: 0.01 * x ** 2,             # y:0.01x^2
'cubic'           : lambda x: 0.0002 * x ** 3,           # y:0.0002x^3
'asymptotal'      : lambda x: 50 * atan( 0.1 * x + 10 ), # y:50*atan(0.1x+10)
'triangle'        : lambda x: abs( -1 * x ),             # y:abs(-x)
}

def logarithmic( x ):
  try:
    return 10 * log( x )
  except ValueError: # math domain error b/c (undef for x < 0)
    return numpy.nan

#signals[ 'logarithmic' ] = logarithmic


def CreateArtificialFeatureSet_Continuous( name="ContinuousArtificialFS", n_samples=100,
    num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
    n_samples_per_group=1, random_state=None ):
    """
    Analogous to sklearn.datasets.make_regression, but simplified.

    Returns an instance of a FeatureSet_Continuous.

    The number of features contained in the feature set will be a multiple of the number of
    signals contained in the dict "signals" (len(signals) == 12 at the time of writing, so
    len( new_fs.featurenames_list ) == 12, 24, 36, etc ).

    The more features the user asks for via the argument "num_features_per_signal_type",
    the more each successive feature generated using that signal will have a greater degree of
    gaussian noise added to it (controlled via args "noise_gradient" and "initial_noise_sigma").
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
        raise ValueError( 'Arg random_state must be an instance of numpy.random.RandomState, an int, or the value True')
    else: # no noise added to feature values
      normal = lambda mu, sigma, n: numpy.zeros( n )

    from .FeatureSet import FeatureSet_Continuous

    new_fs = FeatureSet_Continuous()
    new_fs.name = name
    new_fs.source_path = new_fs.name
    new_fs.feature_vector_version = '2.0'
    new_fs.discrete = False

    if n_samples_per_group > 1:
      # make n_samples evenly divisible
      n_samples = int( n_samples // n_samples_per_group ) * n_samples_per_group
      new_fs.num_samples_per_group = n_samples_per_group
    new_fs.num_images = n_samples
    new_fs.num_features = num_features_per_signal_type * len( signals )
    new_fs.data_matrix = numpy.empty( ( n_samples, new_fs.num_features ) )
    # The function call numpy.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_samples)
    new_fs._contiguous_ground_truths = list( numpy.mgrid[ lbound : ubound : step ] )

    # Generate artificial feature names
    # N.B. The feature generation signals are sorted in alphanum order!
    new_fs.featurenames_list = [ "{0}{1:04d}".format(fname, i)\
                                    for fname in sorted( signals.keys() ) \
                                    for i in xrange( num_features_per_signal_type ) ]
    # Creating sample metadata
    if n_samples_per_group == 1:
      new_fs._contiguous_samplenames_list = [ "FakeContinuousSample{0:03d}".format( i )\
                                             for i in xrange( n_samples ) ]
      new_fs._contiguous_samplegroupid_list = range( n_samples ) # not xrange
      new_fs._contiguous_samplesequenceid_list =  [1] * n_samples
    else:
      # Format: FakeContinuousSample_i<sample index>_g<group>_t<tile#>
      temp1 = [ "FakeContinuousSample_i{0:03d}".format( i ) for i in xrange( n_samples ) ]
      n_samplegroups = n_samples / n_samples_per_group
      temp2 = [ "_g{0:03d}_t{1:02d}".format( samplegroup_index, tile_index ) \
        for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_samplenames_list = [ a + b for a, b in zip( temp1, temp2 ) ]
      new_fs._contiguous_samplegroupid_list = \
          [ samplegroup_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_samplesequenceid_list = \
          [ tile_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]

    # Create features across all classes at the same time
    feat_count = 0
    # N.B. The features are in sort order!
    for func_name in sorted( signals.keys() ):
      f = signals[ func_name ]
      raw_feature_values = map( f, new_fs._contiguous_ground_truths )
      for feat_index in xrange( num_features_per_signal_type ):
        # Add noise proportional to the feature index
        noise_vector = normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples )
        new_fs.data_matrix[:,feat_count] = numpy.add( noise_vector, raw_feature_values )
        feat_count += 1

    new_fs._RebuildViews()
    return new_fs

def CreateArtificialFeatureSet_Discrete( name="DiscreteArtificialFS", n_samples=100,
    n_classes=2, num_features_per_signal_type=10, noise_gradient=5,
    initial_noise_sigma=10, n_samples_per_group=1, interpolatable=True, random_state=None ):
    """
    Analogous to sklearn.datasets.make_classification, but simplified.

    Number of samples is reduced to be evenly divisible by number of classes.

    The number of features contained in the feature set will be a multiple of the number of
    signals contained in the dict "signals" (len(signals) == 12 at the time of writing, so
    len( new_fs.featurenames_list ) == 12, 24, 36, etc ).

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
        raise ValueError( 'Arg random_state must be an instance of numpy.random.RandomState, an int, or the value True')
    else: # no noise added to feature values
      normal = lambda mu, sigma, n: numpy.zeros( n )

    from .FeatureSet import FeatureSet_Discrete

    if n_samples_per_group < 1 or type(n_samples_per_group) is not int:
      raise ValueError( "n_samples_per_group has to be an integer, and at least 1." )

    # Initialize the basic data members
    new_fs = FeatureSet_Discrete()
    new_fs.name = name
    new_fs.source_path = new_fs.name
    new_fs.feature_vector_version = '2.0'
    new_fs.discrete = True
    new_fs.num_samples_per_group = n_samples_per_group

    # num_samples must be evenly divisible by the number of samples per sample group
    # and the number of classes - therefore the final number of total samples may be
    # less than the user asked for.
    n_samples_per_class = int( n_samples // n_classes )
    n_samplegroups_per_class = int( n_samples_per_class // n_samples_per_group )
    n_samples_per_class = int( n_samplegroups_per_class * n_samples_per_group ) # change number inputted
    n_samples = int( n_samples_per_class * n_classes ) # changes number inputted!
    if n_samples <= 0:
      raise ValueError( "Specify n_samples to be a multiple of n_classes ({0}) * n_samples_per_group ({1}) >= {2}".format( n_classes, n_samples_per_group, n_classes* n_samples_per_group ) )
    new_fs.num_images = n_samples
    new_fs.num_features = num_features_per_signal_type * len( signals ) 
    new_fs.data_matrix = numpy.empty( ( n_samples, new_fs.num_features ) )
    new_fs.num_classes = n_classes
    new_fs.classsizes_list = [ n_samples_per_class for i in xrange( n_classes ) ]
 
    # The function call numpy.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_classes)

    # The artificial class names are chosen to be numbers evenly spaced between 
    # lbound=-100 and ubound=100, which is the interval within which the functions
    # chosen to create the artificial features behave "interestingly."
    new_fs.interpolation_coefficients = list( numpy.mgrid[ lbound : ubound : step ] )

    new_fs.classnames_list = []
    for val in new_fs.interpolation_coefficients:
      # Try to use integers for each individual class name
      if float(int(val)) == val:
        new_fs.classnames_list.append( "FakeClass{0:02d}".format( int(val) ) )
      else:
        new_fs.classnames_list.append( "FakeClass{0:.02f}".format( val )  )

    # Generate artificial feature names
    # N.B. The feature generation signals are sorted in alphanum order!
    new_fs.featurenames_list = [ "{0}{1:04d}".format(fname, i)\
                                    for fname in sorted( signals.keys() ) \
                                    for i in xrange( num_features_per_signal_type ) ]

    if n_samples_per_group >= 1:
      group_index = 0

    # Creating sample metadata
    if n_samples_per_group == 1:
      # Format: <class name>_<sample index within class>
      new_fs._contiguous_samplenames_list = \
        [ "{0}_{1:03d}".format( class_name, i ) \
             for class_name in new_fs.classnames_list for i in xrange( n_samples_per_class ) ]
      new_fs._contiguous_samplegroupid_list = range( n_samples ) # not xrange
      new_fs._contiguous_samplesequenceid_list =  [1] * n_samples
    else:
      # Format: <class name>_i<sample index within class>_g<group>_t<tile#>
      temp1 = [ "{0}_i{1:03d}".format( class_name, i ) \
        for class_name in new_fs.classnames_list for i in xrange( n_samples_per_class ) ]
      n_samplegroups = n_samples / n_samples_per_group
      temp2 = [ "_g{0:03d}_t{1:02d}".format( samplegroup_index, tile_index ) \
        for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_samplenames_list = [ a + b for a, b in zip( temp1, temp2 ) ]
      new_fs._contiguous_samplegroupid_list = \
          [ samplegroup_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]
      new_fs._contiguous_samplesequenceid_list = \
          [ tile_index for samplegroup_index in xrange( n_samplegroups ) \
          for tile_index in xrange( n_samples_per_group ) ]

    # Create features across all classes at the same time
    feat_count = 0
    # N.B. The features are in sort order!
    for func_name in sorted( signals.keys() ):
      f = signals[ func_name ]
      raw_class_feature_values = map( f, new_fs.interpolation_coefficients )
      raw_feature_values = numpy.empty( n_samples, )
      for i, val in enumerate( raw_class_feature_values ):
        raw_feature_values[ i * n_samples_per_class: (i+1) * n_samples_per_class].fill( val )

      for feat_index in xrange( num_features_per_signal_type ):
        # Add noise proportional to the feature index
        noise_vector = normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples )
        new_fs.data_matrix[:,feat_count] = numpy.add( noise_vector, raw_feature_values )
        feat_count += 1

    if not interpolatable:
      # delete the coefficients if user asks for a pure classification problem feat. set.
      new_fs.interpolation_coefficients = None
      # Use class labels instead of class values
      new_fs._contiguous_ground_truths = [ new_fs.classnames_list[i] \
          for i in xrange( n_classes ) \
            for j in xrange( n_samples_per_class ) ]
    else:
      new_fs._contiguous_ground_truths = [ new_fs.interpolation_coefficients[i] \
        for i in xrange( n_classes ) for j in range( n_samples_per_class ) ]

    new_fs._RebuildViews()
    return new_fs
