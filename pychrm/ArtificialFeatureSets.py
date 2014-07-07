from math import sin, tan, e, atan, log
import numpy as np

lbound = -100
ubound = 100

# signal types (picked to have good range of output for inputs [-100,100])
signals = {
'constant'        : lambda x: 100,                       # y:100
'positive_linear' : lambda x: x,                         # y:x
'negative_linear' : lambda x: -1*x,                      # y:-x
'step'            : lambda x: -100 if x < 0 else 100,    # y:-100, x<0; y:100 x>0
'sinusoidal'      : lambda x: 100 * sin( 0.1 * x ),      # y:100*sin(0.1x)
'tangent'         : lambda x: 100 * tan( 0.1 * x ),      # y:100*tan(0.1x)
'exponential'     : lambda x: e ** x,                    # y:e^(0.1x)
'quadratic'       : lambda x: 0.01 * x ** 2,             # y:0.01x^2
'cubic'           : lambda x: 0.0002 * x ** 3,           # y:0.0002x^3
'asymptotal'      : lambda x: 50 * atan( 0.1 * x + 10 ), # y:50*atan(0.1x+10)
'triangle'        : lambda x: abs( -1 * x ),             # y:abs(-x)
}

def logarithmic( x ):
  try:
    return 10 * log( x )
  except ValueError: # math domain error b/c (undef for x < 0)
    return np.nan

signals[ 'logarithmic' ] = logarithmic


def CreateArtificialFeatureSet_Continuous( name="ContinuousArtificialFS", n_samples=100, num_features_per_signal_type=10, noise_gradient=5,
        initial_noise_sigma = 10):
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

    from .FeatureSet import FeatureSet_Continuous

    new_fs = FeatureSet_Continuous()
    new_fs.name = name
    new_fs.source_path = new_fs.name
    new_fs.feature_vector_version = '2.0'
 
    new_fs.data_matrix = np.empty( ( n_samples, num_features_per_signal_type * len( signals ) ) )
    # The function call np.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_samples)
    new_fs.ground_truths = list( np.mgrid[ lbound : ubound : step ] )
    new_fs.imagenames_list = [ "ArtificialSample{0:03d}".format( i ) for i in xrange( n_samples ) ]
    new_fs.num_images = n_samples
    new_fs.num_features = num_features_per_signal_type * len( signals ) 
    new_fs.contiguous_imagenames_list = new_fs.imagenames_list

    for func_index, func_name in enumerate( signals.keys() ):
        f = signals[ func_name ]
        raw_feature_values = map( f, new_fs.ground_truths )

        for feat_index in xrange( num_features_per_signal_type ):
            new_fs.featurenames_list.append( "{0}{1:04d}".format(func_name, feat_index) )
            # Add noise proportional to the feature index
            noise_vector = np.random.normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples )
            feat_col_index = feat_index + func_index * num_features_per_signal_type 
            new_fs.data_matrix[:,feat_col_index] = np.add( noise_vector, raw_feature_values )

    return new_fs

def CreateArtificialFeatureSet_Discrete( name="DiscreteArtificialFS", n_samples=100, n_classes=2, num_features_per_signal_type=10, noise_gradient=5,
    initial_noise_sigma = 10):
    """
    Analogous to sklearn.datasets.make_classification, but simplified.

    Number of samples is reduced to be evenly divisible by number of classes.

    The number of features contained in the feature set will be a multiple of the number of
    signals contained in the dict "signals" (len(signals) == 12 at the time of writing, so
    len( new_fs.featurenames_list ) == 12, 24, 36, etc ).

    The more features the user asks for via the argument "num_features_per_signal_type",
    the more each successive feature generated using that signal will have a greater degree of
    gaussian noise added to it (controlled via args "noise_gradient" and "initial_noise_sigma").
    """

    from .FeatureSet import FeatureSet_Discrete

    new_fs = FeatureSet_Discrete()
    new_fs.name = name
    new_fs.source_path = new_fs.name
    new_fs.feature_vector_version = '2.0'

    # num_samples must be evenly divisible by the number of classes
    n_samples_per_class = int( n_samples // n_classes )
    new_fs.num_images = n_samples = n_samples_per_class * n_classes
    new_fs.num_features = num_features_per_signal_type * len( signals ) 
    new_fs.num_classes = n_classes
    new_fs.classsizes_list = [ n_samples_per_class for i in range( n_classes ) ]
 
    # The function call np.mgrid() requires for the value indicating the number of steps
    # to be imaginary, for some inexplicable reason.
    step = complex(0, n_classes)
    new_fs.interpolation_coefficients = list( np.mgrid[ lbound : ubound : step ] )
    new_fs.classnames_list = []
    for val in new_fs.interpolation_coefficients:
        # Try to use integers
        if float(int(val)) == val:
            new_fs.classnames_list.append( "ArtificialClass_{0:02d}".format( int(val) ) )
        else:
            new_fs.classnames_list.append( "ArtificialClass_{0:.02f}".format( val )  )

    new_fs.featurenames_list = [ "{0}{1:04d}".format(fname, i)\
                                    for fname in sorted( signals.keys() ) \
                                    for i in xrange( num_features_per_signal_type ) ]

    for class_index in xrange( n_classes ):
        new_fs.imagenames_list.append( \
            [ "{0}_{1:03d}:".format( new_fs.classnames_list[ class_index ], i ) for i in range( n_samples_per_class ) ] )
        new_fs.data_list.append( np.empty( ( n_samples_per_class, new_fs.num_features ) ) )
        for func_name in sorted( signals.keys() ):
            f = signals[ func_name ]
            raw_feature_values = [ f( new_fs.interpolation_coefficients[ class_index ] ) ] * n_samples_per_class 
            for feat_index in xrange( num_features_per_signal_type ):
                # Add noise proportional to the feature index
                noise_vector = np.random.normal( 0, initial_noise_sigma + feat_index * noise_gradient, n_samples_per_class )
                new_fs.data_list[ class_index ][:,feat_index] = \
                        np.add( noise_vector, raw_feature_values )

    new_fs.ContiguousDataMatrix()
    return new_fs
