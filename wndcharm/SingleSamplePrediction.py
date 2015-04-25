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
import wndcharm # for ImageMatrix
from .utils import output_railroad_switch
from .FeatureVector import FeatureVector
from .FeatureWeights import FeatureWeights
from .FeatureSpace import FeatureSpace

#=================================================================================
class SingleSamplePrediction( object ):
    """Base class to contain prediction results for a single image/ROI (a.k.a "sample"),
    which includes predicted class, marginal probabilities, etc."""

    def __init__( self ):
        """Constructor"""

        self.name = None
        self.source_filepath = None
        self.ground_truth_value = None
        self.predicted_value = None
        self.batch_number = None
        self.samplegroupid = None
        self.samplesequenceid = None

        #: Indicates the position of the ROI within an image
        self.tile_index = None

    #==============================================================
    def __repr__( self ):
        return str( self )

#=================================================================================
class SingleSampleClassification( SingleSamplePrediction ):
    """Classification result for a single image/ROI (a.k.a "sample"),
    which includes predicted class, marginal probabilities, etc."""

    def __init__( self ):
        """Constructor"""
        super( SingleSampleClassification, self ).__init__()

        self.marginal_probabilities = []
        self.normalization_factor = None
        self.marginal_probabilities = None

        #: predicted_class_name will always be a string
        #: the interpolated value, if applicable, gets stored in self.predicted_vlaue
        self.predicted_class_name = 'UNKNOWN'
        self.ground_truth_class_name = 'UNKNOWN'

    #==============================================================
    def __eq__( self, other ):

        for a, b in zip( self.marginal_probabilities, other.marginal_probabilities ):
            # FIXME: For now, define equal as marg_probs agreeing within 5%
            if abs( a-b ) > 0.005:
                return False
        # The following is pretty useless since normalization factors are often in the 10e-20
        # range.
        #if abs( self.normalization_factor - other.normalization_factor ) > 0.001:
        #    return False
        return True

    #==============================================================
    def __ne__( self, other ):
        return not self == other

    #==============================================================
    @output_railroad_switch
    def Print( self, line_item=False ):
        """Output classification line item data, including predicted class and marginal
        probabilities"""
        
        if line_item:
            # img name:
            outstr = self.source_filepath if self.source_filepath else ""
            if self.tile_index is not None:
                if self.tile_index == 'AVG':
                    outstr += " (AVG)"
                else:
                    outstr += " ({0}/{1})".format( self.tile_index + 1, self.num_samples_in_group )
            # normalization factor:
            if self.normalization_factor is None:
                # no normalization factor means this is a non-call
                print outstr + "\t--COLLISION--"
                return
            outstr += "\t{0:0.3g}\t".format( self.normalization_factor )

            # marginal probabilities:
            outstr += "\t".join(\
                     [ "{0:0.3f}".format( prob ) for prob in self.marginal_probabilities ] )
            outstr += "\t"
            # actual class:
            if self.ground_truth_class_name:
                outstr += "{0}\t".format( self.ground_truth_class_name )
            else:
                outstr += "*\t"
            # predicted class:
            outstr += self.predicted_class_name + "\t"
            # interpolated value, if applicable
            if self.predicted_value is not None:
                outstr += "{0:0.3f}".format( self.predicted_value )
            print outstr
        else:
            print str(self)

    #==============================================================
    def __str__( self ):
    
        outstr = '<' + self.__class__.__name__

        if self.name:
            samp_name = self.name
        elif self.source_filepath:
            samp_name = self.source_filepath
        else:
            samp_name = ""

        if len( samp_name ) > 25:
            samp_name = '...' + samp_name[ -25: ]

        outstr += ' "' + samp_name + '"'
        if self.samplegroupid is not None:
            outstr += ' grp=' + str( self.samplegroupid )
        if self.samplesequenceid is not None:
            outstr += ' seq=' + str( self.samplesequenceid )
        if self.predicted_class_name:
            outstr += ' pred="' + self.predicted_class_name + '"'
        if self.ground_truth_class_name:
            outstr += ' act="' + self.ground_truth_class_name + '"'
        if self.marginal_probabilities is not None:
            outstr += ' marg probs={'
            for val in self.marginal_probabilities:
                outstr += "{0:0.3f},".format( val )
            outstr += '}'
        if self.predicted_value is not None:
            outstr += " interp={0:0.2f}".format( self.predicted_value )
        return outstr + '>'

    #==============================================================
    @classmethod
    def _WND5( cls, trainingset, testimg, feature_weights ):
        """
        Don't call this function directly, use the wrapper functions 
        FeatureSpaceClassification.New() (for FeatureSets) or
        SingleSampleClassification.NewWND5() (for single images/ROIs).
        Both of these functions have dummyproofing.

        For N images and M features:
            trainingset is list of length L of N x M numpy matrices
            testtile is a 1 x M list of feature values
        NOTE: the trainingset and test image must have the same number of features!!!
        AND: the features must be in the same order!!
        Returns an instance of the class SingleSampleClassification
        """

        #print "classifying..."
        epsilon = np.finfo( np.float ).eps

        num_features_in_testimg = len( testimg ) 
        weights_squared = np.square( feature_weights )

        # initialize
        class_similarities = [0] * trainingset.num_classes

        for class_index in range( trainingset.num_classes ):
            #print "Calculating distances to class {0}".format( class_index )
            num_tiles, num_features = trainingset.data_list[ class_index ].shape
            assert num_features_in_testimg == num_features,\
            "num features {0}, num features in test img {1}".format( num_features, num_test_img_features )

            # create a view
            sig_matrix = trainingset.data_list[ class_index ]
            wnd_sum = 0
            num_collisions = 0

            #print "num tiles: {0}, num_test_img_features {1}".format( num_tiles, num_test_img_features )
            for tile_index in range( num_tiles ):
                # epsilon checking for each feature is too expensive
                # FIXME: Do quick and dirty summation check until we can figure something else out
                dists = np.absolute( sig_matrix[ tile_index ] - testimg )
                w_dist = np.sum( dists )
#                print "train img {0} dist : {1}".format( tile_index, w_dist )
#                if (np.isinf(w_dist)):
#                    print "dists: "+str(dists)
                if w_dist < epsilon:
                    num_collisions += 1
                    continue
                dists = np.multiply( weights_squared, np.square( dists ) )
                w_dist = np.sum( dists )
                # The exponent -5 is the "5" in "WND5"
                class_similarities[ class_index ] += w_dist ** -5
            #print "\n"

            denom = num_tiles - num_collisions
            if denom == 0:
                # This sample collided with every sample in the test set
                # return a non-call
                return cls()
            class_similarities[ class_index ] /= denom
#            print "class_similarities: "+str(class_similarities)

        result = cls()
        norm_factor = sum( class_similarities )
        result.normalization_factor = norm_factor 
        result.marginal_probabilities = [ x / norm_factor for x in class_similarities ]
        return result

    #=================================================================================
    @classmethod
    def NewWND5( cls, training_set, feature_weights, test_samp, quiet = False ):
        """@brief: A wrapper function for _ClassifyOneImageWND5 that does dummyproofing
        @return: An instance of a FeatureSpaceClassification"""

        if not isinstance( training_set, FeatureSpace ):
            raise ValueError( 'First argument to NewWND5 must be of type "FeatureSpace", you gave a {0}'.format( type( training_set ).__name__ ) )
        
        if not isinstance( feature_weights, FeatureWeights ):
            raise ValueError( 'Second argument to NewWND5 must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )

        if not isinstance( test_samp, FeatureVector ):
            raise ValueError( 'Third argument to NewWND5 must be of type "FeatureVector", you gave a {0}'.format( type( test_samp ).__name__ ) )

        train_set_len = len( training_set.feature_names )
        test_set_len = len( test_samp.feature_names )
        feature_weights_len = len( feature_weights.feature_names )

        if test_samp.feature_names != feature_weights.feature_names:
            raise ValueError("Can't classify, features in signature don't match features in weights." )

        if test_samp.feature_names != training_set.feature_names:
            raise ValueError("Can't classify, features in signature don't match features in training_set." )

        if not quiet:
            print "Classifying image '{0}' ({1} features) against test set '{2}' ({3} features)".\
             format( test_samp.name, train_set_len, training_set.name, test_set_len )

        result = cls._WND5( training_set, test_samp.values, feature_weights.values )

        if isinstance( test_samp.source_filepath, wndcharm.ImageMatrix ) and \
                test_samp.source_filepath.source:
            result.source_filepath = test_samp.source_filepath.source
        else:
            result.source_filepath = test_samp.name

        if test_samp.samplegroupid is not None:
            result.samplegroupid = test_samp.samplegroupid
        if test_samp.samplesequenceid is not None:
            resultsamplesequenceid = test_samp.samplesequenceid

        marg_probs = np.array( result.marginal_probabilities )
        result.predicted_class_name = training_set.class_names[ marg_probs.argmax() ]
        # interpolated value, if applicable

        if training_set.interpolation_coefficients is not None and \
                len( training_set.interpolation_coefficients ) == len( marg_probs ):
            interp_val = np.sum( marg_probs * training_set.interpolation_coefficients )
            result.predicted_value = interp_val

        if not quiet:
            column_header = "image\tnorm. fact.\t"
            column_header +=\
             "".join( [ "p(" + class_name + ")\t" for class_name in training_set.class_names ] )
            column_header += "act. class\tpred. class\tpred. val."
            print column_header
            result.Print( line_item=True )
        return result

#=================================================================================
class SingleSampleRegression( SingleSamplePrediction ):
    """Predicted value result of a regression of a single image/ROI (a.k.a. "sample")."""

    #==============================================================
    def __init__( self ):
        """Constructor"""
        super( SingleSampleRegression, self ).__init__()

    #==============================================================
    @output_railroad_switch
    def Print( self, line_item = False ):
        """Output results."""

        if line_item:
            # img name:
            output_str = str( self.source_filepath )
            output_str += "\t"
            # actual class:
            if self.ground_truth_value is not None:
                output_str += str( self.ground_truth_value ) + "\t"
            else:
                output_str += "*\t"
            # predicted class:
            output_str += str( self.predicted_value )
            print output_str
        else:
            str( self )

    #==============================================================
    def __str__( self ):
        outstr = '<' + self.__class__.__name__

        if self.name:
            samp_name = self.name
        elif self.source_filepath:
            samp_name = self.source_filepath
        else:
            samp_name = ""

        if len( samp_name ) > 25:
            samp_name = '...' + samp_name[ -25: ]

        outstr += ' "' + samp_name + '"'
        if self.predicted_value is not None:
            outstr += " interp={0.2f}".format( self.predicted_value )
        if self.ground_truth_class_name:
            outstr += ' act="' + self.ground_truth_val + '"'
        return outstr + '>'

    #==============================================================
    @classmethod
    def _MultivariateLinear( cls, one_image_features, feature_weights ):
        """Produce a predicted value for a single image based on the regression parameters
        contained in the PearsonFeatureWeights argument "feature_weights".
        
        Don't call this function directly, but instead use the member function
        New() on the FeatureSpaceRegression class which has dummyproofing."""

        per_feature_predicted_vals = []

        for i in range( len( one_image_features ) ):
            feat_val = one_image_features[i]
            weight = feature_weights.values[i]
            slope = feature_weights.slopes[i]
            intercept = feature_weights.intercepts[i]

            # y = mx+b
            # feature value = slope * age + intercept
            # solve for the abscissa:
            per_feature_predicted_vals.append( weight * ( feat_val - intercept ) / slope )

        result = cls()
        result.predicted_value = sum( per_feature_predicted_vals )

        return result

    #==============================================================
    @classmethod
    def _LeastSquares( cls, one_image_features, training_set ):
        """Produce a predicted value for a single image based on numpy.linalg.lstsq().
        
        Don't call this function directly, but instead use the member function
        NewLeastSquares() on the FeatureSpaceRegression class
        which has dummyproofing and performs the necessary matrix augmentation."""

        from numpy.linalg import lstsq
        from numpy import dot

        A = lstsq( training_set.data_matrix, np.array( training_set.ground_truth_values ) )[0]
        
        result = cls()
        result.predicted_value = dot( one_image_features, A )
        return result
