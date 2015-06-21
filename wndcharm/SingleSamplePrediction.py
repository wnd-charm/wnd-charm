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
class _SingleSamplePrediction( object ):
    """Base class, do not instantiate."""
    def __init__( self ):
        self.name = None
        self.source_filepath = None
        self.ground_truth_value = None
        self.predicted_value = None
        # Labels really only apply to classifications, but we establish the member
        # on the base class and initialize it to something on the daughter class
        # to facilitate polymorphism/sorting in PerSampleStatistics.
        self.ground_truth_label = None
        self.predicted_label = None
        self.split_number = None
        self.sample_group_id = None
        self.num_samples_in_group = 1
        #: Indicates the position of the ROI within an image
        self.sample_sequence_id = None
        self.discrete = None

    #==============================================================
    def __repr__( self ):
        return str( self )

    #==============================================================
    def __str__( self ):
        """I like to use getattr in here just in case you're debugging the constructor"""

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
        if self.sample_group_id is not None:
            outstr += ' grp=' + str( self.sample_group_id )
        if self.sample_sequence_id is not None:
            outstr += ' seq=' + str( self.sample_sequence_id )

        # Classification specific attribs,
        # a.k.a if self.__class__ == SingleSampleClassification:
        if self.discrete == True:
            if self.predicted_label is not None:
                outstr += ' pred="' + self.predicted_label + '"'
            gt_label = getattr( self, 'ground_truth_label', None )
            if gt_label is not None:
                outstr += ' act="' + gt_label + '"'
            mps = getattr( self, 'marginal_probabilities', None )
            if mps is not None:
                outstr += ' marg probs={'
                for val in mps:
                    outstr += "{0:0.3f},".format( val )
                outstr += '}'

        if self.predicted_value is not None:
            outstr += " pred={0:0.2f}".format( self.predicted_value )

        # Averaging specific attribs,
        # a.k.a if self.__class__ == AveragedSingleSamplePrediction:
        std = getattr( self, 'std', None )
        if std is not None:
            outstr += " std={0:0.2f}".format( std )

        res_list = getattr( self, 'individual_results', None )
        if res_list is not None:
            outstr += " n={0}".format( len( res_list ) )

        return outstr + '>'

    #==============================================================
    @output_railroad_switch
    def Print( self, line_item=False, include_name=True, include_split_number=False,
            include_col_header=False, training_set_class_names=None ):
        """Output tab-separated prediction results data"""

        if not line_item:
            print str(self)
            return

        outstr = ""
        col_header = ""

        # sample name:
        if include_name and self.source_filepath:
            outstr += self.source_filepath + '\t'
            if include_col_header:
                col_header += "Samp. Name\t"

        if self.sample_sequence_id is not None:
            if self.sample_sequence_id == 'AVG':
                outstr += "(AVG)"
            elif self.num_samples_in_group != 1:
                outstr += "({0}/{1})".format( self.sample_sequence_id + 1, self.num_samples_in_group )
            outstr += '\t'
            if include_col_header:
                col_header += "ROI Index\t"

        # split number
        if include_split_number and self.split_number is not None:
            outstr += '{0}\t'.format( self.split_number )
            if include_col_header:
                col_header += "Split\t"

        # Classification stuff
        if self.discrete:
            # normalization factor:
            if self.normalization_factor is None:
                # no normalization factor/ marg probs means this is a non-call
                # FIXME: this is a very WND-CHARM specific thing
                # Maybe people can just put 'N/A' in there if they want to use this object?
                outstr += "--COLLISION--\t"
            else:
                outstr += "{0:0.3g}\t".format( self.normalization_factor )
            if include_col_header:
                col_header += "Norm. Fact.\t"

            # marginal probabilities:
            if self.marginal_probabilities is not None:
                outstr += "\t".join(\
                    [ "{0:0.3f}".format( prob ) for prob in self.marginal_probabilities ] )

                if include_col_header:
                    if training_set_class_names:
                        col_header += "\t".join( [ "p(" + class_name + ")" \
                            for class_name in training_set_class_names ] )
                    else:
                        col_header += "\t".join( [ "p(Class{0})".format( val ) \
                            for val in xrange( len( self.marginal_probabilities ) ) ] )
                    col_header += "\t"
            outstr += "\t"

            # actual class:
            if self.ground_truth_label is not None:
                outstr += self.ground_truth_label + '\t'
            else:
                outstr += "*\t"
            if include_col_header:
                col_header += "Act. Class\t"

            # predicted class:
            outstr += self.predicted_label + "\t"
            if include_col_header:
                col_header += "Pred. Class\t"

        # predicted value
        if self.predicted_value is not None:
            outstr += "{0:0.3f}\t".format( self.predicted_value )
            if include_col_header:
                col_header += "Pred. Val.\t"

        # Averaging specific attribs,
        # a.k.a if self.__class__ == AveragedSingleSamplePrediction:
        res_list = getattr( self, 'individual_results', None )
        if res_list is not None:
            outstr += "{0}\t".format( len( res_list ) )
            if include_col_header:
                col_header += "N calls\t"

            std = getattr( self, 'std', None )
            if std is not None:
                outstr += "{0:0.2f}\t".format( std )
                if include_col_header:
                    col_header += "Pred. Val. StDev\t"

        if include_col_header:
            print col_header

        print outstr

#=================================================================================
class SingleSampleClassification( _SingleSamplePrediction ):
    """Classification result for a single image/ROI (a.k.a "sample"),
    which includes predicted class, marginal probabilities, etc."""

    def __init__( self ):
        super( SingleSampleClassification, self ).__init__()

        self.marginal_probabilities = []
        self.normalization_factor = None
        self.marginal_probabilities = None

        #: predicted_label will always be a string
        #: the interpolated value, if applicable, gets stored in self.predicted_vlaue
        self.predicted_label = 'UNKNOWN'
        self.ground_truth_label = 'UNKNOWN'
        self.discrete = True

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

        if test_samp.sample_group_id is not None:
            result.sample_group_id = test_samp.sample_group_id
        if test_samp.sample_sequence_id is not None:
            result.sample_sequence_id = test_samp.sample_sequence_id

        marg_probs = np.array( result.marginal_probabilities )
        result.predicted_label = training_set.class_names[ marg_probs.argmax() ]
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
class SingleSampleRegression( _SingleSamplePrediction ):
    """Predicted value result of a regression of a single image/ROI (a.k.a. "sample")."""

    #==============================================================
    def __init__( self ):
        """Constructor"""
        super( SingleSampleRegression, self ).__init__()
        self.discrete = False

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

#=================================================================================
class AveragedSingleSamplePrediction( _SingleSamplePrediction ):
    """Classification result representing a single sample classified
    multiple times.

    Marginal probability lists have to be the same length."""

    def __init__( self, list_of_results, training_set_class_names=None ):
        super( AveragedSingleSamplePrediction, self ).__init__()

        #alias:
        reslist = list_of_results

        # This is some ghetto polymorphism right here!
        self.discrete = isinstance( reslist[0], SingleSampleClassification )

        # save 'em; useful when you want to print out the individual calls
        # to see how predictions change over the splits, ala PerSampleStatistics()
        self.individual_results = reslist

        self.source_filepath = reslist[0].source_filepath
        self.ground_truth_value = reslist[0].ground_truth_value
        self.ground_truth_label = reslist[0].ground_truth_label

        sample_sequence_ids = set( [ res.sample_sequence_id for res in reslist ] )
        if sample_sequence_ids != set( (None,) ) and len( sample_sequence_ids ) == 1:
            self.sample_sequence_id = sample_sequence_ids.pop()
        else:
            self.sample_sequence_id = 'AVG'

        if self.discrete:
            # Sometimes the result comes back with a non-call, like when
            # the sample image collides with every test image.
            # In that case no norm factor and no marg probs.
            norm_facts = [ res.normalization_factor for res in reslist \
                                                 if res.normalization_factor ]
            self.normalization_factor = np.array( norm_facts ).mean()

            # To take per-class averages of marginal probabilities,
            # collect all m.p.'s in a long list, reshape, and do a numpy per-column mean.
            mp_list_of_lists = [ res.marginal_probabilities for res in reslist \
                                                if res.marginal_probabilities  ]

            self.num_samples_in_group = n_res = len( mp_list_of_lists )
            mp_reject_count = len( reslist ) - len( mp_list_of_lists )
            if mp_reject_count == len( reslist ):
                self.predicted_label = "*Collided with every training image*"
            else:
                all_mps = np.array( [ mp for mp_list in mp_list_of_lists for mp in mp_list ] )
                all_mps.resize( ( n_res, len(all_mps) / n_res ) )
                avg_mps = all_mps.mean( axis=0 )
                self.marginal_probabilities = list( avg_mps )
                if training_set_class_names is not None:
                    self.predicted_label = training_set_class_names[ avg_mps.argmax() ]

        # There may be no predicted values
        pred_vals = [ res.predicted_value for res in reslist ]
        pred_vals_set = set( pred_vals )
        self.std = None
        if pred_vals_set != set( (None,) ):
            if len( pred_vals_set ) == 1:
                self.predicted_value = pred_vals_set.pop()
                self.std = 0
            else:
                pred_vals = np.array( pred_vals )
                self.predicted_value = pred_vals.mean()
                self.std = pred_vals.std()
