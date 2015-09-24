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
import wndcharm
from .utils import output_railroad_switch

#############################################################################
# class definition of FeatureWeights
#############################################################################
class FeatureWeights( object ):
    """FeatureWeights is an abstract base class
    that comprises one-half of a WND-CHARM classifier (the other being a FeatureSpace).

    It is a list of strings which are the names of
    individual image descriptors (features) and a corresponding list of doubles which
    are the weights assigned to those features. Since WND-CHARM is a generalized
    pattern recognition algorithm that calculates the same image descriptors for all
    images, it is through features weights that trained classifiers can zero-in on only
    those features which provide distinctiveness across classes and ignore noisy features.
    Thus any instance of a FeatureWeights class is context-specific."""

    def __init__( self, name=None, size=None ):
        self.name = name
        self.associated_feature_space = None
        if size is not None:
            self.feature_names = [None] * size
            self.values = [None] * size
        else:
            self.feature_names = None
            self.values = None
    #================================================================
    def __len__( self ):
        try:
            return len( self.feature_names )
        except:
            return 0
    #================================================================
    def __str__( self ):
        s = '<' + self.__class__.__name__
        if self.name:
            s += '"' + str( self.name ) + '"'
        s += " n_features=" + str( len(self.feature_names) )
        if len( self.feature_names ) > 0:
            s+= ' feat0="{0}" val0={1}'.format( self.feature_names[0], self.values[0] )
        s += '>'
        return s
    #================================================================
    def __repr__( self ):
        return str( self )
    #================================================================
    @classmethod
    def NewFromFile( cls, weights_filepath ):
        """Load feature weights from a C++ WND-CHARM-sytle feature weights text file,
        i.e., the one created by using the command line "wndchrm train -vw/path/to/weights.txt"""
        raise NotImplementedError

    #================================================================
    def Threshold( self, num_features_to_be_used  ):
        """@breif Returns an instance of a FeatureWeights class with the top n relevant features in that order"""
        raise NotImplementedError

    #================================================================
    @classmethod
    def NewFromFeatureSpace( cls, num_features_to_be_used  ):
        """@breif Calculate FeatureWeights from a FeatureSpace"""
        raise NotImplementedError

    #================================================================
    @output_railroad_switch
    def Print( self, display=None ):
        """@breif Prints out feature values and statistics"""
        raise NotImplementedError


#############################################################################
# class definition of FisherFeatureWeights
#############################################################################
class FisherFeatureWeights( FeatureWeights ):
    """
    FisherFeatureWeights is a concrete class using feature weighting that uses the
    Fisher Discriminant as a weight for image descriptors. Features that have a high
    Fisher score have high inter-class variability and low intra-class variability.

    This type of feature weight is one result of training an image classifier from a
    training set that has distinct discrete classes, e.g., is a cell mono- or binucleate.
    Fisher classifiers can be used on a continuous morphology regime if the user creates
    groups or bins for the training images.
    """

    def __init__( self, name=None, size=None ):
        """Simply calls parent constructor"""
        super( FisherFeatureWeights, self ).__init__( name=name, size=size )

    #================================================================
    @classmethod
    def NewFromFile( cls, weights_filepath ):
        """Written to read in files created by wndchrm -vw/path/to/weightsfile.txt"""

        weights = cls( name=weights_filepath )
        with open( weights_filepath ) as weights_file:
            # split line "number <space> name"
            raw_vals, raw_names = \
              zip( *[ line.strip().split( None, 1 ) for line in weights_file.read().splitlines() ] )

        weights.values = [ float( val ) for val in raw_vals ]
        weights.feature_names = [None] * len( raw_names )

        for i, name_raw_str in enumerate( raw_names ):
            # getFeatureInfoByName does some checking, returns a None if it can't parse it
            retval = wndcharm.FeatureNames.getFeatureInfoByName( name_raw_str )
            if retval:
                weights.feature_names[i] = retval.name
            else:
                weights.feature_names[i] = name_raw_str

        return weights

    #================================================================
    @classmethod
    def NewFromFeatureSpace( cls, fs ):
        """Takes a FeatureSpace as input and calculates a Fisher score for
        each feature. Returns a newly instantiated instance of FisherFeatureWeights.

        For:
        N = number of classes
        F = number of features
        It = total number of images in training set
        Ic = number of images in a given class
        """

        if fs.normalized_against is None:
            raise ValueError( "Before generating feature weights, call Normalize() of the feature space." )

        # we deal with NANs/INFs separately, so turn off numpy warnings about invalid floats.
        oldsettings = np.seterr(all='ignore')

        # 1D matrix 1 * F
        population_means = np.mean( fs.data_matrix, axis = 0 )

        # WARNING, this only works in python27:
        # ====================================
        # If 'fs' is a balanced training set (i.e., same number of images
        # in each class), you can use pure matrix calls without iteration:

        # 3D matrix N * Ic * F
        #all_images_classes_separate = np.array( fs.data_list )

        #if len( all_images_classes_separate.shape ) == 3:

            # 2D matrix N * F
        #    intra_class_means = np.mean( all_images_classes_separate, axis = 1 )
            # 2D matrix N * F
        #    intra_class_variances = np.var( all_images_classes_separate, axis = 1 )

        #else:
        # ====================================

        # 2D matrix shape N * F
        intra_class_means = np.empty( (fs.num_classes, len(fs.feature_names)) )
        # 2D matrix shape N * F
        intra_class_variances = np.empty( (fs.num_classes, len(fs.feature_names)) )

        class_index = 0
        for class_feature_matrix in fs.data_list:
            intra_class_means[ class_index ] = np.mean( class_feature_matrix, axis=0 )
            # Note that by default, numpy divides by N instead of the more common N-1, hence ddof=1.
            intra_class_variances[ class_index ] = np.var( class_feature_matrix, axis=0, ddof=1 )
            class_index += 1

        # 1D matrix 1 * F
        # we deal with NANs/INFs separately, so turn off numpy warnings about invalid floats.
        # for the record, in numpy:
        #     1./0. = inf, 0./inf = 0., 1./inf = 0. inf/0. = inf, inf/inf = nan
        #     0./0. = nan,  nan/0. = nan, 0/nan = nan, nan/nan = nan, nan/inf = nan, inf/nan = nan
        # We can't deal with NANs only, must also deal with pos/neg infs
        # The masked array allows for dealing with "invalid" floats, which includes nan and +/-inf
        denom = np.mean( intra_class_variances, axis=0 )
        denom[denom == 0] = np.nan
        feature_weights_m =  np.ma.masked_invalid (
            ( np.square( population_means - intra_class_means ).sum( axis = 0 ) /
            ( fs.num_classes - 1 ) ) / denom )
        # return numpy error settings to original
        np.seterr(**oldsettings)

        new_fw = cls()
        new_fw.feature_names = fs.feature_names[:]
        # the filled(0) method of the masked array sets all nan and infs to 0
        new_fw.values = feature_weights_m.filled(0).tolist()
        new_fw.associated_feature_space = fs

        return new_fw

    #================================================================
    def EliminateZeros( self ):
        """Eliminates any features with a weight of zero, and returns a new instance of
        FisherFeatureWeights without those features."""

        new_weights = FisherFeatureWeights()
        scores = zip( self.feature_names, self.values )
        nonzero_scores = [ (name, weight) for name, weight in scores if weight != 0 ]
        new_weights.feature_names, new_weights.values = zip( *nonzero_scores )
        return new_weights

    #================================================================
    def Threshold( self, num_features_to_be_used=None, _all=False ):
        """Returns an instance of a FisherFeatureWeights class with the top n relevant features
        in order.

        num_features_to_be_used can be integer n where 0 < n <= len( self.feature_names )
        or num_features_to_be_used can be float n where 0 < n <= 1

        if _all == True: simple reorder by feature rank returning even 0-weighted features
        If _all == 'nonzero': returns non-zero weighted features ranked by weight."""

        if _all:
            num_features_to_be_used = len( self.values )
        # Default is top 15% of features
        elif num_features_to_be_used is None:
            num_features_to_be_used = int( len( self.values ) * 0.15 )
        elif type( num_features_to_be_used ) is float:
            if num_features_to_be_used <= 0 or num_features_to_be_used > 1:
                raise ValueError('Choose feature reduction fraction on interval (0,1] (got "{0}"'.format( num_features_to_be_used  ) )
            num_features_to_be_used = int( round( num_features_to_be_used * len( self ) ) )
        elif num_features_to_be_used > len( self.values ) or num_features_to_be_used <= 0:
            raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
                                  format( len( self.values ), num_features_to_be_used ) )

        new_weights = self.__class__()
        raw_featureweights = zip( self.values, self.feature_names )
        # raw_featureweights is now a list of tuples := [ (value1, name1), (value2, name2), ... ]

        sorted_featureweights = sorted( raw_featureweights, key=lambda a: a[0], reverse = True )
        # take top N features
        from itertools import islice
        use_these_feature_weights = list( islice( sorted_featureweights, num_features_to_be_used ) )

        if _all is not True:
            # You have a problem if any of the features have corellation coefficients of 0
            for i, feat_fig_of_merit in enumerate( [ line[0] for line in use_these_feature_weights ] ):
                if feat_fig_of_merit == 0:
                    if _all == 'non-zero':
                        print "Features rank index {0} and below have a correllation coefficient of 0. ".format( i )
                        print 'Using {0} features'.format( i )
                        use_these_feature_weights = use_these_feature_weights[ : i ]
                        break
                    else:
                        err_msg = "Can't reduce feature weights \"{0}\" to {1} features. ".format( self.name, num_features_to_be_used )
                        err_msg += "Features ranked {0} and below have a Fisher score of 0. ".format( i )
                        err_msg += "Request less features. "
                        raise ValueError( err_msg )

        # we want lists, not tuples!
        new_weights.values, new_weights.feature_names =\
          [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

        new_weights.associated_feature_space = self.associated_feature_space

        return new_weights

    #================================================================
    def Slice( self, start_index=0, stop_index ):
        """Return a new instance of FisherFeatureWeights containing a chunk
        of middle-ranked features."""

        min_index = None
        max_index = None

        if stop_index > start_index:
            min_index = start_index
            max_index = stop_index
        else:
            min_index = stop_index
            max_index = start_index

        if (min_index < 0) or ( max_index > len( self.values ) ):
            raise ValueError( 'Cannot slice, check your start and stop indices.' )

        new_weights = self.__class__()
        raw_featureweights = zip( self.feature_names, self.values )

        from itertools import islice
        use_these_feature_weights = \
                list( islice( raw_featureweights, min_index, max_index ) )

        # we want lists, not tuples!
        new_weights.feature_names, new_weights.values =\
          [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

        new_weights.associated_feature_space = self.associated_feature_space

        return new_weights

    #================================================================
    @output_railroad_switch
    def Print( self, display=None ):
        """Prints out feature weight values and statistics.
        display (int) - number of features you want printed, from beginning"""

        if display != None:
            from itertools import islice
            features = islice( zip( self.values, self.feature_names ), display )
            remainder = len( self ) - display
        else:
            features = zip( self.values, self.feature_names )
            remainder = None

        s = self.__class__.__name__
        if self.name:
            s += ' "{0}:"'.format( self.name )
        s += " ({0} features)".format( len( self ) )
        print s
        print "Rank\tValue\tName"
        print "====\t=====\t===="
        for i, (val, name) in enumerate( features, start=1 ):
            print "{0}\t{1:.6f}\t{2}".format( i, val, name )

        if remainder:
            print "<output truncated by user via \"display\" arg, {0} more feature weights>".format( remainder )

#############################################################################
# class definition of PearsonFeatureWeights
#############################################################################
class PearsonFeatureWeights( FeatureWeights ):
    """A concrete class that calculates correlation coefficients as well as
    regression parameters for each feature. Features are weighted based on how well
    they linearly correlate (i.e., high Pearson correlation coefficient) with an experimental
    variable.

    An example system where a continuous classifier could be used could be
    would be defining a spectrum of morphology across age or dose response."""

    def __init__( self, name=None, size=None ):
        """Constructor"""
        super( PearsonFeatureWeights, self ).__init__( name=name, size=size )
        if size is not None:
            self.slopes = [None] * size
            self.intercepts = [None] * size
            self.pearson_coeffs = [None] * size
            self.pearson_stderrs = [None] * size
            self.pearson_p_values = [None] * size
            self.spearman_coeffs = [None] * size
            self.spearman_p_values = [None] * size
        else:
            self.slopes = None
            self.intercepts = None
            self.pearson_coeffs = None
            self.pearson_stderrs = None
            self.pearson_p_values = None
            self.spearman_coeffs = None
            self.spearman_p_values = None

    #================================================================
    @classmethod
    def NewFromFeatureSpace( cls, fs ):
        """Calculate regression parameters and correlation statistics that fully define
        a continuous classifier.

        At present the feature weights are proportional the Pearson correlation coefficient
        for each given feature."""

        if fs.normalized_against is None:
            raise ValueError( "Before generating feature weights, call Normalize() of the feature space." )

        from scipy.stats import linregress, spearmanr

        # Known issue: running stats.linregress() with np.seterr (all='raise') has caused
        # arithmetic underflow (FloatingPointError: 'underflow encountered in stdtr' )
        # I think this is something we can safely ignore in this function, and return settings
        # back to normal at the end. -CEC
        np.seterr (under='ignore')

        if fs.name:
            name = cls.__name__ + ' from training set "' + fs.name + '"'
        else:
            name = None

        new_fw = cls( name=name, size=fs.num_features )
        new_fw.associated_feature_space = fs

        #r_val_sum = 0
        r_val_squared_sum = 0
        #r_val_cubed_sum = 0

        ground_truths = np.array( [ float(val) for val in fs.ground_truth_values ] )

        new_fw.feature_names = fs.feature_names[:]

        for feature_index in range( fs.num_features ):

            feature_values = fs.data_matrix[ :, feature_index ]

            slope, intercept, pearson_coeff, p_value, std_err = \
                         linregress( ground_truths, feature_values )

            new_fw.pearson_coeffs[ feature_index ] = pearson_coeff
            new_fw.slopes[ feature_index ] = slope
            new_fw.intercepts[ feature_index ] = intercept
            new_fw.pearson_stderrs[ feature_index ] = std_err
            new_fw.pearson_p_values[ feature_index ] = p_value

            #from math import fabs
            #r_val_sum += fabs(pearson_coeff)
            r_val_squared_sum += pearson_coeff * pearson_coeff
            #r_val_cubed_sum += pearson_coeff * pearson_coeff * pearson_coeff

            try:
                spearman_coeff, spearman_p_val = spearmanr( ground_truths, feature_values )
            except FloatingPointError:
                # to avoid: "FloatingPointError: invalid value encountered in true_divide"
                spearman_coeff, spearman_p_val = ( 0, 1 )

            new_fw.spearman_coeffs[ feature_index ] = spearman_coeff
            new_fw.spearman_p_values[ feature_index ] = spearman_p_val

        #new_fw.values = [ fabs(val) / r_val_sum for val in new_fw.pearson_coeffs ]
        new_fw.values = [val*val / r_val_squared_sum for val in new_fw.pearson_coeffs ]
        #new_fw.values = [val*val*val / r_val_cubed_sum for val in new_fw.pearson_coeffs ]

        # Reset numpy
        np.seterr (all='raise')

        return new_fw

    #================================================================
    def Threshold( self, num_features_to_be_used=None, _all=False, use_spearman=False,
                 min_corr_coeff=None ):
        """Returns a new instance of a PearsonFeatureWeights class derived from this
        instance where the number of features has been reduced to only the top n features,
        where n is specified by the num_features_to_be_used argument.

        If min_corr_coeff is specified, the argument num_features_to_be_used is ignored.

        if _all == True implies simple reorder by feature rank, including 0-weighted features
        If _all == 'nonzero', returns non-zero weighted features ranked by weight."""

        if _all:
            num_features_to_be_used = len( self.values )
        elif num_features_to_be_used is None:
            # Default is top 15% of features
            num_features_to_be_used = int( len( self.values ) * 0.15 )
        elif type( num_features_to_be_used ) is float:
            if num_features_to_be_used <= 0 or num_features_to_be_used > 1:
                raise ValueError('Choose feature reduction fraction on interval (0,1] (got "{0}"'.format( num_features_to_be_used  ) )
            num_features_to_be_used = int( round( num_features_to_be_used * len( self ) ) )
        elif num_features_to_be_used <= 0 or num_features_to_be_used > len( self.values ):
            raise ValueError('Cannot reduce a set of {0} feature weights to requested {1} features.'.\
                                  format( len( self.values ), num_features_to_be_used ) )

        new_weights = self.__class__()
        if self.name:
            if num_features_to_be_used == len( self.feature_names ):
                new_weights.name = self.name + " (rank-ordered)"
            else:
                new_weights.name = self.name + " (top {0} features)".format( num_features_to_be_used )

        if use_spearman:
            abs_corr_coeffs = [ abs( val ) for val in self.spearman_coeffs ]
        else:
            abs_corr_coeffs = [ abs( val ) for val in self.pearson_coeffs ]

        raw_featureweights = zip( abs_corr_coeffs, self.feature_names, self.pearson_coeffs, \
            self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
            self.spearman_coeffs, self.spearman_p_values )

        # sort from max to min
        sorted_featureweights = sorted( raw_featureweights, key=lambda r: r[0], reverse = True )

        # take most correlated features, both positive and negative
        if min_corr_coeff is not None:
            try:
                val = abs( float( min_corr_coeff ) )
            except:
                raise ValueError( 'Cannot convert {0} to a float.'.format( min_corr_coeff ) )
            if val <= 0 or val > 1:
                raise ValueError( 'Abs val of min correlation coefficient must be between 0 and 1.' )

            from itertools import takewhile
            use_these_feature_weights = list( takewhile( lambda x: x[0]>min_corr_coeff, sorted_featureweights ) )
        else:
            from itertools import islice
            use_these_feature_weights = list( islice( sorted_featureweights, num_features_to_be_used ) )
            if _all is not True:
                # You have a problem if any of the features have corellation coefficients of 0
                for i, feat_fig_of_merit in enumerate( [ line[0] for line in use_these_feature_weights ] ):
                    if feat_fig_of_merit == 0:
                        if _all == 'nonzero':
                            print "Features rank index {0} and below have a correllation coefficient of 0. ".format( i )
                            print 'Using {0} features'.format( i )
                            use_these_feature_weights = use_these_feature_weights[ : i ]
                            break
                        else:
                            err_msg = "Can't reduce feature weights \"{0}\" to {1} features. ".format( self.name, num_features_to_be_used )
                            err_msg += "Features ranked {0} and below have a correllation coefficient of 0. ".format( i )
                            err_msg += "Request less features. "
                            raise ValueError( err_msg )

        # we want lists, not tuples!
        abs_corr_coeffs, new_weights.feature_names, new_weights.pearson_coeffs, new_weights.slopes, \
            new_weights.intercepts, new_weights.pearson_stderrs, new_weights.pearson_p_values,\
            new_weights.spearman_coeffs, new_weights. spearman_p_values =\
              [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

        r_val_sum = 0
        for val in abs_corr_coeffs:
            r_val_sum += val * val
        new_weights.values = [ ( (val*val) / r_val_sum ) for val in abs_corr_coeffs ]

        new_weights.associated_feature_space = self.associated_feature_space

        return new_weights

    #================================================================
    def Slice( self, start_index, stop_index ):
        """Return a new instance of PearsonFeatureWeights populated with a
        chunk of middle-ranked features specified by arguments start_index and stop_index."""

        min_index = None
        max_index = None

        if stop_index > start_index:
            min_index = start_index
            max_index = stop_index
        else:
            min_index = stop_index
            max_index = start_index

        if (min_index < 0) or ( max_index > len( self.values ) ):
            raise ValueError( 'Cannot slice, check your start and stop indices.' )

        new_weights = self.__class__()
        if self.name:
            new_weights.name = self.name + " (sliced {0}-{1})".format( min_index, max_index )

        abs_val_pearson_coeffs = [ abs( val ) for val in self.pearson_coeffs ]
        raw_featureweights = zip( self.feature_names, abs_val_pearson_coeffs, self.pearson_coeffs, \
            self.slopes, self.intercepts, self.pearson_stderrs, self.pearson_p_values, \
            self.spearman_coeffs, self.spearman_p_values )

        from itertools import islice
        use_these_feature_weights = list( islice( raw_featureweights, min_index, max_index ) )
 
        new_weights.feature_names, abs_pearson_coeffs, new_weights.pearson_coeffs, new_weights.slopes, \
            new_weights.intercepts, new_weights.pearson_stderrs, new_weights.pearson_p_values,\
            new_weights.spearman_coeffs, new_weights. spearman_p_values =\
              [ list( unzipped_tuple ) for unzipped_tuple in zip( *use_these_feature_weights ) ]

        r_val_sum = 0
        for val in abs_pearson_coeffs:
            r_val_sum += val
        new_weights.values = [ val / r_val_sum for val in abs_pearson_coeffs ]

        new_weights.associated_feature_space = self.associated_feature_space

        return new_weights

    #================================================================
    @output_railroad_switch
    def Print( self, display=None, print_legend=True ):
        """@brief Prints out feature values and statistics"""

        if not display:
            display = len(self)
            remainder = None
        else:
            remainder = len(self) - display

        s = self.__class__.__name__
        if self.name:
            s += ' "{0}:"'.format( self.name )
        s += " ({0} features)".format( len(self ) )
        print s

        print "-----------------------------------"

        if print_legend:
            print "Legend:"
            print "IFW - Feature weight applied to the individual feature"
            print "IPC - Pearson correlation coefficient of feature values vs ground truth"
            print "IPE - Standard Error of IPC"
            print "IPP - P-value of IPC"
            print "ISC - Spearman correlation coefficient of feature values vs ground truth"
            print "IPP - P-value of ISC"
            print ""
        print "NUM\tIFW\tIPC\tIPE\tIPP\tISC\tIPP\tNAME"
        print "===\t===\t===\t===\t===\t===\t===\t===="
        for i in range( display ):
            line_item = "{0}\t".format( i + 1 )
            line_item += "{0:2.4f}\t".format( self.values[i] )
            line_item += "{0:2.4f}\t".format( self.pearson_coeffs[i] )
            line_item += "{0:2.4f}\t".format( self.pearson_stderrs[i] )
            line_item += "{0:2.4f}\t".format( self.pearson_p_values[i] )
            line_item += "{0:2.4f}\t".format( self.spearman_coeffs[i] )
            line_item += "{0:2.4f}\t".format( self.spearman_p_values[i] )
            if len( self.feature_names[i] ) < 50:
                line_item += self.feature_names[i]
            else:
                line_item += self.feature_names[i][:50] + '... (truncated)'
            print line_item
        if remainder:
            print "<output truncated by user, {0} more feature weights>".format( remainder )
