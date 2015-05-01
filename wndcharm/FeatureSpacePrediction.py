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
from .utils import output_railroad_switch
from .FeatureSpace import FeatureSpace
from .FeatureWeights import FeatureWeights, FisherFeatureWeights, PearsonFeatureWeights
from .SingleSamplePrediction import SingleSampleClassification, SingleSampleRegression

#=================================================================================
class FeatureSpacePrediction( object ):
    """Base class container for individual SingleSamplePrediction instances.
    Use for generating predicted values for a FeatureSpace all at once.
    
    The equivalent concept in the C++ implementation is the "split,"
    i.e., command line arg of -n10 would generate 10 train/test splits."""

    #==============================================================
    def __init__( self, training_set=None, test_set=None, feature_weights=None, name=None,
        batch_number=None ):
        """FeatureSpacePrediction constructor"""

        self.training_set = training_set
        self.test_set = test_set
        self.feature_weights = feature_weights
        self.name = name
        if self.name is None and training_set and training_set.name:
            self.name = training_set.name
        self.individual_results = []
        self.tiled_results = None
        self.tiled_ground_truth_values = None
        self.tiled_predicted_values = None

        # Give myself a number so that it looks good when I print out results
        if not batch_number:
            batch_number = self.__class__.obj_count
            self.__class__.obj_count += 1

        self.batch_number = batch_number

        self.num_classifications = 0

        self.std_err = None
        self.predicted_values = None
        self.ground_truth_values = None

        #: If there was randomness associated with generating results
        #: set use_error_bars = True to calculate confidence intervals for
        #: resulting figures of merit
        self.use_error_bars = None

        # This stuff is for correllation analysis
        self.pearson_coeff = None
        self.pearson_p_value = None
        self.pearson_std_err = None
        self.spearman_coeff = None
        self.spearman_p_value = None

        # used for pretty printing
        self.column_headers = None
        self.column_separators = None

    #==============================================================
    def __len__( self ):
        if self.tiled_results:
            return len( self.tiled_results )
        try:
            return len( self.individual_results )
        except:
            return 0

    #==============================================================
    def GenerateStats( self ):
        """Base method; Calculates statistics describing how well predicted values
        in self.individual_results correlate with their ground truth values."""
        #FIXME: Calculate p-value of our standard error figure of merit.

        if self.tiled_results:
            classification_results = self.tiled_results
            ground_truth_values = self.tiled_ground_truth_values
            predicted_values = self.tiled_predicted_values
        else:
            classification_results = self.individual_results
            ground_truth_values = self.ground_truth_values
            predicted_values = self.predicted_values

        if len( classification_results ) == 0:
            raise ValueError( 'No individual results in this batch with which to generate statistics.' )

        self.num_classifications = len( classification_results )

        if ground_truth_values and predicted_values:
            assert( len(ground_truth_values) == len(predicted_values) )

            gt = np.array( ground_truth_values )
            pv = np.array( predicted_values )

            diffs = gt - pv
            diffs = np.square( diffs )
            err_sum = np.sum( diffs )

            from math import sqrt
            self.std_err = sqrt( err_sum / self.num_classifications )

            # no point in doing regression stuff if there's only 1 individual result:
            if len( classification_results ) > 1:
                from scipy.stats import linregress, spearmanr

                # For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
                np.seterr (under='ignore')
                slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
                         linregress( ground_truth_values, predicted_values )

                try:
                    self.spearman_coeff, self.spearman_p_value =\
                   spearmanr( ground_truth_values, predicted_values )
                except FloatingPointError:
                    # to avoid: "FloatingPointError: invalid value encountered in true_divide"
                    self.spearman_coeff, self.spearman_p_value = ( 0, 1 )

                np.seterr (all='raise')
        return self

    #==============================================================
    def NewNFold( self, num_folds=5, *args, **kwargs ):
        """Base method, implemented in daughter classes."""
        raise NotImplementedError()

    #==============================================================
    def RankOrderSort( self ):
        """Rank-order sorts ground truth/predicted value data points for purposes of 
        being graphed."""

        if self.tiled_results:
            classification_results = self.tiled_results
            ground_truth_values = self.tiled_ground_truth_values
            predicted_values = self.tiled_predicted_values
        else:
            classification_results = self.individual_results
            ground_truth_values = self.ground_truth_values
            predicted_values = self.predicted_values

        if not self.ground_truth_values or not self.predicted_values:
            self.GenerateStats()

        # Check again:
        if not ground_truth_values:
            raise ValueError( "Can't rank-order sort: no numeric ground truth values for predicted results." )
        if not predicted_values:
            raise ValueError( "Can't rank-order sort: no sample predicted values" )

        value_pairs = zip( ground_truth_values, predicted_values )

        # sort by ground_truth value first, predicted value second
        sorted_pairs = sorted( sorted( value_pairs, key=lambda x: x[0] ), key=lambda x: x[1] )
        
        # we want lists, not tuples!
        ground_truth_values, predicted_values =\
            [ list( unzipped_tuple ) for unzipped_tuple in zip( *sorted_pairs ) ]

        if self.tiled_results:
            self.tiled_ground_truth_values = ground_truth_values
            self.tiled_predicted_values = predicted_values
        else:
            self.ground_truth_values = ground_truth_values
            self.predicted_values = predicted_values
        return self

    #==============================================================
    def ConfusionMatrix( self ):
        """Only applies for classification problems"""

        outstr = "Confusion Matrix:\n"
        if not self.confusion_matrix:
            outstr += '<empty>'
            return outstr

        # These are the row headers:
        gt_names = sorted( self.confusion_matrix )
        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.confusion_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

        pred_class_names = sorted( pred_class_names )

        # print column headers and horizontal table bar separator
        outstr += ' '*first_col_width + '\t' + "\t".join( pred_class_names ) + '\t|\ttotal\tacc.\n'
        outstr += " "*first_col_width + "\t" + "\t".join( [ '-'*len(name) for name in pred_class_names ] ) + '\t|\t-----\t----\n'

        for gt_class in gt_names:
            gt_row_dict = self.confusion_matrix[ gt_class ]
            outstr += gt_class + '\t'
            for pred_class in pred_class_names:
                if pred_class not in gt_row_dict:
                    outstr += '0\t'
                else:
                    outstr += '{0}\t'.format( gt_row_dict[ pred_class ] )
            outstr += '|\t{0}\t{1:0.2f}%'.format( self.num_classifications_per_class[ gt_class ],
                100.0 * self.num_correct_classifications_per_class[ gt_class ] / self.num_classifications_per_class[ gt_class ] )
            outstr += '\n'
        return outstr

    #==============================================================
    def SimilarityMatrix( self ):
        """Only applies for classification problems
        Doesn't make sense to do this unless the matrix is square"""
        # if row labels == column labels:
        if self.test_set.class_names != self.training_set.class_names:
            return ValueError( "Can't do similarity matrix if self.test_set.class_names != self.training_set.class_names" )

        outstr = "Similarity Matrix:\n"
        if not self.similarity_matrix:
            outstr += '<empty>'
            return outstr

        # These are the row headers:
        gt_names = sorted( self.similarity_matrix )
        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.similarity_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

        pred_class_names = sorted( pred_class_names )

        # print column headers and horizontal table bar separator
        outstr += ' '*first_col_width + '\t' + "\t".join( pred_class_names ) + '\n'
        outstr += " "*first_col_width + "\t" + "\t".join( [ '-'*len(name) for name in pred_class_names ] ) + '\n'

        for gt_class in gt_names:
            gt_row_dict = self.similarity_matrix[ gt_class ]
            outstr += gt_class + '\t'
            for pred_class in pred_class_names:
                if pred_class not in gt_row_dict:
                    outstr += '\t'
                else:
                    outstr += '{0:0.2f}\t'.format( gt_row_dict[ pred_class ] )
            outstr += '\n'
        return outstr

    #==============================================================
    def AvgClassProbMatrix( self ):
        """Only applies for classification problems.
        Meant to be able to be called while results are coming in,
        showing results so far."""

        outstr = "Average Class Probability Matrix:\n"
        if not self.average_class_probability_matrix:
            outstr += '<empty>'
            return outstr

        # These are the row headers:
        gt_names = sorted( self.average_class_probability_matrix )
        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.average_class_probability_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

        pred_class_names = sorted( pred_class_names )

        # print column headers and horizontal table bar separator
        outstr += ' '*first_col_width + '\t' + "\t".join( pred_class_names ) + '\n'
        outstr += " "*first_col_width + "\t" + "\t".join( [ '-'*len(name) for name in pred_class_names ] ) + '\n'

        for gt_class in gt_names:
            gt_row_dict = self.average_class_probability_matrix[ gt_class ]
            outstr += gt_class + '\t'
            for pred_class in pred_class_names:
                if pred_class not in gt_row_dict:
                    outstr += '\t'
                else:
                    outstr += '{0:0.4f}\t'.format( gt_row_dict[ pred_class ] )
            outstr += '\n'
        return outstr

#=================================================================================
class FeatureSpaceClassification( FeatureSpacePrediction ):
    """Container for SingleSampleClassification instances.
    Use for classifying all samples in a FeatureSpace in one sitting."""

    #: obj_count class attribute
    obj_count = 0

    #==============================================================
    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

        super( FeatureSpaceClassification, self ).__init__( *args, **kwargs )

        self.num_correct_classifications = None
        self.classification_accuracy = None

        self.num_classifications_per_class = None
        self.num_correct_classifications_per_class = None

        self.confusion_matrix = None
        self.similarity_matrix = None
        self.average_class_probability_matrix = None

    #==============================================================    
    def __str__( self ):
        outstr = '<' + self.__class__.__name__
        if self.batch_number is not None:
            outstr += ' #' + str( self.batch_number )
        if self.individual_results:
            outstr += ' n=' + str( len( self.individual_results ) )
        if self.num_correct_classifications:
            outstr += ' n_corr=' + str( self.num_correct_classifications )
        if self.classification_accuracy is not None:
            outstr += ' acc={0:0.2f}%'.format( self.classification_accuracy * 100 )
        if self.std_err is not None:
            outstr += ' std_err={0:0.2f}'.format( self.std_err )

        return outstr + '>'
    #==============================================================
    def __repr__( self ):
        return str(self)

    #==============================================================
    def GenerateStats( self ):
        """Fills out the confusion matrix, similarity matrix, and class probability matrix."""

        # Run a standard error analysis if ground_truth/predicted vals exist (call base method):
        # This also sets self.num_classifications
        super( FeatureSpaceClassification, self ).GenerateStats()

        num_classes = self.training_set.num_classes

        # Initialize the matrices:

        # Remember! Dicts are not guaranteed to maintain key order but lists are
        # When cycling through the matrix, iterate over the lists, and not the keys of the dict
        from collections import defaultdict # introduced Python 2.5

        # These are dicts of dicts in the form:
        # self.confusion_matrix[ <Ground Truth Class> ][ <Predicted Class> ] == count
        self.confusion_matrix = defaultdict( lambda: defaultdict( int ) )
        self.average_class_probability_matrix = defaultdict( lambda: defaultdict( float ) )

        self.num_correct_classifications = 0

        self.num_classifications_per_class = defaultdict( int )
        self.num_correct_classifications_per_class = defaultdict( int )

        classification_results = self.tiled_results if self.tiled_results else self.individual_results

        for indiv_result in classification_results:
            gt_class = indiv_result.ground_truth_class_name
            if gt_class == None:
                gt_class = "UNKNOWN"
            pred_class = indiv_result.predicted_class_name

            self.num_classifications_per_class[ gt_class ] += 1

            if gt_class == pred_class:
                self.num_correct_classifications += 1
                self.num_correct_classifications_per_class[ gt_class ] += 1
            
            self.confusion_matrix[ gt_class ][ pred_class ] += 1

            # FIXME: is there any possibility that the order of the values in the marginal
            # probability array don't correspond with the order of the training set classes?
            if indiv_result.marginal_probabilities is not None:
                mp = zip( self.training_set.class_names, indiv_result.marginal_probabilities )
                for putative_class, marg_prob in mp:
                    self.average_class_probability_matrix[ gt_class ][ putative_class ] += marg_prob

        # Finalize the Average Class Probability Matrix by dividing each marginal probability
        # sum by the number of marginal probabilities for that ground truth:
        for row in self.test_set.class_names:
            for col in self.training_set.class_names:
                self.average_class_probability_matrix[ row ][ col ] /= \
                   self.num_classifications_per_class[ row ]

        # The similarity matrix is just the average class probability matrix
        # normalized to have 1's in the diagonal.
        # Doesn't make sense to do this unless the matrix is square
        # if row labels == column labels:
        if self.test_set.class_names == self.training_set.class_names:
            from copy import deepcopy
            self.similarity_matrix = deepcopy( self.average_class_probability_matrix )
            for row in self.test_set.class_names:
                denom = self.similarity_matrix[ row ][ row ]
                for col in self.training_set.class_names:
                    self.similarity_matrix[ row ][ col ] /= denom

        self.classification_accuracy = float( self.num_correct_classifications) / float( self.num_classifications )
        return self

    #==============================================================
    @output_railroad_switch
    def Print( self ):
        """Prints out statistics from this batch's image classifications, which include 
        classification accuracy, confusion matrix, similarity matrix, and average class 
        probability matrix."""

        if self.classification_accuracy == None:
            self.GenerateStats()

        classification_results = self.tiled_results if self.tiled_results else self.individual_results

        print '='*50
        s = self.__class__.__name__
        if self.name:
            s += ' "' + self.name + '"'
        s += " (" + str( len( classification_results ) ) + " classifications)"
        print s

        acc = self.classification_accuracy
        n = self.num_classifications
        n_correct = self.num_correct_classifications

        if not self.use_error_bars:
            print "{0}/{1} correct = {2:0.2f}%".format( n_correct, n, acc * 100 )
        else:
            # Using either normal approximation of binomial distribution or the Wilson score interval
            # to calculate standard error of the mean, depending on the situation.
            # For more info, see http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
            # The confidence interval is S.E.M. * quantile for your chosen accuracy
            # The quantile for 95% accuracy is ~ 1.96.
            z = 1.95996;
            z2 = 3.84144 # z^2

            from math import sqrt

            # This is a rule of thumb test to check whecther sample size is large enough
            # to use normal approximation of binomial distribution:
            if ((n * acc) > 5) and ((n * (1 - acc)) > 5):
                # Using normal approximation:
                std_error_of_mean = sqrt( acc * (1-acc) / n )
                conf_interval = z * std_error_of_mean
                print "{0}/{1} correct = {2:0.2f} +/- {3:0.2f}% w/ 95% conf. (normal approx. interval)".format(
                    n_correct, n, acc * 100, conf_interval * 100 )
            else:
                # Using Wilson approximation:
                # This term goes to 1 as number of classifications gets large:
                coeff = 1 / (1+(z2/n))
                raw_acc = acc
                # Wilson accuracy modifies the raw accuracy for low n:
                acc = coeff * (raw_acc + z2/(2*n))
                conf_interval = coeff * z * sqrt( (raw_acc*(1-raw_acc)/n) + (z2/(4*n**2)) )

                outstr = "{0}/{1} correct = {2:0.1f}% raw accuracy".format(
                    n_correct, n, raw_acc * 100, conf_interval * 100 )
                outstr += " ({0:0.2f} +/- {1:0.2f}% w/ 95% conf. (Wilson score interval))".format(
                        acc * 100, conf_interval * 100)
                print outstr

        if self.std_err is not None:
            print "Standard Error: {0:0.4f}".format( self.std_err)
        if self.pearson_coeff is not None:
            print "Pearson Coefficient: {0:0.4f}".format( self.pearson_coeff )
        if self.spearman_coeff is not None:
            print "Spearman Coefficient: {0:0.4f}".format( self.spearman_coeff )
        print "\n"

        print self.ConfusionMatrix()
        print self.SimilarityMatrix()
        print self.AvgClassProbMatrix()
    
    #==============================================================
    @output_railroad_switch
    def PrintDistanceMatrix( self, method='mps' ):
        """Generate a distance matrix indicating similarity of image classes to one another.
        The output can be fed into the PHYLIP or other program to create a dendrogram.

        Methods can be {'max', 'mean', 'top', 'bottom', 'mps'}.
        
        To generate a distance matrix using marginal probability space ('mps'), it is legal for
        the test set and training set to have a different number of classes, since the coefficients
        in the average class probability matrix are used as coordinates, but otherwise
        the the input (similarity) matrix must be square."""

        # The source matrix is a numpy matrix of coefficients from which the 
        # distances will be derived
        source_matrix = np.empty( (self.test_set.num_classes, self.test_set.num_classes) ) #initialize
        row_index = 0
        for row in self.test_set.class_names:
            col_index = 0
            for col in self.training_set.class_names:
                source_matrix[ row_index ][ col_index ] = self.average_class_probability_matrix[row][col]
                col_index += 1
            row_index += 1

        output_matrix = {}

        # initialize the rows
        for test_set_class_name in self.test_set.class_names:
            output_matrix[ test_set_class_name ] = {}
        # now the columns:
        for training_set_class_name in self.training_set.class_names:
            for test_set_class_name in self.test_set.class_names:
                output_matrix[ test_set_class_name ][ training_set_class_name ] = 0

        if method == 'mps':
            if self.average_class_probability_matrix == None:
                self.GenerateStats()

            for row in range( self.test_set.num_classes ):
                for col in range( self.test_set.num_classes ):
                    if row == col:
                        output_matrix[ row, col ] = 0
                    else:
                        output_matrix[ row, col ] = np.sum( np.square( source_matrix[row]-source_matrix[col] ) )
        else:
            if self.similarity_matrix == None:
                self.GenerateStats()
            # Have to check again, since for some experiments it doesn't make sense to have a 
            # similarity matrix, e.g., when test set ground truth is not known.
            if self.similarity_matrix == None:
                raise ValueError( 'Similarity matrix is not possible with this data set. ' + \
                                  'Possible reason for this is that your test set does not have ' + \
                                  'ground truth defined.' )
            if method == 'max':
                raise NotImplementedError
            elif method == 'mean':
                raise NotImplementedError
            elif method == 'top':
                raise NotImplementedError
            elif method == 'bottom':
                raise NotImplementedError
            else:
                raise ValueError( "Distance matrix method {0} not recognized. Valid options are ".format(\
                                                    method ) + "'max', 'mean', 'top', 'bottom', 'mps'" )

        #column_headers = "\t".join( self.test_set.class_names )
        #column_headers += "\n"
        #column_headers += "\t".join( [ '-'*len(name) for name in self.test_set.class_names ] )
        #print "Distance Matrix (method = '{0}'):".format( method )
        #print column_headers
        print self.test_set.num_classes
        for row in range( self.test_set.num_classes ):
            line = "{0}\t".format( self.test_set.class_names[ row ] )
            for col in range( self.test_set.num_classes ):
                line += '{0:0.4f}\t'.format( output_matrix[ row, col ] )
            print line

    #==============================================================
    @classmethod
    @output_railroad_switch
    def NewWND5( cls, training_set, test_set, feature_weights, name=None, batch_number=None,
                quiet=False, norm_factor_threshold=None, error_bars=False ):
        """The equivalent of the "wndcharm classify" command in the command line implementation
        of WND-CHARM. Input a training set, a test set, and feature weights, and returns a
        new instance of a FeatureSpaceClassification, with self.individual_results
        filled with a new instances of SingleSampleClassification.

        FIXME: What happens when the ground truth is not known? Currently they would all be shoved
        into class 1, might not be a big deal since class name should be something
        like "UNKNOWN"
        """

        # type checking
        if not isinstance( training_set, FeatureSpace ):
            raise ValueError( 'First argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
        if not isinstance( test_set, FeatureSpace ):
            raise ValueError( 'Second argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
        if not isinstance( feature_weights, FeatureWeights ):
            raise ValueError( 'Third argument to New must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )
    
        # feature comparison
        if test_set.feature_names != feature_weights.feature_names:
            raise ValueError( "Can't classify, features in test set don't match features in weights. Try translating feature names from old style to new, or performing a FeatureReduce()" )
        if test_set.feature_names != training_set.feature_names:
            raise ValueError( "Can't classify, features in test set don't match features in training set. Try translating feature names from old style to new, or performing a FeatureReduce()" )

        np.seterr( under='ignore' )

        train_set_len = len( training_set.feature_names )
        test_set_len = len( test_set.feature_names )
        feature_weights_len = len( feature_weights.feature_names )

        # instantiate myself
        batch_result = cls( training_set, test_set, feature_weights, name, batch_number )
        batch_result.use_error_bars = error_bars

        # Are the samples to be classified tiled?
        if test_set.num_samples_per_group > 1:
            batch_result.tiled_results = []
            batch_result.tiled_ground_truth_values = []
            batch_result.tiled_predicted_values = []

        # Say what we're going to do
        if not quiet:
            print "Classifying test set '{0}' ({1} features) against training set '{2}' ({3} features)".\
                    format( test_set.name, test_set_len, training_set.name, train_set_len )
            if batch_result.tiled_results:
                print "Performing tiled classification."
            column_header = "image\tnorm. fact.\t"
            column_header +=\
                "".join( [ "p(" + class_name + ")\t" for class_name in training_set.class_names ] )
            column_header += "act. class\tpred. class\tpred. val."
            print column_header

        # Will there be a numeric predicted value associated with this classification?
        train_set_interp_coeffs = None
        if training_set.interpolation_coefficients != None and len( training_set.interpolation_coefficients) != 0:
            train_set_interp_coeffs = np.array( training_set.interpolation_coefficients )
            batch_result.predicted_values = []

        # Are there numeric ground truth values associated with the input samples?
        test_set_interp_coeffs = None
        if test_set.interpolation_coefficients != None and len( test_set.interpolation_coefficients ) != 0:
            test_set_interp_coeffs = np.array( test_set.interpolation_coefficients )
            batch_result.ground_truth_values = []

        for test_class_index in range( test_set.num_classes ):
            num_class_imgs, num_class_features = test_set.data_list[ test_class_index ].shape

            # Get tiling ready if need be
            if test_set.num_samples_per_group > 1:
                tile_results_in_this_sample_group = []

            for test_image_index in range( num_class_imgs ):
                one_image_features = test_set.data_list[ test_class_index ][ test_image_index,: ]
                result = SingleSampleClassification._WND5( training_set, one_image_features, feature_weights.values )
                
                if norm_factor_threshold and (result.normalization_factor > norm_factor_threshold):
                    continue
                result.source_filepath = test_set.sample_names[ test_class_index ][ test_image_index ]
                result.ground_truth_class_name = test_set.class_names[ test_class_index ]
                result.batch_number = batch_number
                result.name = name
                if result.marginal_probabilities:
                    # Sometimes the result comes back with a non-call, like when the sample image
                    # collides with every test image
                    marg_probs = np.array( result.marginal_probabilities )
                    result.predicted_class_name = training_set.class_names[ marg_probs.argmax() ]
                else:
                    marg_probs = np.array( [np.nan] * training_set.num_classes )
                    result.predicted_class_name = "*Collided with every training sample*"

                # interpolated value, if applicable
                if train_set_interp_coeffs is not None:
                    interp_val = np.sum( marg_probs * train_set_interp_coeffs )
                    result.predicted_value = interp_val
                    batch_result.predicted_values.append( interp_val )

                if test_set_interp_coeffs is not None:
                    result.ground_truth_value = test_set_interp_coeffs[ test_class_index ]
                    batch_result.ground_truth_values.append( test_set_interp_coeffs[ test_class_index ] )

                # Helps to identify which results correspond with which sample
                result.samplegroupid = test_set.sample_group_ids[ test_class_index ][ test_image_index ]
                result.samplesequenceid = test_set.sample_sequence_ids[ test_class_index ][ test_image_index ]

                if test_set.num_samples_per_group > 1:
                    result.tile_index = len( tile_results_in_this_sample_group )
                    result.num_samples_in_group = test_set.num_samples_per_group

                if not quiet:
                    result.Print( line_item=True )
                batch_result.individual_results.append( result )

                # TILING SECTION:
                # Create a whole image classification result that
                # is the average of all the calls from all the tiles
                if test_set.num_samples_per_group > 1:
                    tile_results_in_this_sample_group.append( result )
                    if len( tile_results_in_this_sample_group ) >= test_set.num_samples_per_group:
                        aggregated_result = SingleSampleClassification()
                        # Use the last result from above
                        aggregated_result.source_filepath = result.source_filepath
                        aggregated_result.tile_index = 'AVG'
                        aggregated_result.ground_truth_class_name = result.ground_truth_class_name
                        marg_prob_lists = [ [] for i in xrange( training_set.num_classes ) ]
                        norm_factor_list = []
                        for tile_result in tile_results_in_this_sample_group:
                            if tile_result.marginal_probabilities:
                                # Sometimes the result comes back with a non-call, like when the sample image
                                # collides with every test image
                                for class_index, val in enumerate( tile_result.marginal_probabilities ):
                                    marg_prob_lists[ class_index ].append( val )
                                norm_factor_list.append( tile_result.normalization_factor )

                        if any( [ len( class_marg_prob_list ) <= 0 for class_marg_prob_list in marg_prob_lists ] ):
                            aggregated_result.marginal_probabilities = np.array( [np.nan] * training_set.num_classes )
                            aggregated_result.predicted_class_name = "*Collided with every training image*"
                        else:
                            marg_probs = [ sum( l ) / float( len( l ) ) for l in marg_prob_lists ]
                            aggregated_result.marginal_probabilities = marg_probs
                            np_marg_probs = np.array( marg_probs )
                            aggregated_result.predicted_class_name = training_set.class_names[ np_marg_probs.argmax() ]
                            aggregated_result.normalization_factor = sum(norm_factor_list)/float(len(norm_factor_list))

                            # interpolated value, if applicable
                            if train_set_interp_coeffs is not None:
                                aggregated_result.predicted_value = np.sum( np_marg_probs * train_set_interp_coeffs )
                                batch_result.tiled_predicted_values.append( aggregated_result.predicted_value )
                            if test_set_interp_coeffs is not None:
                                aggregated_result.ground_truth_value = test_set_interp_coeffs[ test_class_index ]
                                batch_result.tiled_ground_truth_values.append( aggregated_result.ground_truth_value )
                        # Save references to the tiled results on the aggregated result
                        # To facilitate creation of graphs, heatmaps, etc
                        aggregated_result.tiled_results = tile_results_in_this_sample_group
                        batch_result.tiled_results.append( aggregated_result )
                        if not quiet:
                            aggregated_result.Print( line_item=True )

                        # now reset for the next sample group
                        tile_results_in_this_sample_group = []
                        tile_count = 0
        np.seterr (all='raise')
        return batch_result

#=================================================================================
class FeatureSpaceRegression( FeatureSpacePrediction ):
    """Container for SingleSampleRegression instances.
    Use for regressing all samples in a FeatureSpace in one sitting."""

    #: obj_count class attribute
    obj_count = 0

    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

        super( FeatureSpaceRegression, self ).__init__( *args, **kwargs )
        self.predicted_values = []

    #==============================================================    
    def __str__( self ):
        outstr = '<' + self.__class__.__name__
        if self.batch_number is not None:
            outstr += ' #' + str( self.batch_number )
        if self.individual_results:
            outstr += ' n=' + str( len( self.individual_results ) )
        if self.std_err is not None:
            outstr += ' R={0:0.2f}'.format( self.std_err )
        return outstr + '>'
    #==============================================================
    def __repr__( self ):
        return str(self)

    #=====================================================================
    def GenerateStats( self ):
        """Calls base class method to run a standard error analysis if ground_truth/
        predicted vals exist."""
        super( FeatureSpaceRegression, self ).GenerateStats()
        return self

    #=====================================================================
    @output_railroad_switch
    def Print( self ):
        """Calculates and outputs batch-level statistics based on the
        SingleSampleRegressions contained in self.individualresults."""
        if self.std_err == None:
            self.GenerateStats()

        print "==========================================="
        print "Number of observations: {0}".format( self.num_classifications )
        if self.std_err != None:
            print "Standard error of predicted vs. ground truth values: {0}".format( self.std_err )
        #print "p-value for this split: {0}".format( self.p_value )

    #=====================================================================
    @classmethod
    def NewMultivariateLinear( cls, test_set, feature_weights, name=None, batch_number=None,
            quiet=False ):
        """Uses Pearson-coefficient weighted Multivatiate Linear classifier."""

        # type checking
        if not isinstance( test_set, FeatureSpace ):
            raise ValueError( 'First argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )    
        if not isinstance( feature_weights, PearsonFeatureWeights ):
            raise ValueError( 'Second argument to New must be of type "PearsonFeatureWeights", you gave a {0}'.format( type( feature_weights ).__name__ ) )

        # feature comparison
        if test_set.feature_names != feature_weights.feature_names:
            raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )

        # say what we're gonna do
        if not quiet:
            out_str = 'Classifying test set "{0}" ({1} images, {2} features)\n\tagainst training set "{3}" ({4} images)'
            print out_str.format( test_set.name, test_set.num_samples, \
              len( test_set.feature_names ), feature_weights.associated_feature_space.name, \
              feature_weights.associated_feature_space.num_samples )

        if not quiet:
            column_header = "image\tground truth\tpred. val."
            print column_header

        batch_result = cls( feature_weights.associated_feature_space, test_set, feature_weights,
                name, batch_number )

        if test_set.ground_truth_values is not None and \
                len( test_set.ground_truth_values ) != 0:
            batch_result.ground_truth_values = test_set.ground_truth_values

        for test_image_index in range( test_set.num_samples ):
            one_image_features = test_set.data_matrix[ test_image_index,: ]
            result = SingleSampleRegression._MultivariateLinear( one_image_features, feature_weights )
            result.batch_number = batch_result.batch_number
            result.name = name
            result.source_filepath = test_set.sample_names[ test_image_index ]
            result.ground_truth_value = test_set.ground_truth_values[ test_image_index ]
            batch_result.predicted_values.append( result.predicted_value )

            if not quiet:
                result.Print( line_item = True )
            batch_result.individual_results.append( result )

        batch_result.GenerateStats()
        return batch_result

    #=====================================================================
    @classmethod
    def NewLeastSquares( cls, training_set, test_set, feature_weights, name=None,
            batch_number=None, leave_one_out=False, quiet=False ):
        """Uses Linear Least Squares Regression classifier in a feature space filtered/weighed
        by Pearson coefficients.
        
        Use cases:
        1. if training_set != test_set and both not none: straight classification
        2. if training_set is not None and test_set is None = Leave one out cross validation
        3. if training_set == test_set == not None: Shuffle/Split cross validation
        """

        # Type checking
        # First arg must be a FeatureSpace
        if not isinstance( training_set, FeatureSpace ):
            raise ValueError( 'First argument to NewLeastSquares must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
        # Second arg must be a FeatureSpace or a None
        if (not isinstance( test_set, FeatureSpace )) and (test_set is not None):
            raise ValueError( 'Second argument to NewLeastSquares must be of type "FeatureSpace" or None, you gave a {0}'.format( type( test_set ).__name__ ) )
        if not isinstance( feature_weights, PearsonFeatureWeights ):
            raise ValueError( 'Third argument to New must be of type "PearsonFeatureWeights", you gave a {0}'.format( type( feature_weights ).__name__ ) )

        # If there's both a training_set and a test_set, they both have to have the same features
        if training_set and test_set:
            if training_set.feature_names != training_set.feature_names:
                raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )
        # Check feature_weights
        if training_set.feature_names != feature_weights.feature_names:
            raise ValueError("Can't classify, features don't match. Try a FeatureReduce()" )

        # Least squares regression requires a featres matrix augmented with ones
        # signifying the constant, i.e., the y-intercept
        # Make a copy of the feature sets to augment while leaving original matrices unchanged.
        from copy import deepcopy
        augmented_train_set = deepcopy( training_set )

        # figure out what we're gonna do
        if training_set and not test_set:
            leave_one_out = True
            cross_validation = True
            augmented_test_set = augmented_train_set 
        elif training_set is test_set:
            cross_validation = True
            augmented_test_set = augmented_train_set
        else:
            augmented_test_set = deepcopy( test_set )
            cross_validation = False

        # say what we're gonna do
        if not quiet:
            if cross_validation:
                out_str = 'Cross validation of training set "{0}" ({1} images, {2} features)'.format(
                      training_set.name, training_set.num_samples, len( training_set.feature_names ) )

                out_str += "\nWITH"
                if not leave_one_out:
                    out_str += "OUT"
                out_str += " using leave-one-out analysis.\n"
            
            else:
                out_str = 'Classifying test set "{0}" ({1} images, {2} features)'.format(
                      test_set.name, test_set.num_samples, len( test_set.feature_names ) )
                out_str += '\n\tagainst training set "{0}" ({1} images)'.format(
                            training_set.name, training_set.num_samples )
            print out_str

        # Now, build the augmented feature matrices, which includes multiplying the feature
        # space by the weights, and augmenting the matrices with 1's

        oldsettings = np.seterr(all='ignore')

        augmented_train_set.data_matrix *= feature_weights.values
        augmented_train_set.data_matrix = np.hstack( 
              [ augmented_train_set.data_matrix, np.ones( ( augmented_train_set.num_samples, 1 ) ) ] )
        augmented_train_set.num_features += 1 # Tell the object it has a new feature column
        if not cross_validation:
            augmented_test_set.data_matrix *= feature_weights.values
            augmented_test_set.data_matrix = np.hstack( 
              [ augmented_test_set.data_matrix, np.ones( ( augmented_test_set.num_samples, 1 ) ) ] )
            augmented_test_set.num_features += 1 # Tell the object it has a new feature column

        if not quiet:
            print "image\tground truth\tpred. val."

        batch_result = cls( training_set, test_set, feature_weights, name, batch_number )

        if augmented_test_set.ground_truth_values is not None and \
                len( augmented_test_set.ground_truth_values ) != 0:
            batch_result.ground_truth_values = augmented_test_set.ground_truth_values

        intermediate_train_set = augmented_train_set
        for test_image_index in range( augmented_test_set.num_samples ):
            if leave_one_out:
                sample_group_id = augmented_train_set.sample_group_ids[ test_image_index ]
                intermediate_train_set = augmented_train_set.SampleReduce(\
                        leave_out_sample_group_ids = sample_group_id, override=True)
            one_image_features = augmented_test_set.data_matrix[ test_image_index,: ]
            result = SingleSampleRegression._LeastSquares( \
                           one_image_features, intermediate_train_set )
            result.batch_number = batch_result.batch_number
            result.name = name
            result.source_filepath = augmented_test_set.sample_names[ test_image_index ]
            result.ground_truth_value = augmented_test_set.ground_truth_values[ test_image_index ]
            batch_result.predicted_values.append( result.predicted_value )

            if not quiet:
                result.Print( line_item = True )
            batch_result.individual_results.append( result )

        # return settings to original
        np.seterr(**oldsettings)

        batch_result.GenerateStats()
        return batch_result
