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
from .SingleSamplePrediction import SingleSampleClassification, SingleSampleRegression,\
    AveragedSingleSamplePrediction

#=================================================================================
class _FeatureSpacePrediction( object ):
    """Base class container for individual SingleSamplePrediction instances.
    Use for generating predicted values for a FeatureSpace all at once.
    
    The equivalent concept in the C++ implementation is the "split,"
    i.e., command line arg of -n10 would generate 10 train/test splits."""

    #==============================================================
    def __init__( self, training_set=None, test_set=None, feature_weights=None, name=None,
        split_number=None ):
        """_FeatureSpacePrediction constructor"""

        self.training_set = training_set
        self.test_set = test_set
        self.feature_weights = feature_weights
        self.name = name
        if self.name is None and training_set and training_set.name:
            self.name = training_set.name
        self.individual_results = []
        self.averaged_results = None
        self.averaged_ground_truth_values = None
        self.averaged_predicted_values = None

        # Give myself a number so that it looks good when I print out results
        if not split_number:
            split_number = self.__class__.obj_count
            self.__class__.obj_count += 1

        self.split_number = split_number

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
        if self.averaged_results:
            return len( self.averaged_results )
        try:
            return len( self.individual_results )
        except:
            return 0

    #==============================================================
    def GenerateStats( self ):
        """Base method; Calculates statistics describing how well predicted values
        in self.individual_results correlate with their ground truth values."""
        #FIXME: Calculate p-value of our standard error figure of merit.

        if self.averaged_results:
            classification_results = self.averaged_results
            ground_truth_values = self.averaged_ground_truth_values
            predicted_values = self.averaged_predicted_values
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

                np.seterr( all='raise' )
        return self

#    #==============================================================
#    def NewNFold( self, num_folds=5, *args, **kwargs ):
#        """Base method, implemented in daughter classes."""
#        raise NotImplementedError()

    #==============================================================
    def RankOrderSort( self, use_averaged_results=True ):
        """Rank-order sorts ground truth/predicted value data points for purposes of 
        being graphed.

        use_averaged_results - bool - If this object has averaged results (due to tiling or
            "per sample" aggregation across splits, use those results instead of
            individual results.

        returns - (list, list) - ground truth values, predicted values"""

        self.GenerateStats()

        if use_averaged_results and self.averaged_results:
            classification_results = self.averaged_results
            ground_truth_values = self.averaged_ground_truth_values
            predicted_values = self.averaged_predicted_values
        else:
            classification_results = self.individual_results
            ground_truth_values = self.ground_truth_values
            predicted_values = self.predicted_values
            
        # Check again:
        if not ground_truth_values:
            raise ValueError( "Can't rank-order sort: no numeric ground truth values for predicted results." )
        if not predicted_values:
            raise ValueError( "Can't rank-order sort: no sample predicted values" )

        value_pairs = zip( ground_truth_values, predicted_values )

        # sort by ground_truth value first, predicted value second
        value_pairs.sort( key=lambda x: x[0] )
        value_pairs.sort( key=lambda x: x[1] )
        
        # we want lists, not tuples!
        ground_truth_values, predicted_values =\
            [ list( unzipped_tuple ) for unzipped_tuple in zip( *value_pairs ) ]

        if use_averaged_results and self.averaged_results:
            self.averaged_ground_truth_values = ground_truth_values
            self.averaged_predicted_values = predicted_values
        else:
            self.ground_truth_values = ground_truth_values
            self.predicted_values = predicted_values

        return ground_truth_values, predicted_values

    #==============================================================
    def ConfusionMatrix( self ):
        """Only applies for classification problems"""

        outstr = "Confusion Matrix:\n"
        if not self.confusion_matrix:
            outstr += '<empty>'
            return outstr

        # These are the row headers:

        # do your best to sort the column and row headers by numeric value
        if self.test_set.interpolation_coefficients is not None:
            gt_coeffs = self.test_set.interpolation_coefficients
            gt_class_names = self.test_set.class_names
            gt_lut = { name : val for name, val in zip( gt_class_names, gt_coeffs ) }
            gt_sortfunc = lambda n: gt_lut[n]
        else:
            gt_sortfunc = lambda n: n
        gt_names = sorted( self.confusion_matrix.keys(), key=gt_sortfunc)

        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.confusion_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

        if self.training_set.interpolation_coefficients is not None:
            pred_coeffs = self.training_set.interpolation_coefficients
            pred_class_names = self.training_set.class_names
            pred_lut = { name : val for name, val in zip( pred_class_names, pred_coeffs ) }
            pred_sortfunc = lambda n: pred_lut[n]
        else:
            pred_sortfunc = lambda n: n

        pred_class_names = sorted( pred_class_names, key=pred_sortfunc )

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
        if self.test_set.interpolation_coefficients is not None:
            gt_coeffs = self.test_set.interpolation_coefficients
            gt_class_names = self.test_set.class_names
            gt_lut = { name : val for name, val in zip( gt_class_names, gt_coeffs ) }
            gt_sortfunc = lambda n: gt_lut[n]
        else:
            gt_sortfunc = lambda n: n
        gt_names = sorted( self.similarity_matrix.keys(), key=gt_sortfunc)
        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.similarity_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

#        if self.training_set.interpolation_coefficients is not None:
#            pred_coeffs = self.training_set.interpolation_coefficients
#            pred_class_names = self.training_set.class_names
#            pred_lut = { name : val for name, val in zip( pred_class_names, pred_coeffs ) }
#            pred_sortfunc = lambda n: pred_lut[n]
#        else:
#            pred_sortfunc = lambda n: n

        pred_class_names = sorted( pred_class_names, key=gt_sortfunc )

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
        # do your best to sort the column and row headers by numeric value
        if self.test_set.interpolation_coefficients is not None:
            gt_coeffs = self.test_set.interpolation_coefficients
            gt_class_names = self.test_set.class_names
            gt_lut = { name : val for name, val in zip( gt_class_names, gt_coeffs ) }
            gt_sortfunc = lambda n: gt_lut[n]
        else:
            gt_sortfunc = lambda n: n
        gt_names = sorted( self.average_class_probability_matrix.keys(), key=gt_sortfunc)
        first_col_width = max( [ len( class_name ) for class_name in gt_names ] )

        # Collect all the predicted value class names that have appeared this far:
        # These become the column headers.
        pred_class_names = set()
        for gt_class in gt_names:
            gt_row_dict = self.average_class_probability_matrix[ gt_class ]
            for pred_class in gt_row_dict.keys():
                pred_class_names.add( pred_class )

        if self.training_set.interpolation_coefficients is not None:
            pred_coeffs = self.training_set.interpolation_coefficients
            pred_class_names = self.training_set.class_names
            pred_lut = { name : val for name, val in zip( pred_class_names, pred_coeffs ) }
            pred_sortfunc = lambda n: pred_lut[n]
        else:
            pred_sortfunc = lambda n: n

        pred_class_names = sorted( pred_class_names, key=pred_sortfunc )

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
class FeatureSpaceClassification( _FeatureSpacePrediction ):
    """Container for SingleSampleClassification instances.
    Use for classifying all samples in a FeatureSpace in one sitting."""

    #: obj_count class attribute
    obj_count = 0

    #==============================================================
    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, split_number=None"""

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
        if self.split_number is not None:
            outstr += ' #' + str( self.split_number )
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

        classification_results = self.averaged_results if self.averaged_results else self.individual_results

        for indiv_result in classification_results:
            gt_class = indiv_result.ground_truth_label
            if gt_class == None:
                gt_class = "UNKNOWN"
            pred_class = indiv_result.predicted_label

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

        classification_results = self.averaged_results if self.averaged_results else self.individual_results

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
    def NewWND5( cls, training_set, test_set, feature_weights=None, name=None, split_number=None,
                quiet=False, norm_factor_threshold=None, error_bars=False ):
        """The equivalent of the "wndcharm classify" command in the command line implementation
        of WND-CHARM. Input a training set, a test set, and feature weights, and returns a
        new instance of a FeatureSpaceClassification, with self.individual_results
        filled with a new instances of SingleSampleClassification.

        If feature_weights == None: use 1's as weights."""

        # type checking
        if not isinstance( training_set, FeatureSpace ):
            raise ValueError( 'First argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
        if not isinstance( test_set, FeatureSpace ):
            raise ValueError( 'Second argument to New must be of type "FeatureSpace", you gave a {0}'.format( type( test_set ).__name__ ) )
        if feature_weights is not None and not isinstance( feature_weights, FeatureWeights ):
                raise ValueError( 'Third argument to New must be of type "FeatureWeights" or derived class, you gave a {0}'.format( type( feature_weights ).__name__ ) )

        # feature comparison
        if test_set.feature_names != training_set.feature_names:
            raise ValueError( "Can't classify, features in test set don't match features in training set. Try translating feature names from old style to new, or performing a FeatureReduce()" )
        if feature_weights is not None and test_set.feature_names != feature_weights.feature_names:
            raise ValueError( "Can't classify, features in test set don't match features in weights. Try translating feature names from old style to new, or performing a FeatureReduce()" )

        # ignore divides => dist matrices with 0's down the diagonal will have nan's
        # which should be ahndled by the mask
        np.seterr( under='ignore', divide='ignore' )

        n_feats = len( training_set.feature_names )
        if feature_weights is None:
            feature_weights = FisherFeatureWeights( size=n_feats )
            feature_weights.values = np.ones( (n_feats,) )
            feature_weights.feature_names = training_set.feature_names

        # instantiate myself
        split_result = cls( training_set, test_set, feature_weights, name, split_number )
        split_result.use_error_bars = error_bars

        # Say what we're going to do
        if not quiet:
            print "Classifying test set '{0}' ({1} samples) against training set '{2}' ({3} samples)".\
                    format( test_set.name, test_set.num_samples, training_set.name, training_set.num_samples )
            if test_set.num_samples_per_group > 1:
                print "Performing tiled classification."

        # Any collisions? (i.e., where the distance from test sample to training sample
        # is below machine epsilon, i.e., 2.2204e-16
        epsilon = np.finfo( np.float ).eps
        import numpy.ma as ma

        # Create slicer for training set class boundaries
        slice_list = []
        start_index = 0
        for n_class_train_samps in training_set.class_sizes:
            end_index = start_index + n_class_train_samps
            slice_list.append( slice( start_index, end_index ) )
            start_index = end_index

        # Create distance matrix:
        # result dist_mat where rows => train samps and cols => test_samps
        wts = np.array( feature_weights.values )
        w_train_featspace = training_set.data_matrix * wts
        if training_set is test_set:
            from scipy.spatial.distance import pdist
            from scipy.spatial.distance import squareform
            # collisions defined by L1 norm, aka taxicab, aka Manhattan, aka cityblock
            # in UNWEIGHTED FEATURE SPACE
            L1_dists = np.absolute( squareform( pdist( training_set.data_matrix, 'cityblock' ) ) )
            L1_dists_ma = ma.masked_less_equal( L1_dists, epsilon, False )
            # distances use L2 norm
            raw_dist_mat = squareform( pdist( w_train_featspace, 'sqeuclidean' ) )
            dist_mat = ma.masked_array( raw_dist_mat, mask=L1_dists_ma.mask)
        else:
            from scipy.spatial.distance import cdist
            w_test_featspace = test_set.data_matrix * wts
            # collisions defined by L1 norm, aka taxicab, aka Manhattan, aka cityblock
            # in UNWEIGHTED FEATURE SPACE
            L1_dists = np.absolute( cdist( training_set.data_matrix, test_set.data_matrix, 'cityblock' ) )
            L1_dists_ma = ma.masked_less_equal( L1_dists, epsilon, False )
            # distances use L2 norm
            raw_dist_mat = cdist( w_train_featspace, w_test_featspace, 'sqeuclidean' )
            dist_mat = ma.masked_array( raw_dist_mat, mask=L1_dists_ma.mask )

        # Create marginal probabilities from distance matrix:
        similarity_mat = np.power( dist_mat, -5 )

        first_time_through = True
        for test_samp_index, test_samp_sims in enumerate( similarity_mat.T ):
            result = SingleSampleClassification()
            per_class_sims_list = [ test_samp_sims[ class_slice ] \
                    for class_slice in slice_list ]
            class_siml_means = []
            for class_sims in per_class_sims_list:
                try:
                    val = class_sims.compressed().mean()
                except FloatingPointError:
                    # mean of empty slice raises floating point error
                    val = np.nan
                class_siml_means.append( val )

            class_siml_means = np.array( class_siml_means )
            if not np.any( np.isnan( class_siml_means ) ):
                result.normalization_factor = class_siml_means.sum()
                result.marginal_probabilities = \
                    class_siml_means / result.normalization_factor
                result.predicted_label = \
                    training_set.class_names[ result.marginal_probabilities.argmax() ]

            result.sample_group_id = test_set._contiguous_sample_group_ids[ test_samp_index ]
            # FIXME: better to explicitly set s_seq_ids rather than trust the user to set it
            result.sample_sequence_id = test_set._contiguous_sample_group_ids[ test_samp_index ]
            result.source_filepath = test_set._contiguous_sample_names[ test_samp_index ]
            result.ground_truth_label = test_set._contiguous_ground_truth_labels[ test_samp_index ]
            result.ground_truth_value = test_set._contiguous_ground_truth_values[ test_samp_index ]
            result.split_number = split_number
            result.name = name
            split_result.individual_results.append( result )
            if not quiet:
                result.Print( line_item=True, include_col_header=first_time_through,
                        training_set_class_names=training_set.class_names )
            first_time_through = False

        # Predicted value via class coefficients, if applicable
        if training_set.interpolation_coefficients is not None:
            predicted_values = []
            ground_truth_values = []
            for result in split_result.individual_results:
                if result.marginal_probabilities is not None:
                    result.predicted_value = np.sum( result.marginal_probabilities * \
                                                 training_set.interpolation_coefficients )
                    predicted_values.append( result.predicted_value )
                    ground_truth_values.append( result.ground_truth_value )
            if predicted_values:
                split_result.predicted_values = predicted_values
                split_result.ground_truth_values = ground_truth_values

        # TILING SECTION:
        # Create a whole image classification result that
        # is the average of all the calls from all the tiles
        if test_set.num_samples_per_group > 1:
            split_result.averaged_results = []
            averaged_predicted_values = []
            averaged_ground_truth_values = []
            for start_index in xrange( 0, test_set.num_samples, test_set.num_samples_per_group ):
                end_index = start_index + test_set.num_samples_per_group
                tiles = split_result.individual_results[ start_index: end_index ]
                avg_result = AveragedSingleSamplePrediction( tiles, training_set.class_names )

                if avg_result.predicted_value is not None:
                    averaged_predicted_values.append( avg_result.predicted_value )
                    averaged_ground_truth_values.append( avg_result.ground_truth_value )

                split_result.averaged_results.append( avg_result )
                if not quiet:
                    avg_result.Print( line_item=True )
            if averaged_predicted_values:
                split_result.averaged_predicted_values = averaged_predicted_values
                split_result.averaged_ground_truth_values = averaged_ground_truth_values

        np.seterr (all='raise')
        return split_result

#=================================================================================
class FeatureSpaceRegression( _FeatureSpacePrediction ):
    """Container for SingleSampleRegression instances.
    Use for regressing all samples in a FeatureSpace in one sitting."""

    #: obj_count class attribute
    obj_count = 0

    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, split_number=None"""

        super( FeatureSpaceRegression, self ).__init__( *args, **kwargs )
        self.predicted_values = []

    #==============================================================    
    def __str__( self ):
        outstr = '<' + self.__class__.__name__
        if self.split_number is not None:
            outstr += ' #' + str( self.split_number )
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
    def NewMultivariateLinear( cls, test_set, feature_weights, name=None, split_number=None,
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

        split_result = cls( feature_weights.associated_feature_space, test_set, feature_weights,
                name, split_number )

        if test_set.ground_truth_values is not None and \
                len( test_set.ground_truth_values ) != 0:
            split_result.ground_truth_values = test_set.ground_truth_values

        for test_image_index in range( test_set.num_samples ):
            one_image_features = test_set.data_matrix[ test_image_index,: ]
            result = SingleSampleRegression._MultivariateLinear( one_image_features, feature_weights )
            result.split_number = split_result.split_number
            result.name = name
            result.source_filepath = test_set.sample_names[ test_image_index ]
            result.ground_truth_value = test_set.ground_truth_values[ test_image_index ]
            split_result.predicted_values.append( result.predicted_value )

            if not quiet:
                result.Print( line_item = True )
            split_result.individual_results.append( result )

        split_result.GenerateStats()
        return split_result

    #=====================================================================
    @classmethod
    def NewLeastSquares( cls, training_set, test_set, feature_weights, name=None,
            split_number=None, leave_one_out=False, quiet=False ):
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

        split_result = cls( training_set, test_set, feature_weights, name, split_number )

        if augmented_test_set.ground_truth_values is not None and \
                len( augmented_test_set.ground_truth_values ) != 0:
            split_result.ground_truth_values = augmented_test_set.ground_truth_values

        intermediate_train_set = augmented_train_set
        for test_image_index in range( augmented_test_set.num_samples ):
            if leave_one_out:
                sample_group_id = augmented_train_set.sample_group_ids[ test_image_index ]
                intermediate_train_set = augmented_train_set.SampleReduce(\
                        leave_out_sample_group_ids = sample_group_id, override=True)
            one_image_features = augmented_test_set.data_matrix[ test_image_index,: ]
            result = SingleSampleRegression._LeastSquares( \
                           one_image_features, intermediate_train_set )
            result.split_number = split_result.split_number
            result.name = name
            result.source_filepath = augmented_test_set.sample_names[ test_image_index ]
            result.ground_truth_value = augmented_test_set.ground_truth_values[ test_image_index ]
            split_result.predicted_values.append( result.predicted_value )

            if not quiet:
                result.Print( line_item = True )
            split_result.individual_results.append( result )

        # return settings to original
        np.seterr(**oldsettings)

        split_result.GenerateStats()
        return split_result
