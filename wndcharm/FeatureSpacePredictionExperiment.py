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
from .FeatureSpace import FeatureSpace, CheckIfClassNamesAreInterpolatable
from .FeatureSpacePrediction import FeatureSpacePrediction, FeatureSpaceClassification, \
        FeatureSpaceRegression
from .FeatureWeights import FisherFeatureWeights, PearsonFeatureWeights
from .SingleSamplePrediction import SingleSampleClassification

#============================================================================
class FeatureSpacePredictionExperiment( FeatureSpacePrediction ):
    """Base class container for FeatureSpacePrediction instances,
    i.e., when classifying/regressing a FeatureSpace multiple times, as in train/test splits.
    Methods here aggregate results and calculate statistics across splits."""

    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

        super( FeatureSpacePredictionExperiment, self ).__init__( *args, **kwargs )

        #: A dictionary where the name is the key, and the value is a list of individual results
        self.accumulated_individual_results = None

        self.feature_weight_statistics = None

        #: keep track of stats related to predicted values for reporting purposes
        self.individual_stats = None

    #=====================================================================
    def GenerateStats( self ):
        """Only partially implemented.
        
        Completed functionality:
        1. Simple aggregation of ground truth->predicted value pairs across splits.
           Use the function PerSampleStatistics() to average results for specific
             images across splits.
        2. Calculation of aggregated feature weight statistics

        Considerations for future implementation.
        1. The test set may or may not have ground truth (discrete and continuous)
        2. The results may not have a predicted value (discrete only)
        3. Continuous classifications do not have marginal probabilities
        4. Hybrid test sets (discrete test sets loaded into a continuous test set)
           have "pseudo-classes," i.e., easily binnable ground truth values."""

        # Simple result aggregation.
        # Calls GenerateStats() on the individual batches if the 
        # ground truth->predicted value pairs haven't been scraped
        # from the batch's list of individual ImageClassifications.
        from itertools import chain

        lists_of_ground_truths = []
        lists_of_predicted_values = []

        self.figure_of_merit = 0
        for batch_result in self.individual_results:
            self.num_classifications += len( batch_result.individual_results )
            if batch_result.figure_of_merit == None:
                batch_result.GenerateStats()
            if batch_result.ground_truth_values:
                lists_of_ground_truths.append( batch_result.ground_truth_values )
            if batch_result.predicted_values:
                lists_of_predicted_values.append( batch_result.predicted_values )

        if lists_of_ground_truths:
            self.ground_truth_values = list( chain( *lists_of_ground_truths ) )
        if lists_of_predicted_values:
            self.predicted_values = list( chain( *lists_of_predicted_values ) )

        # Aggregate feature weight statistics across splits, if any:
        feature_weight_lists = {}
        for batch_result in self.individual_results:
            if not batch_result.feature_weights:
                continue
            weight_names_and_values = zip( batch_result.feature_weights.featurenames_list, 
                                                            batch_result.feature_weights.values)
            for name, weight in weight_names_and_values:
                if not name in feature_weight_lists:
                    feature_weight_lists[ name ] = []
                feature_weight_lists[ name ].append( weight )

        if feature_weight_lists is not {}:
            for feature_name in feature_weight_lists:
                feature_weight_lists[ feature_name ] = np.array( feature_weight_lists[ feature_name ] )

            feature_weight_stats = []
            for fname in feature_weight_lists:
                fwl = feature_weight_lists[ fname ]
                count = len( fwl )
                fwl_w_zeros = np.zeros( len( self.individual_results ) )
                fwl_w_zeros[0:count] = fwl
                feature_weight_stats.append( ( np.mean( fwl_w_zeros ),
                                             count, np.std(fwl), np.min(fwl), np.max(fwl), fname ) )

            # Sort on mean values, i.e. index 0 of tuple created above
            self.feature_weight_statistics = sorted( feature_weight_stats, key=lambda a: a[0], reverse=True )

        # Remember, there's no such thing as a confusion matrix for a continuous class
        return self

    #=====================================================================
    @output_railroad_switch
    def WeightsOptimizationAnalysis( self ):
        """For use in only those classification experiments where predicted values are generated
        as part of the classification.
        
        This function is used when the user has varied the number/composition of image features
        in each of the batches and wants to figure out which classifier results in the best
        figure of merit. This function evaluates which of the batches was created using
        the optimal classifier configuration and outputs that."""

        print "Calculating optimized continuous classifier for training set {0}".\
                format( self.training_set.source_filepath )

        best_classification_result = None
        best_weights = None
        opt_number_features = None
        min_std_err = float( "inf" )

        for batch_result in self.individual_results:
            weights = batch_result.feature_weights
            num_features = len( weights.values )
            if batch_result.figure_of_merit < min_std_err:
                min_std_err = batch_result.figure_of_merit
                best_classification_result = batch_result
                best_weights = weights
                opt_number_features = num_features

        print "Optimum number of features: {0}".format( opt_number_features )
        best_classification_result.Print()
        best_weights.Print()
        print ""
        print "Aggregate feature weight analysis:"
        print "-----------------------------------"
        print "Legend:"
        print "NUM - Number of features used in aggregate / Individual feature rank"
        print "ASE - Standard Error of Final Predicted Value (using aggregated feature) vs ground truth"
        print "APC - Pearson correlation coefficient of Final Predicted Values vs ground truth"
        print "APE - Standard Error of APC"
        print "APP - P-value of APC"
        print "ASC - Spearman correlation coefficient of Final Predicted Values vs ground truth"
        print "APP - P-value of ASC"
        print ""

        print "NUM\tASE\tAPC\tAPE\tAPP\tASC\tAPP"
        print "===\t===\t===\t===\t===\t===\t==="
        for result in self.individual_results:
            line_item = "{0}\t".format( len( result.feature_weights.values ) ) # NUM
            line_item += "{0:.4f}\t".format( result.figure_of_merit ) # ASE
            line_item += "{0:.4f}\t".format( result.pearson_coeff ) # APC
            line_item += "{0:.4f}\t".format( result.pearson_std_err ) # APE
            line_item += "{0:.4f}\t".format( result.pearson_p_value ) # APP
            line_item += "{0:.4f}\t".format( result.spearman_coeff ) # ASC
            line_item += "{0:.4f}\t".format( result.spearman_p_value ) # ASP
            print line_item

    #=====================================================================
    @classmethod
    def NewShuffleSplit( cls, feature_set, n_iter=5, name=None, features_size=0.15,
                           train_size=None, test_size=None, random_state=True, classifier=None,
                           quiet=False):
        """args train_size, test_size, and random_state are all passed through to Split()
        feature_size if a float is feature usage fraction, if in is top n features."""

        experiment = cls( training_set=feature_set, test_set=feature_set, name=name )
        if isinstance( features_size, float ):
            if features_size < 0 or features_size > 1.0:
                raise ValueError('Arg "features_size" must be on interval [0,1] if a float.')
            num_features = int( round( features_size * feature_set.num_features ) )
        elif isinstance( features_size, int ):
            if features_size < 0 or features_size > feature_set.num_features:
                raise ValueError( 'must specify num_features or feature_usage_fraction in kwargs')
            num_features = features_size
        else:
            raise ValueError( 'Arg "features_size" must be valid float or int.' )

        if not quiet:
            print "using top "+str (num_features)+" features"

        # If you passed the same random_state into Split, you'd get the same exact split for
        # all n_iter. Therefore use the seed passed in here to predictably generate seeds
        # for the Split() iterations.
        if random_state:
            from numpy.random import RandomState
            from functools import partial
            maxint = 2 ** 32 - 1
            if random_state is True:
                from numpy.random import randint as np_randint
                randint = partial( np_randint, low=0, high=maxint )
            elif type( random_state ) is int:
                randint = partial( RandomState( random_state ).randint, low=0, high=maxint )
            elif type( random_state ) is RandomState:
                randint = partial( random_state.randint, low=0, high=maxint )
            else:
                raise ValueError( "arg random_state must be an int, instance of numpy.random.RandomState, or True")
        else:
            randint = lambda: None # samples split the same way all iterations - not recommended!

        for split_index in range( n_iter ):
            train_set, test_set = feature_set.Split(
                train_size, test_size, random_state=randint(), quiet=quiet )
            train_set.Normalize( quiet=quiet )
            
            if feature_set.discrete:
                weights = \
                  FisherFeatureWeights.NewFromFeatureSpace( train_set ).Threshold( num_features )
            else:    
                weights = \
                  PearsonFeatureWeights.NewFromFeatureSpace( train_set ).Threshold( num_features )

            if not quiet:
                weights.Print()
            reduced_train_set = train_set.FeatureReduce( weights )
            reduced_test_set = test_set.FeatureReduce( weights )
            reduced_test_set.Normalize( reduced_train_set, quiet=quiet )

            if feature_set.discrete:
                batch_result = FeatureSpaceClassification.NewWND5( reduced_train_set, \
                 reduced_test_set, weights, batch_number=split_index, quiet=quiet )
            else:
                if classifier == 'linear':
                    batch_result = FeatureSpaceRegression.NewMultivariateLinear(
                            reduced_train_set, weights, batch_number=split_index, quiet=quiet )
                else: # default classifier == 'lstsq':
                    batch_result = FeatureSpaceRegression.NewLeastSquares(
                        reduced_train_set, reduced_test_set, weights, batch_number=split_index, quiet=quiet )

            experiment.individual_results.append( batch_result )

        return experiment

#============================================================================
class FeatureSpaceClassificationExperiment( FeatureSpacePredictionExperiment ):
    """Container for FeatureSpaceClassifications instances,
    i.e., when classifying a FeatureSpace multiple times, as in train/test splits.
    Methods here aggregate results and calculate statistics across splits.

    The information contained here comprises everything that would appear in an
    HTML file generated by the C++ implementation of WND-CHARM.

    Member "self.figure_of_merit" represents the mean classification accuracy across
    all FeatureSpaceClassification instances in "self.individual_results"."""

    def __init__( self, *args, **kwargs ):
        """Possible kwargs, with defaults:
        training_set=None, test_set=None, feature_weights=None, name=None, batch_number=None"""

        super( FeatureSpaceClassificationExperiment, self ).__init__( *args, **kwargs )

        self.num_correct_classifications = None

        self.confusion_matrix = None
        self.average_similarity_matrix = None
        self.average_class_probability_matrix = None

    #=====================================================================
    def GenerateStats( self ):
        """Not fully implemented yet. Need to implement generation of confusion, similarity, and
        average class probability matrices from constituent batches."""

        # Base class does feature weight analysis, ground truth-pred. value aggregation
        super( FeatureSpaceClassificationExperiment, self ).GenerateStats()
        
        self.num_classifications = 0
        self.num_correct_classifications = 0
        for batch_result in self.individual_results:
            for indiv_result in batch_result.individual_results:
                self.num_classifications += 1
                if indiv_result.ground_truth_class_name == indiv_result.predicted_class_name:
                    self.num_correct_classifications += 1

        self.figure_of_merit = float( self.num_correct_classifications) / float( self.num_classifications )
        self.classification_accuracy = self.figure_of_merit

        #FIXME: Create confusion, similarity, and class probability matrices
        return self

    #=====================================================================
    @output_railroad_switch
    def Print( self ):
        """Generate and output statistics across all batches, as well as the figures of merit
        for each individual batch."""
        
        if self.figure_of_merit == None:
            self.GenerateStats()

        print "==========================================="
        print "Experiment name: {0}".format( self.name )
        print "Summary:"
        print "Number of batches: {0}".format( len( self.individual_results ) )
        print "Total number of classifications: {0}".format( self.num_classifications )
        print "Total number of CORRECT classifications: {0}".format( self.num_correct_classifications )
        print "Total classification accuracy: {0:0.4f}\n".format( self.classification_accuracy )

        print "Batch Accuracies:"
        print "#\tAcc.\tName"
        print "------------------------------------"

        for batch_result in sorted( self.individual_results, key=lambda X: X.classification_accuracy, reverse=True ):
            if batch_result.figure_of_merit == None:
                batch_result.GenerateStats()
            print "{0}\t{1:0.4f}\t{2}".format( batch_result.batch_number,
                              batch_result.classification_accuracy, batch_result.name )

        outstr = "{0}\t{1:0.3f}\t{2:>3}\t{3:0.3f}\t{4:0.3f}\t{5:0.3f}\t{6}"
        print "Feature Weight Analysis (top 50 features):"
        print "Rank\tmean\tcount\tStdDev\tMin\tMax\tName"
        print "----\t----\t-----\t------\t---\t---\t----"
        for count, fw_stat in enumerate( self.feature_weight_statistics, 1 ):
            print outstr.format( count, *fw_stat )
            if count >= 50:
                    break

    #=====================================================================
    @classmethod
    @output_railroad_switch
    def NewFromHTMLReport( cls, path_to_html ):
        """Takes as input an HTML report generated by C++ implementation wndchrm
        Parses the report and builds up a Pychrm representation of the results,
        which facilitates analysis, graphing, etc."""

        import re
        row_re = re.compile( r'<tr>(.+?)</tr>' )
        name_re = re.compile( r'"(.+?)"' )
        num_re = re.compile( r'(\d*\.?\d+)' )

        # FIXME: This should fail if there isn't some part of the class names that are interpretable
        # as a number, specifically when it tries to calculate an "interpolated" (predicted) value
        # for the sample based on marginal probabilities.

        def ParseClassSummaryHTML( the_html ):
            rows = row_re.findall( the_html )
            ts = FeatureSpace()
            ts.num_classes = 0
            ts.interpolation_coefficients = []
            ts.classnames_list = []
            for rownum, row in enumerate( rows ):
                if rownum == 0:
                    continue # skip column header
                ts.num_classes += 1
                classname = re.search( r'<th>(.+?)</th>', row ).group(1)
                ts.classnames_list.append( classname )
            ts.interpolation_coefficients = CheckIfClassNamesAreInterpolatable( ts.classnames_list )
            return ts

        # The following will be determined once the number of classes has been ascertained
        normalization_col = 1
        mp_col = 2
        ground_truth_col = None
        predicted_col = None
        interp_val_col = None # We don't even use this
        name_col = None

        _training_set = None
        _test_set = None
        exp = cls()
        exp.name = path_to_html

        trainingset_definition = False
        trainingset_html = ""
        testset_definition = False
        testset_html = ""

        insidesplit = False
        split = None
        splitcount = 0
        split_linecount = None
        with open( path_to_html ) as file:
            for line in file:
                if 'trainset_summary' in line:
                    trainingset_definition = True
                elif trainingset_definition == True:
                    trainingset_html += line.strip()
                    if '</table>' in line:
                        trainingset_definition = False
                        ts = _training_set = ParseClassSummaryHTML( trainingset_html )
                        ground_truth_col = ts.num_classes + 3
                        predicted_col = ts.num_classes + 4
                        interp_val_col = ts.num_classes + 6
                        name_col = ts.num_classes + 7

                elif 'testset_summary' in line:
                    testset_definition = True
                elif testset_definition == True:
                    testset_html += line.strip()
                    if '</table>' in line:
                        testset_definition = False
                        _test_set = ParseClassSummaryHTML( testset_html )

                elif line.startswith( '<TABLE ID="IndividualImages_split' ):
                    # If we haven't seen a test set definition by now, we ain't gonna see one period.
                    if not _test_set:
                        _test_set = _training_set
                    insidesplit = True
                    split = FeatureSpaceClassification( training_set=_training_set, test_set=_test_set )
                    split.predicted_values = []
                    split.ground_truth_values = []
                    splitcount += 1
                    split_linecount = 0
                elif insidesplit and line.startswith( '</table><br><br>' ):
                    insidesplit = False
                    exp.individual_results.append( split )
                elif insidesplit:
                    split_linecount += 1
                    if split_linecount == 1:
                        # First line in individual results is column headers
                        # Pure clasification without interpolation
                        if 'Interpolated Value' not in line:
                            interp_val_col = None
                            name_col = ts.num_classes + 6 # one less than set above -- that column won't exist
                        continue
                    noends = line.strip( '<trd/>\n' ) # take the tr and td tags off front end
                    values = noends.split( '</td><td>' )
                    result = SingleSampleClassification()

                    result.normalization_factor = float( values[ normalization_col ] )
                    result.marginal_probabilities = \
                            [ float( val.strip( '</b>' ) ) for val in values[ mp_col : mp_col + _training_set.num_classes ] ]
                    result.predicted_class_name = values[ predicted_col ]
                    # Sometimes c-chrm labels classes with a * to say it's not part of the training set
                    result.ground_truth_class_name = values[ ground_truth_col ].strip('*')
                    result.name = name_re.search( values[ name_col ] ).groups()[0]
                    result.source_filepath = result.name
                    if ts.interpolation_coefficients is not None:
                        result.ground_truth_value = \
                        ts.interpolation_coefficients[ ts.classnames_list.index(result.ground_truth_class_name ) ]
                        #result.predicted_value = float( values[ interp_val_col ] )
                        result.predicted_value = \
                        sum( [ x*y for x,y in zip( result.marginal_probabilities, _training_set.interpolation_coefficients ) ] )
                        split.predicted_values.append( result.predicted_value )
                        split.ground_truth_values.append( result.ground_truth_value )
                    #result.Print( line_item = True )
                    result.batch_number = splitcount
                    split.individual_results.append(result)

        exp.training_set = _training_set
        exp.test_set = _test_set

        exp.GenerateStats()
        return exp

    #=====================================================================
    @output_railroad_switch
    def PerSampleStatistics( self ):
        """Characterizes variability of regressed predicted values across batches.
        SingleSamplePrediction info is aggregated for each individual sample,
        statistics calculated and printed out."""

        if self.individual_results == 0:
            raise ValueError( 'No batch results to analyze' )

        #self.predicted_values = []
        self.ground_truth_values = []

        self.accumulated_individual_results = {}
        self.individual_stats = {}

        for batch in self.individual_results:
            for result in batch.individual_results:
                if not result.source_filepath in self.accumulated_individual_results:
                    # initialize list of individual results for this file
                    self.accumulated_individual_results[ result.source_filepath ] = []
                self.accumulated_individual_results[ result.source_filepath ].append( result )

        for filename in self.accumulated_individual_results:

            # Get marginal probability averages
            mp_totals = None
            for result in self.accumulated_individual_results[filename]:
                if not mp_totals:
                    mp_totals = result.marginal_probabilities[:]
                else:
                    new_total = []
                    for class_total, new_mp in zip( mp_totals, result.marginal_probabilities ):
                        new_total.append( class_total + new_mp )
                    mp_totals = new_total

            mp_avgs = [ float(mp_totals[i]) / len( self.accumulated_individual_results[filename] ) for i in range( len( mp_totals ) ) ]
            #vals = np.array ([result.predicted_value for result in self.accumulated_individual_results[filename] ])
            vals = [result.predicted_class_name for result in self.accumulated_individual_results[filename] ]
            #self.ground_truth_values.append( self.accumulated_individual_results[filename][0].ground_truth_value )
            gt_class = self.accumulated_individual_results[filename][0].ground_truth_class_name
            self.ground_truth_values.append( gt_class )
            #self.predicted_values.append( np.mean(vals) )
            self.individual_stats[filename] = ( len(vals), float( vals.count( gt_class ) ) / len(vals), mp_avgs, gt_class )

        print "==========================================="
        print 'Experiment name: "{0}"'.format( self.name ) + ' Individual results\n'

        mp_delim = "  "
        discrlineoutstr = "\tsplit {split_num:02d}: pred: {pred_class}\tact: {actual_class}\tnorm factor: {norm_factor:0.3g},\tmarg probs: ( {norm_dists} )"
        outstr = "\t---> Tested {0} times, avg correct: {1:0.3f}, avg marg probs ( {2} )"

        #create view
        res_dict = self.accumulated_individual_results

        # sort by ground truth, then alphanum
        sort_func = lambda A, B: cmp( A, B ) if res_dict[A][0].ground_truth_class_name == res_dict[B][0].ground_truth_class_name else cmp( res_dict[A][0].source_filepath, res_dict[B][0].source_filepath  ) 
        sorted_images = sorted( self.accumulated_individual_results.iterkeys(), sort_func )

        for samplename in sorted_images:
            print 'File "' + samplename + '"'
            for result in self.accumulated_individual_results[ samplename ]:
                marg_probs = [ "{0:0.3f}".format( num ) for num in result.marginal_probabilities ]
                print discrlineoutstr.format( split_num = result.batch_number, \
                                         pred_class = result.predicted_class_name, \
                                         actual_class = result.ground_truth_class_name, \
                                         norm_factor = result.normalization_factor, \
                                         norm_dists = mp_delim.join( marg_probs ) )

            marg_probs = [ "{0:0.3f}".format( num ) for num in self.individual_stats[ samplename ][2] ]
            print outstr.format( self.individual_stats[ samplename ][0], self.individual_stats[ samplename ][1], mp_delim.join( marg_probs ) )

        # If 2 or 3 class problem, plot individuals in marginal probability space
# END class definition for FeatureSpaceClassificationExperiment

#============================================================================
class FeatureSpaceRegressionExperiment( FeatureSpacePredictionExperiment ):
    """Container for FeatureSpaceRegression instances,
    i.e., when regressing a FeatureSpace multiple times, as in train/test splits.
    Methods here aggregate results and calculate statistics across splits.

    Member "self.figure_of_merit" represents the mean Standard Error across
    FeatureSpaceRegression instances in "self.individual_results"."""

    def __init__( self, *args, **kwargs):
        super( FeatureSpaceRegressionExperiment, self ).__init__( *args, **kwargs )

    #=====================================================================
    def GenerateStats( self ):
        """Calculates statistics describing how well predicted values
        correlate with ground truth across all batches.

        Requires scipy.stats package to be installed"""

        # Base class does feature weight analysis, ground truth-pred. value aggregation
        super( FeatureSpaceRegressionExperiment, self ).GenerateStats()
    
        gt = np.array( self.ground_truth_values )
        pv = np.array( self.predicted_values )
        
        diffs = gt - pv
        diffs = np.square( diffs )
        err_sum = np.sum( diffs )

        from math import sqrt
        from scipy.stats import linregress, spearmanr

        self.figure_of_merit = sqrt( err_sum / self.num_classifications )

        # For now, ignore "FloatingPointError: 'underflow encountered in stdtr'"
        np.seterr (under='ignore')
        slope, intercept, self.pearson_coeff, self.pearson_p_value, self.pearson_std_err = \
                                 linregress( self.ground_truth_values, self.predicted_values )

        try:
            self.spearman_coeff, self.spearman_p_value =\
                     spearmanr( self.ground_truth_values, self.predicted_values )
        except FloatingPointError:
            self.spearman_coeff, self.spearman_p_value = ( 0, 1 )

        np.seterr (all='raise')
        return self

    #=====================================================================
    @output_railroad_switch
    def Print( self ):
        """Output statistics from this experiment."""
        if self.figure_of_merit == None:
            self.GenerateStats()

        print "==========================================="
        print "Experiment name: {0}".format( self.name )
        print "Summary:"
        print "Number of batches: {0}".format( len( self.individual_results ) )
        print "Total number of classifications: {0}".format( self.num_classifications )
        print "Total standard error: {0}".format( self.figure_of_merit )
        print "Pearson corellation coefficient: {0}".format( self.pearson_coeff )
        print "Pearson p-value: {0}".format( self.pearson_p_value )        

        outstr = "{0}\t{1:0.3f}\t{2:>3}\t{3:0.3f}\t{4:0.3f}\t{5:0.3f}\t{6}"
        print "Feature Weight Analysis (top 50 features):"
        print "Rank\tmean\tcount\tStdDev\tMin\tMax\tName"
        print "----\t----\t-----\t------\t---\t---\t----"
        for count, fw_stat in enumerate( self.feature_weight_statistics, 1 ):
            print outstr.format( count, *fw_stat )
            if count >= 50:
                    break

    #=====================================================================
    @output_railroad_switch
    def PerSampleStatistics( self ):
        """Characterizes variability of classifications across batches.
        SingleSamplePrediction info is aggregated for each individual sample,
        statistics calculated and printed out."""

        if self.individual_results == 0:
            raise ValueError( 'No batch results to analyze' )

        self.predicted_values = []
        self.ground_truth_values = []

        self.accumulated_individual_results = {}
        self.individual_stats = {}

        for batch in self.individual_results:
            for result in batch.individual_results:
                if not result.source_filepath in self.accumulated_individual_results:
                    # initialize list of individual results for this file
                    self.accumulated_individual_results[ result.source_filepath ] = []
                self.accumulated_individual_results[ result.source_filepath ].append( result )

        for filename in self.accumulated_individual_results:
            vals = np.array( [result.predicted_value for result in self.accumulated_individual_results[filename] ])
            self.ground_truth_values.append( self.accumulated_individual_results[filename][0].ground_truth_value )
            self.predicted_values.append( np.mean(vals) )
            self.individual_stats[filename] = ( len(vals), np.min(vals), np.mean(vals), \
                                                                                            np.max(vals), np.std(vals) ) 

        print "==========================================="
        print 'Experiment name: "{0}"'.format( self.name ) + ' Individual results\n'
        mp = "  "
        contlineoutstr = "\tsplit {split_num:02d} '{batch_name}': actual: {actual_class}. Predicted val: {pred_val:0.3f}"
        outstr = "\t---> Tested {0} times, low {1:0.3f}, mean {2:0.3f}, high {3:0.3f}, std dev {4:0.3f}"

        #create view
        res_dict = self.accumulated_individual_results

        # sort by ground truth, then alphanum
        sort_func = lambda A, B: cmp( A, B ) if res_dict[A][0].ground_truth_value == res_dict[B][0].ground_truth_value else cmp( res_dict[A][0].ground_truth_value, res_dict[B][0].ground_truth_value  ) 
        sorted_images = sorted( self.accumulated_individual_results.iterkeys(), sort_func )

        for samplename in sorted_images:
            print 'File "' + samplename + '"'
            for result in self.accumulated_individual_results[ samplename ]:
                print contlineoutstr.format( split_num = result.batch_number, \
                                         batch_name = result.name, \
                                         actual_class = result.ground_truth_value, \
                                         pred_val = result.predicted_value )
            print outstr.format( *self.individual_stats[ samplename ] )

