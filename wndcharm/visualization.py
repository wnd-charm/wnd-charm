#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (C) 2015 National Institutes of Health
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Written by:  Christopher Coletta (github.com/colettace)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from .utils import output_railroad_switch
from .FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment

#============================================================================
class _BaseGraph( object ):
    """An abstract base class that is supposed to hold onto objects on which to call
    matplotlib.pyplot API methods."""

    def __init__( self ):

        # general stuff:
        self.chart_title = None
        self.file_name = None
        self.split_result = None

        # pyplot-specific stuff
        self.figure = None
        self.main_axes = None

    def SaveToFile( self, filepath ):
    
        if self.figure == None:
            raise ValueError( 'No figure to save!' )
        self.figure.savefig( filepath )
        print 'Wrote chart "{0}" to file "{1}"'.format( self.chart_title, filepath )
            
#============================================================================
class PredictedValuesGraph( _BaseGraph ):
    """This is a concrete class that can produce two types of graphs that are produced
    from SingleSamplePrediction data stored in a FeatureSpacePrediction."""

    #=================================================================
    def __init__( self, result, name=None, use_averaged_results=True, figsize=(15, 10), fontsize=16, titlefontsize=18 ):
        """Constructor sorts ground truth values contained in FeatureSpacePrediction
        and loads them into self.grouped_coords
        
        use_averaged_results - bool - If this object has averaged results (due to tiling or 
            "per sample" aggregation across splits, use those results instead of 
            individual results."""

        #FIXME: implement user-definable bin edges

        self.split_result = result
        if name is None:
            name = result.name

        self.chart_title = name
        self.figsize = figsize
        self.fontsize = fontsize
        self.titlefontsize = titlefontsize

        gt_vals, pred_vals = result.RankOrderSort( use_averaged_results=use_averaged_results )
        whole_list = zip( gt_vals, pred_vals )

        self.grouped_coords = {}

        if result.test_set.discrete:
            self.class_names = result.test_set.class_names
            self.class_values = result.test_set.interpolation_coefficients
            self.num_classes = result.test_set.num_classes
            for class_val, class_name in zip( self.class_values, self.class_names ):
                self.grouped_coords[ class_name ] = \
                        [ xy for xy in whole_list if xy[0] == class_val ]
        else:
            class_name = result.test_set.name
            self.class_names = [ class_name ]
            self.class_values = [ 1 ]
            self.num_classes = 1
            self.grouped_coords[ class_name ] = whole_list

        _min = min( self.class_values )
        ampl = max( self.class_values ) - _min
     
        import matplotlib.pyplot as plt
        self.class_colors = plt.cm.jet( [ float(val -_min)/ampl for val in self.class_values ] )

    #=================================================================
    def show( self ):
        self.figure.show()

    #=====================================================================
    @classmethod
    @output_railroad_switch
    def NewFromHTMLReport( cls, filepath, **kwargs):
        """Helper function to facilitate the fast generation of graphs from C++-generated
        HTML Report files."""

        exp = FeatureSpaceClassificationExperiment.NewFromHTMLReport( filepath )
        exp.GenerateStats()
        exp.PerSampleStatistics( quiet=True )
        newgraphobj = cls( exp, **kwargs )
        return newgraphobj

    #=====================================================================
    def RankOrderedPredictedValuesGraph( self, chart_title=None, fig=None, ax=None ):
        """Visualizes the distribution of predicted values generated by classification.

        For each individual ImageClassification with ground truth value (i.e., class id) and
        predicted value, all results are grouped within their class, sorted by predicted value
        in ascending order, then ploted side-by-side.

        Requires the package matplotlib to be installed.

        fig - matplotlib.figure.Figure, or None (default):
            If provided, generate axes on the passed-in figure obj, otherwise
            instantiate a new one.

        ax - matplotlib.axes._subplots.AxesSubplot or None (default):
            Plot data on these axes, if provided, otherwise create new axes on figure."""

        import matplotlib.pyplot as plt

        if fig is None and ax is None:
            fig = plt.figure( figsize=self.figsize )
            ax = fig.add_subplot(111)
        elif fig is None:
            fig = ax.get_figure()
        elif ax is None:
            ax = fig.gca()

        self.figure = fig
        self.main_axes = ax

        if chart_title:
            self.chart_title = chart_title

        self.main_axes.set_title( self.chart_title, size=self.titlefontsize )
        self.main_axes.set_xlabel( 'count', size=self.fontsize )
        self.main_axes.set_ylabel( 'Predicted Value Scores', size=self.fontsize )

        abscissa_index = 1

        for class_name, class_color in zip( self.class_names, self.class_colors ):
            ground_truth_vals, predicted_vals = zip( *self.grouped_coords[ class_name ] )
            x_vals = [ i + abscissa_index for i in range( len( ground_truth_vals ) ) ]
            self.main_axes.scatter( x_vals, predicted_vals, c=class_color, marker='o',
                    s=80, edgecolor='none', alpha=0.3, label=class_name )
            abscissa_index += len( ground_truth_vals )

        #self.main_axes.legend( loc = 'lower right')
        self.main_axes.legend( loc = 'lower right', fontsize=self.fontsize )
        return self
        
    #=====================================================================
    def KernelSmoothedDensityGraph( self, chart_title=None, bw_method=None,
            fig=None, ax=None  ):
        """This graph visualizes the distribution of predicted values generated by classification.
        A kernel-smoothed probability density function is plotted for each image class on
        the same chart allowing comparison of distribution of predicted values amoung image class.

        Requires the packages matplotlib and scipy. Uses scipy.stats.gaussian_kde to
        generate kernel-smoothed probability density functions.

        bw_method -
            Arg passed directly to bw_method arg in scipy.stats.gaussian_kde

        fig - matplotlib.figure.Figure, or None (default):
            If provided, generate axes on the passed-in figure obj, otherwise
            instantiate a new one.

        ax - matplotlib.axes._subplots.AxesSubplot or None (default):
            Plot data on these axes, if provided, otherwise create new axes on figure."""

        import matplotlib.pyplot as plt

        if fig is None and ax is None:
            fig = plt.figure( figsize=self.figsize )
            ax = fig.add_subplot(111)
        elif fig is None:
            fig = ax.get_figure()
        elif ax is None:
            ax = fig.gca()

        self.figure = fig
        self.main_axes = ax

        if chart_title:
            self.chart_title = chart_title

        self.main_axes.set_title( self.chart_title, size=self.titlefontsize )
        self.main_axes.set_xlabel( 'Age score', size=self.fontsize )
        self.main_axes.set_ylabel( 'Probability density', size=self.fontsize )

        from scipy.stats import gaussian_kde

        # For now, ignore "FloatingPointError: underflow encountered in exp'"
        np.seterr( under='ignore' )

        for class_name, class_color in zip( self.class_names, self.class_colors ):
            ground_truth_vals, predicted_vals = zip( *self.grouped_coords[ class_name ] )

            pred_vals = np.array( predicted_vals )
            lobound = pred_vals.min()
            hibound = pred_vals.max()
            kernel_smoother = gaussian_kde( pred_vals, bw_method=bw_method )
            intervals = np.mgrid[ lobound:hibound:100j ]
            density_estimates = kernel_smoother.evaluate( intervals )
            self.main_axes.plot( intervals, density_estimates, c=class_color,
                linewidth=4, label=class_name )

        np.seterr( all='raise' )

        self.main_axes.legend( fontsize=self.fontsize )
        return self

#============================================================================
class FeatureTimingVersusAccuracyGraph( _BaseGraph ):
    """A cost/benefit analysis of the number of features used and the time it takes to calculate
    that number of features for a single image"""

    #FIXME: Add ability to do the first 50 or 100 features, make the graph, then
    #       ability to resume from where it left off to do the next 50.

    def __init__( self, training_set, feature_weights, test_image_path,
        chart_title=None, max_num_features=300 ):

        self.timing_axes = None
        import time
        timings = []

        from .FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment
        from .SingleSamplePrediction import SingleSampleClassification
        from .FeatureSpacePrediction import FeatureSpaceClassification
        from .FeatureVector import FeatureVector

        experiment = FeatureSpaceClassificationExperiment( training_set, training_set, feature_weights )
        for number_of_features_to_use in range( 1, max_num_features + 1 ):

            reduced_ts = None
            reduced_fw = None
            three_timings = []
            # Take the best of 3
            for timing in range( 3 ):
                # Time the creation and classification of a single signature
                t1 = time.time()
                reduced_fw = feature_weights.Threshold( number_of_features_to_use )
                sig = FeatureVector( source_filepath=test_image_path, feature_names=reduced_fw.feature_names ).GenerateFeatures()
                reduced_ts = training_set.FeatureReduce( reduced_fw )
                sig.Normalize( reduced_ts )
        
                result = SingleSampleClassification.NewWND5( reduced_ts, reduced_fw, sig )
                result.Print()
                # FIXME: save intermediates just in case of interruption or parallization
                # result.PickleMe()
                t2 = time.time()
                three_timings.append( t2 - t1 )

            timings.append( min( three_timings ) )

            # now, do a fit-on-fit test to measure classification accuracy
            split_result = FeatureSpaceClassification.NewWND5( reduced_ts, reduced_ts, reduced_fw )
            split_result.Print()
            experiment.individual_results.append( split_result )

        import matplotlib.pyplot as plt

        x_vals = list( range( 1, max_num_features + 1 ) )

        self.figure = plt.figure()
        self.main_axes = self.figure.add_subplot(111)
        if chart_title == None:
            self.chart_title = "Feature timing v. classification accuracy"    
        else:
            self.chart_title = chart_title
        self.main_axes.set_title( self.chart_title )
        self.main_axes.set_xlabel( 'Number of features' )
        self.main_axes.set_ylabel( 'Classification accuracy (%)', color='b' )
        classification_accuracies = \
          [ split_result.classification_accuracy * 100 for split_result in experiment.individual_results ]

        self.main_axes.plot( x_vals, classification_accuracies, color='b', linewidth=2 )
        for tl in self.main_axes.get_yticklabels():
            tl.set_color('b')    

        self.timing_axes = self.main_axes.twinx()
        self.timing_axes.set_ylabel( 'Time to calculate features (s)', color='r' )
        self.timing_axes.plot( x_vals, timings, color='r' )
        for tl in self.timing_axes.get_yticklabels():
            tl.set_color('r')    

#============================================================================
class AccuracyVersusNumFeaturesGraph( _BaseGraph ):
    """Hyper-parameter optimization of number of features.
    Graphs the figure of merit as a function of number of top-ranked features."""

    def __init__( self, param_space=None, param_scale='log', lda_comparison=True,
            chart_title=None, figsize=(12, 8), y_min=None, y_max=None, quiet=True, **kwargs ):
        """Creates graph for Classifier/Regressor figure of merit as a function of
        number of top-ranked features used in classification.

        Calls wndcharm.FeatureSpacePredictionExperiment.NumFeaturesGridSearch, which is
        itself a wrapper for NewShuffleSplit, to which kwargs gets passed through.

        Args:
            param_space - iterable or int or None
                iterable of ints specifying num features to be used each iteration
                int - uses n intervals evenly spaced numbers along log scale from
                    1 to kwargs['feature_space'].num_features
                None - same as int above but specifying 20 intervals, results in param_space
                    1, 2, 3, 4, 7, 10, 16, 24, 36, 54, 80, 119, 178, 266, ... 2919
            param_scale - {'log', 'linear' }, pass through to matplotlib.ax.set_xscale()
            lda_comparison - bool - Do two runs with and without Linear Discriminant Analysis
                and graph the results on the same axes.
            y_min, y_max - float - chart param.
            **kwargs - passed directly through to NewShuffleSplit()
        """

        if kwargs['feature_space'].discrete:
            from .FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment as Experiment
        else:
            from .FeatureSpacePredictionExperiment import FeatureSpaceRegressionExperiment as Experiment

        if lda_comparison and 'lda' in kwargs:
            del kwargs['lda']

        if lda_comparison:
            main_coords = Experiment.NumFeaturesGridSearch( param_space=param_space,
                    quiet=quiet, lda=False, **kwargs )
            X, Y = zip( *main_coords )
            lda_coords = Experiment.NumFeaturesGridSearch( param_space=param_space,
                    quiet=quiet, lda=True, **kwargs )
            X_lda, Y_lda = zip( *lda_coords )
            all_ys = Y + Y_lda
            if y_min == None:
                y_min = min( all_ys )
            y_max = max( all_ys )
        else:
            main_coords = Experiment.NumFeaturesGridSearch( param_space=param_space,
                    quiet=quiet, **kwargs )
            X, Y = zip( *main_coords )
            if y_min == None:
                y_min = min( Y )
            y_max = max( Y )

        horiz_buffer = 0.05 * ( y_max - y_min )
        y_min -= horiz_buffer
        y_max += horiz_buffer
        total_n_feats = kwargs['feature_space'].num_features

        import matplotlib.pyplot as plt

        self.figure = plt.figure( figsize=figsize, facecolor='white' )
        self.main_axes = self.figure.add_subplot(111)
        if chart_title == None:
            self.chart_title = "Feature Space Predicton figure of merit vs. # features"
        else:
            self.chart_title = chart_title
        self.main_axes.set_title( self.chart_title, size=18 )

        self.main_axes.set_xlabel( '# top-ranked features', size=16 )
        self.main_axes.set_xscale( param_scale )
        self.main_axes.set_ylabel( 'Figure of Merit', color='b', size=16 )
        self.main_axes.set_ylim( [ y_min, y_max ] )
        self.main_axes.plot( X, Y, color='b', marker='o', linestyle='--' )
        for x, y in zip( X, Y ):
            text = '{:0.03} (n={}, frac={:0.3f})'.format( y, x, float(x)/total_n_feats )
            self.main_axes.annotate( text, xy=(x,y), rotation=-30 )

        for tl in self.main_axes.get_yticklabels():
            tl.set_color('b')

#        self.main_axes.annotate( 'min R={0:.3f} @ {1}'.format(min_ls_yval, optimal_num_feats_ls),
#        color='b',
#        xy=( optimal_num_feats_ls, min_ls_yval ),
#        xytext=( optimal_num_feats_ls, 0.8 * _max ),
#        arrowprops=dict(facecolor='black', shrink=0.05),
#        horizontalalignment='right' )

        if lda_comparison:
            self.lda_axes = self.main_axes.twinx()
            self.lda_axes.set_xscale( param_scale )
            self.lda_axes.set_ylabel( 'Figure of Merit WITH LDA', color='r', size=16 )
            self.lda_axes.set_ylim( [ y_min, y_max ] )
            self.lda_axes.plot( X_lda, Y_lda, color='r', marker='o', linestyle='--' )
            for tl in self.lda_axes.get_yticklabels():
                tl.set_color('r')

#        self.lda_axes.annotate( 'min R={0:.3f} @ {1}'.format(min_voting_yval, optimal_num_feats_voting),
#        color='r',
#        xy=( optimal_num_feats_voting, min_voting_yval ),
#        xytext=( optimal_num_feats_voting, 0.6 * _max ),
#        arrowprops=dict(facecolor='black', shrink=0.05),
#        horizontalalignment='right' )
