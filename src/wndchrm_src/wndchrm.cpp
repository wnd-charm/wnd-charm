/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2007 Open Microscopy Environment                             */
/*         Massachusetts Institue of Technology,                                 */
/*         National Institutes of Health,                                        */
/*         University of Dundee                                                  */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribusplit_numte it and/or     */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:  Lior Shamir <shamirl [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
// variadic
#include <stdarg.h>
// system errors
// isdigit
#include <ctype.h>
#include <algorithm>

#include "TrainingSet.h"
#include "wndchrm_error.h"
#include "Tasks.h"

#define MAX_SPLITS 10000
#define MAX_SAMPLES 190000

/* global variable */
extern int verbosity;

#include <sys/time.h>
void randomize() {
	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
}

void setup_featureset (featureset_t *featureset) {
	int	rotations,tiles_x,tiles_y;
	char sample_name[SAMPLE_NAME_LENGTH],preproc_name[SAMPLE_NAME_LENGTH];
	int sample_name_lngth,preproc_name_lngth;
	int rot_index,tile_index_x,tile_index_y;
	int n_samples = 0;
	preproc_opts_t *preproc_opts = &(featureset->preproc_opts);
	sampling_opts_t *sampling_opts = &(featureset->sampling_opts);
	feature_opts_t *feature_opts = &(featureset->feature_opts);

	rotations = sampling_opts->rotations;
	tiles_x = sampling_opts->tiles_x;
	tiles_y = sampling_opts->tiles_y;

	*preproc_name = '\0';
	preproc_name_lngth = 0;
	if (preproc_opts->bounding_rect.x > -1) {
		preproc_name_lngth += sprintf (preproc_name+preproc_name_lngth,"-%s%d_%d_%d_%d",
			preproc_opts->bounding_rect_base,preproc_opts->bounding_rect.x,preproc_opts->bounding_rect.y,
			preproc_opts->bounding_rect.w,preproc_opts->bounding_rect.h);
	}
	if (preproc_opts->downsample > 0 && preproc_opts->downsample < 100) {
		preproc_name_lngth += sprintf (preproc_name+preproc_name_lngth,"-%s%d",
			preproc_opts->downsample_base,preproc_opts->downsample);
	}
	if (preproc_opts->mean > -1) {
		preproc_name_lngth += sprintf (preproc_name+preproc_name_lngth,"-%s%d",
			preproc_opts->normalize_base,preproc_opts->mean);
		if (preproc_opts->stddev > -1) preproc_name_lngth += sprintf (preproc_name+preproc_name_lngth,"_%d",
			preproc_opts->stddev);
	}

	for (rot_index=0;rot_index<rotations;rot_index++) {
//			printf ("rotation:%d\n",rot_index);
		for (tile_index_x=0;tile_index_x<tiles_x;tile_index_x++) {
			for (tile_index_y=0;tile_index_y<tiles_y;tile_index_y++) {
				strcpy (sample_name,preproc_name);
				sample_name_lngth = strlen (sample_name);
				// this code block should be made more generic with regards to sampling options and params
				// sample opts
				if (rotations && rot_index) { // only add for non-0 rotations
					sample_name_lngth += sprintf (sample_name+sample_name_lngth,"-%s_%d",sampling_opts->rot_base,rot_index);
				}
				if ( (tiles_x > 1 || tiles_y > 1) ) {
					if (tiles_x == tiles_y) {
						sample_name_lngth += sprintf (sample_name+sample_name_lngth,"-%s%d_%d_%d",sampling_opts->tile_base,tiles_x,tile_index_x,tile_index_y);
					} else {
						sample_name_lngth += sprintf (sample_name+sample_name_lngth,"-%s%dx%d_%d_%d",sampling_opts->tile_base,tiles_x,tiles_y,tile_index_x,tile_index_y);
					}
				}
				// feature opts
				if (feature_opts->compute_colors) {
					sample_name_lngth += sprintf (sample_name+sample_name_lngth,"-%s",feature_opts->compute_colors_base);
				}
				if (feature_opts->large_set) {
					sample_name_lngth += sprintf (sample_name+sample_name_lngth,"-%s",feature_opts->large_set_base);
				}
				
				strcpy(featureset->samples[n_samples].sample_name,sample_name);
				featureset->samples[n_samples].rot_index = rot_index;
				featureset->samples[n_samples].tile_index_x = tile_index_x;
				featureset->samples[n_samples].tile_index_y = tile_index_y;
				n_samples++;
			}
		}
	}
	featureset->n_samples = n_samples;

	if (!feature_opts->large_set && !feature_opts->compute_colors)
		featureset->n_features = NUM_DEF_FEATURES;
	else if (!feature_opts->large_set && feature_opts->compute_colors)
		featureset->n_features = NUM_C_FEATURES;
	else if (feature_opts->large_set && !feature_opts->compute_colors)
		featureset->n_features = NUM_L_FEATURES;
	else if (feature_opts->large_set && feature_opts->compute_colors)
		featureset->n_features = NUM_LC_FEATURES;
				
}


/*
check_split_params - checks parameters for consistency with regards to training/testing a given dataset.
Returns 1 on success, 0 upon failure.
*/
int check_split_params (int *n_train_p, int *n_test_p, double *split_ratio, TrainingSet *dataset, TrainingSet *testset, int class_num, int samples_per_image, int balanced_splits, int max_training_images, int max_test_images, int exact_training_images) {
	int class_index, smallest_class=0;
	int max_balanced_samples,max_balanced_i;

	// Initialize what we will be returning
	*n_train_p = 0;
	*n_test_p = 0;

	/*
	  Bounds checking on samples.
	  -i switch is in images, but training happens on samples.
	  -t (and -R) switches provide samples per image (tiles*tiles*rotations).
	  Enforce balanced classes, with option to force unbalanced classes
	    Calculate maximum number of training images (samples) for a balanced classifier (Ib = images in the smallest class)
	    -rN: ratio assumes balanced classes, force unbalanced with -r#.  Meaning is fraction of samples to be used for *testing*.
	    -iN: Number of training images used is min (iN,Ib). i#N drops classes with less than N images, resulting in a balanced classifier.
	         -r is ignored if both -i and -r are specified.
	    -jN: Parameter is ignored if a separate test .fit is supplied.  All samples in the test .fit are used for testing in each split
	         Balanced testing is enforced (jN = min (jN,Ib-iN).  Program terminates if jN <= 0).
	         Unbalanced testing can be done by using unbalanced test samples in the test .fit
	         If -r#, testing is still balanced (jN = min (jN,Ib-iN), where iN is rN*Ib)
	*/
	
	max_balanced_samples=MAX_SAMPLES;
	for (class_index = 1; class_index <= dataset->class_num; class_index++) {
		if (dataset->class_nsamples[class_index] < max_balanced_samples) {
			max_balanced_samples = dataset->class_nsamples[class_index];
			smallest_class = class_index;
		}
	}
	max_balanced_i = max_balanced_samples / samples_per_image;
	if( verbosity >= 2 ) printf ("Max balanced training images: %d\n",max_balanced_i);
	// Check provided parameters against balanced testing/training
	if (max_training_images > 0 && !exact_training_images) { // N.B.: -i overrides -r, except if exact_training_images
		if (max_training_images > max_balanced_i && testset) {
			catError("WARNING: Specified training images (%d) exceeds maximum for balanced training (%d).\n  %d images used for training.\n  Use -r# instead of -i to over-ride balanced training.\n",
				max_training_images,max_balanced_i,max_balanced_i);
			max_training_images = max_balanced_i;
		} else if (max_training_images >= max_balanced_i && testset == NULL) { // No images left for testing unless we have a test .fit
			catError("ERROR: Specified training images (%d) exceeds maximum for balanced training (%d).  No images left for testing.\n",max_training_images,max_balanced_i-1);
			catError("  Use -iN with N < %d.\n",max_balanced_i);
			catError("  Or, use -rN with 1.0 < N > 0.0\n");
			catError("  Or, use -r#N to over-ride balanced training.\n");
			delete dataset;
			return (0);
		}
		*split_ratio = 0.0; // Don't base splits on a ratio - use max_training_images/max_test_images
	} else { // -i unspecified or used for exact_training_images, use split_ratio (default or specified - already set in main)
		max_training_images = (int)floor( (*split_ratio * (float)max_balanced_i) + 0.5 ); // rounding
		if (max_training_images >= max_balanced_i && testset == NULL) { // rounding error left no test images
			catError("ERROR: No images left for testing using specified -r=%f.\n",*split_ratio);
			catError("  Use -rN with N < %f\n", ((float)max_balanced_i - 1.0) / (float)max_balanced_i);
			catError("  Or, use -iN with N < %d.\n",max_balanced_i);
			catError("Exiting - no testing performed\n");
			delete dataset;
			return (0);
		}
		// If an exact split ratio was specified, then use it.
		// Otherwise, the same number of training images is used in each class, specified by max_training_images
		if (balanced_splits) *split_ratio = 0.0;
	}
	// Note that checks above leave max_training_images < max_balanced_i
	
	if (max_test_images > 0) { // -jN specified
		if (testset) {
			catError("WARNING: The -j%d parameter is ignored when a test set is specified (%s).\n",max_test_images,testset->source_path);
		} else if ( max_test_images > (max_balanced_i - max_training_images) ) { // -jN always balanced unless test .fit
			if (max_balanced_i - max_training_images > 0) {
				catError("WARNING: Insufficient images for balanced training (%d) and specified testing (%d).  %d images used for testing.\n",
					max_training_images,max_test_images,max_balanced_i - max_training_images);
				max_test_images = max_balanced_i - max_training_images;
			}
		}
	} else { // -jN not specified
		max_test_images = max_balanced_i - max_training_images;
		if (!balanced_splits) max_test_images = 0;
	}

	// Set the return values
	*n_train_p = max_training_images;
	*n_test_p = max_test_images;
	return (1);

}


int split_and_test(TrainingSet *ts, char *report_file_name, int argc, char **argv, int class_num, int method, featureset_t *featureset, double split_ratio, int balanced_splits, double max_features, double used_mrmr, long split_num,
	int report,int max_training_images, int exact_training_images, int max_test_images, char *phylib_path,int distance_method, int phylip_algorithm,int export_tsv,
	long first_n, char *weight_file_buffer, char weight_vector_action, int N, TrainingSet *testset, int ignore_group, int tile_areas, int max_tile, int image_similarities, int random_splits) {
	TrainingSet *train,*test,**TilesTrainingSets=NULL;
	std::vector<data_split> splits;
	char group_name[64];
	FILE *output_file;
	int split_index,tile_index;
	int n_train,n_test;
	int class_index;
	int res;
	int i;
	
	// set samples per image
	int samples_per_image = featureset->n_samples;

	// Remove classes from the end if N is specified
	if (N>0) while (ts->class_num>N) ts->RemoveClass(ts->class_num);

// Remove classes with less than max_training_images if exact_training_images is true
	if (exact_training_images) {
		class_index=ts->class_num;
		while( class_index > 0 ) {
			if( ts->class_nsamples[ class_index ] * samples_per_image <  max_training_images ) {
				ts->RemoveClass( class_index );
			}
			class_index--;
		}
	}

	// If a testset was specified, make sure its classes are consistent with ts.
	if (testset) {
		testset->train_class = new int[testset->class_num+1];
		// CEC added June 8 2012, initialize to 0 even the zeroth class,
		// since the zeroth class is associated with the unknown class
		for( int i = 0; i < testset->class_num+1; i++ )
			testset->train_class[i] = 0;

		int ts_class_index;
		for (class_index = 1; class_index <= testset->class_num; class_index++) {
			for (ts_class_index = 1; ts_class_index <= ts->class_num; ts_class_index++) {
				if ( ! strcmp(ts->class_labels[ts_class_index],testset->class_labels[class_index]) ) {
					testset->train_class[class_index] = ts_class_index;
					break;
				}
			}
			if (!testset->train_class[class_index]) {
				catError ("WARNING: Test set class label '%s' does not match any training set class.  Marked with '*'.\n",testset->class_labels[class_index]);
			}
		}
	}

// Check that the train and test sets have the same number of features
	if (testset && testset->signature_count != ts->signature_count) {
		catError ("The number of features in the train set '%s' (%d) is inconsistent with test set '%s' (%d).\n",
			ts->source_path,ts->signature_count,testset->source_path,testset->signature_count);
		return(showError(1, NULL));
	}

// Check that the train and test sets have the same feature versions
	if (testset && (testset->feature_vec_version != ts->feature_vec_version ||
		testset->feature_vec_type != ts->feature_vec_type) ) {
			catError ("ERROR: The feature versions in the train set '%s' (%d.%d) is inconsistent with test set '%s' (%d.%d).\n",
				ts->source_path, ts->feature_vec_version, ts->feature_vec_type,
				testset->source_path, testset->feature_vec_version, testset->feature_vec_type
			);
			catError ("Delete .fit and .sig files generated by older versions of wndchrm and try again.\n");
		return(showError(1, NULL));
	}
	if (testset && testset->signature_count != ts->signature_count) {
		catError ("The feature versions in the train set '%s' (%d.%d) is inconsistent with test set '%s' (%d.%d).\n",
			ts->source_path,
			ts->signature_count,testset->source_path,testset->signature_count);
		return(showError(1, NULL));
	}

// Check the parameters and set train/test image numbers
	if (!check_split_params (&n_train, &n_test, &split_ratio, ts, testset,
		class_num, samples_per_image, balanced_splits, max_training_images, max_test_images, exact_training_images))
			return(showError(1, NULL));

	// Instantiate a featuregroup_t for the top level ts
	// to keep track of feature statistics across splits.
	if( split_num > 1 ) 
		ts->aggregated_feature_stats = new featuregroups_t;

/*
	In ts-split(),
      if (test->count > 0) number_of_test_samples = 0; // test already has samples from a file
      // if ratio is 0, use the max_train_samples and max_test_samples (balanced training and testing)
      // if ratio is > 0 and <= 1, use ratio (unbalanced training)      
*/
	if (verbosity>=2) {
		if (split_ratio > 0) printf ("samples per image=%d, UNBALANCED training fraction=%g\n",samples_per_image,split_ratio);
		else printf ("samples per image=%d, training images: %d, testing images %d\n",samples_per_image,n_train,n_test);
	}
	splits.resize(split_num);
	for (split_index=0;split_index<split_num;split_index++)
	{
		double accuracy;
		double feature_weight_distance=-1.0;

		train=new TrainingSet(ts->count,ts->class_num);
		if (testset) test = testset;
		else test=new TrainingSet(ts->count,ts->class_num);
		splits[split_index].confusion_matrix=new unsigned short[(ts->class_num+1)*(ts->class_num+1)];
		splits[split_index].training_images=new unsigned short[(ts->class_num+1)];
		splits[split_index].testing_images=new unsigned short[(ts->class_num+1)];
		splits[split_index].class_accuracies=new double[(ts->class_num+1)];
		splits[split_index].similarity_matrix=new double[(ts->class_num+1)*(ts->class_num+1)];
		splits[split_index].class_probability_matrix=new double[(ts->class_num+1)*(ts->class_num+1)];
		if (tile_areas) {
			splits[split_index].tile_area_accuracy=new double[samples_per_image];
			for (tile_index=0;tile_index<samples_per_image;tile_index++) splits[split_index].tile_area_accuracy[tile_index]=0.0;
		}
		else splits[split_index].tile_area_accuracy=NULL;

		res=ts->split(random_splits,split_ratio,train,test,samples_per_image,n_train,n_test,&(splits[split_index]));
		if ( res < 0) return (res);
		if (image_similarities) splits[split_index].image_similarities=new double[(1+test->count/(samples_per_image))*(1+test->count/(samples_per_image))];
		else splits[split_index].image_similarities=NULL;

		//int temp=train->class_num;
		//train->class_num=1;
		if (tile_areas)  // split into several datasets such that each dataset contains tiles of the same location
		{
			TilesTrainingSets=new TrainingSet*[samples_per_image];
			res = train->SplitAreas(samples_per_image, TilesTrainingSets);
			if (res < 0) return (res);
			for (tile_index=0;tile_index<samples_per_image;tile_index++)
			{
				TilesTrainingSets[tile_index]->normalize();
				TilesTrainingSets[tile_index]->SetFisherScores(max_features,used_mrmr,NULL);
			}
		}
		else
		{
			train->normalize(); // normalize the feature values of the training set
			train->SetFisherScores(max_features,used_mrmr,&(splits[split_index]));  // compute the Fisher Scores for the image features
			if( ts->aggregated_feature_stats ) {
				if( ts->aggregated_feature_stats->empty() ) {
					featuregroup_stats_t temp;
					for( i = 0; i < ts->signature_count; i++ ) {
						temp.name = train->SignatureNames[i];
						temp.min = train->SignatureWeights[i];
						temp.max = train->SignatureWeights[i];
						temp.sum_weight = train->SignatureWeights[i];
						temp.sum_weight2 = train->SignatureWeights[i] * train->SignatureWeights[i];
						temp.mean = 0;
						temp.stddev = 0;
						temp.n_features = 1;	
						ts->aggregated_feature_stats->push_back( temp ); // makes a copy
					}
				}
				else
				{
					for( i = 0; i < ts->signature_count; i++ ) {
						if( train->SignatureWeights[i] < (*(ts->aggregated_feature_stats))[i].min )
							 (*(ts->aggregated_feature_stats))[i].min = train->SignatureWeights[i];
						if( train->SignatureWeights[i] > (*(ts->aggregated_feature_stats))[i].max )
							 (*(ts->aggregated_feature_stats))[i].max = train->SignatureWeights[i];
						(*(ts->aggregated_feature_stats))[i].sum_weight += train->SignatureWeights[i];
						(*(ts->aggregated_feature_stats))[i].sum_weight2 += train->SignatureWeights[i] * train->SignatureWeights[i];
						(*(ts->aggregated_feature_stats))[i].n_features++;
					}
				}
			}
		}
		//train->class_num=temp;
		if (weight_vector_action=='w')
			if(!train->SaveWeightVector(weight_file_buffer))
				showError(1,"Could not write weight vector to '%s'\n",weight_file_buffer);
		if (weight_vector_action=='r' || weight_vector_action=='+' || weight_vector_action=='-')
		{
			feature_weight_distance=train->LoadWeightVector(weight_file_buffer,(weight_vector_action=='+')-(weight_vector_action=='-'));
			if (tile_areas) for (tile_index=0;tile_index<samples_per_image;tile_index++) feature_weight_distance=TilesTrainingSets[tile_index]->LoadWeightVector(weight_file_buffer,(weight_vector_action=='+')-(weight_vector_action=='-'));	   
			if (feature_weight_distance<0) showError(1,"Could not load weight vector from '%s'\n",weight_file_buffer);
		}
		if (report) splits[split_index].individual_images=new char[(int)((test->count/(samples_per_image))*(class_num*15))];
		else splits[split_index].individual_images=NULL;
		if (ignore_group)   /* assign to zero all features of the group */
		{
			if (!(ts->IgnoreFeatureGroup(ignore_group,group_name))) {
				delete train;
				if (!testset) delete test;
				delete splits[split_index].confusion_matrix;
				delete splits[split_index].training_images;
				delete splits[split_index].testing_images;
				delete splits[split_index].class_accuracies;
				delete splits[split_index].similarity_matrix;
				delete splits[split_index].individual_images;
				showError(1,"Errors while trying to ignore group %d '%s'\n",ignore_group,group_name);
				return(0);
			}
		}
		
		// The following is the separator between split results
		printf( "\n----------\n" );

		// Label the columns
		if (verbosity>=1) {
			printf("image\t");
			if (ts->is_continuous) {
				printf("act. val.\tpred. val.\n");
			} else {
				printf ("norm. fact.\t");
				for (class_index=1;class_index<=ts->class_num;class_index++) {
					printf("p(%s)\t",ts->class_labels[class_index]);
				}
				printf("act. class\tpred. class");
				if (ts->is_numeric) printf ("\tpred. val.");
				printf ("\n");
			}
		}

		accuracy=train->Test(test,method,samples_per_image,tile_areas,TilesTrainingSets,max_tile,first_n,&(splits[split_index]));

		splits[split_index].feature_weight_distance=feature_weight_distance;
		splits[split_index].method=method;
		splits[split_index].pearson_coefficient=test->pearson(samples_per_image,&(splits[split_index].avg_abs_dif),&(splits[split_index].pearson_p_value));

		if (!report && !ignore_group && verbosity > 2 )   // print the accuracy and confusion and similarity matrices
		{ 
			printf( "\n" );
			ts->PrintConfusion(stdout,splits[split_index].confusion_matrix,NULL);//,0,0);
			ts->PrintConfusion(stdout,NULL,splits[split_index].similarity_matrix);//,0,0);
			if (ts->is_continuous) printf("Pearson Correlation: %f \n\n",splits[split_index].pearson_coefficient);
			else printf("\nAccuracy: %f \n",accuracy);
		}

		if (TilesTrainingSets)    // delete the training sets allocated for the different areas
		{
			for (tile_index=0;tile_index<samples_per_image;tile_index++)
			delete TilesTrainingSets[tile_index];
			delete TilesTrainingSets;
		}
		delete train;
		if (!testset) delete test;
	} // End for (split_index=0;split_index<split_num;split_index++)

	if( ts->aggregated_feature_stats ) {
		// Finish up computing the averages for the feature weights
		for(featuregroups_t::iterator avgs_it = ts->aggregated_feature_stats->begin(); avgs_it != ts->aggregated_feature_stats->end(); ++avgs_it )
		{
			avgs_it->mean = avgs_it->sum_weight / (double)(avgs_it->n_features);
			avgs_it->stddev = sqrt (
				(avgs_it->sum_weight2 - (avgs_it->sum_weight * avgs_it->mean)) / (double)(avgs_it->n_features - 1));
		}
		// Sort the aggregated feature statistics vector by mean weight
		sort_by_mean_weight_t sort_by_mean_weight_func;
		sort( ts->aggregated_feature_stats->begin(), ts->aggregated_feature_stats->end(), sort_by_mean_weight_func);
		// Lop off the vector after a threshold
		// Right now, we only care about the top 50 features
		ts->aggregated_feature_stats->erase (ts->aggregated_feature_stats->begin() + 50, ts->aggregated_feature_stats->end());
	}

	// if( verbosity >= 2 ) printf("\n\n");
	if (!report && verbosity != 1 )	{
		printf( "\n----------\n" );

		// print the average accuracy
		double avg_accuracy=0,avg_pearson=0;
		int split_index;
		for( split_index = 0; split_index < split_num; split_index++ )
		{
			avg_accuracy+=splits[split_index].accuracy;
			avg_pearson+=splits[split_index].pearson_coefficient;
		}
		if (ignore_group) printf("Accuracy assessment without using feature group '%s' - ",group_name); 
		if (ts->is_continuous) printf("Average Pearson Correlation (%ld splits): %f\n",split_num,avg_pearson/(double)split_num);
		else printf("Average accuracy (%ld splits): %f\n",split_num,avg_accuracy/(double)split_num);
		
		// print out averages across all splits
		int length = (1+ts->class_num) * (1+ts->class_num);
		short unsigned int *confusion_matrix = new short unsigned int[ length ];
		double *avg_similarity_matrix = new double[ length ];
		double *avg_class_probability_matrix = new double[ length ];
		for( i = 0; i < length; ++i ) {
			confusion_matrix[ i ] = 0;
			avg_similarity_matrix[ i ] = 0;
			avg_class_probability_matrix[ i ] = 0;
		}

		for( int row = 1; row <= ts->class_num; ++row )
		{
			for( int col = 1; col <= ts->class_num; ++col )
			{
				int confus_sum = 0;
				double avg_siml = 0.0;
				double avg_class_prob_sum = 0.0;
				for( split_index = 0; split_index < split_num; ++split_index )
				{
					confus_sum += splits[split_index].confusion_matrix[ row * ts->class_num + col ];
					avg_siml += splits[split_index].similarity_matrix[ row * ts->class_num + col ];
					avg_class_prob_sum += splits[split_index].class_probability_matrix[ row * ts->class_num + col ];
				}
				confusion_matrix[ row * ts->class_num + col ] = confus_sum;
				avg_similarity_matrix[ row * ts->class_num + col ] = avg_siml / split_num;
				avg_class_probability_matrix[ row * ts->class_num + col ] = avg_class_prob_sum / split_num;
			}
		}
		printf("\nConfusion Matrix (sum of all splits)\n");
		ts->PrintConfusion(stdout, confusion_matrix, NULL);
		printf("\nAverage Similarity Matrix\n");
		ts->PrintConfusion(stdout,NULL, avg_similarity_matrix);
		printf("\nAverage Class Probability Matrix\n");
		ts->PrintConfusion(stdout,NULL, avg_class_probability_matrix);
		printf("\n----------\n");
	}

	if (report)
	{
		if (report_file_name){
			if (!strchr(report_file_name,'.')) strcat(report_file_name,".html");
			output_file=fopen(report_file_name,"w");
			if (!output_file) showError(1, "Could not open file for writing '%s'\n",report_file_name);
		}
		else output_file=stdout;     
		ts->report(output_file,argc,argv,report_file_name,splits,split_num,featureset,n_train,
				phylib_path, distance_method, phylip_algorithm,export_tsv,
				testset,image_similarities, bool(random_splits) );
		if (output_file!=stdout) fclose(output_file);
		// copy the .ps and .jpg of the dendrogram to the output path of the report and also copy the tsv files
		if (export_tsv || phylib_path) {
			char command_line[512],ps_file_path[512];
			strcpy(ps_file_path,report_file_name);
			if ( strrchr(ps_file_path,'/') ) {
				(strrchr(ps_file_path,'/'))[1]='\0';
				if (phylib_path) {
					sprintf(command_line,"mv ./%s*.ps %s",ts->name,ps_file_path);
					system(command_line);
					sprintf(command_line,"mv ./%s*.jpg %s",ts->name,ps_file_path);
					system(command_line);
				}

				if (export_tsv) {
					sprintf(command_line,"cp -r ./tsv %s",ps_file_path);
					system(command_line);		  
					sprintf(command_line,"rm -r ./tsv");
					system(command_line);
				}
			}
		}
	}

	for (split_index=0;split_index<split_num;split_index++) {
		delete splits[split_index].confusion_matrix;	
		delete splits[split_index].training_images;	
		delete splits[split_index].testing_images;	
		delete splits[split_index].class_accuracies;	
		delete splits[split_index].similarity_matrix;
		if (splits[split_index].individual_images) delete splits[split_index].individual_images;
		if (splits[split_index].tile_area_accuracy) delete splits[split_index].tile_area_accuracy;
		if (splits[split_index].image_similarities) delete splits[split_index].image_similarities;	   
	}

	if( ts->aggregated_feature_stats ) delete ts->aggregated_feature_stats;
	return(1);
}


void ShowHelp()
{
	printf("\n"PACKAGE_STRING".  Laboratory of Genetics/NIA/NIH \n");
	printf("usage: \n======\nwndchrm [ train | test | classify ] [-mtslcdowfrijnpqvNSBACDTh] [<dataset>|<train set>] [<test set>|<feature file>] [<report_file>]\n");
	printf("  <dataset> is a <root directory>, <feature file>, <file of filenames>, <image directory> or <image filename>\n");
	printf("  <root directory> is a directory of sub-directories containing class images with one class per sub-directory.\n");
	printf("      The sub-directory names will be used as the class labels. Currently supported file formats: TIFF, PPM. \n");
	printf("  <feature file> is the file generated by the train command containing all computed image features (should end in .fit).\n");
	printf("       This filename is a required parameter for storing the output of 'train'\n");       
	printf("  <file of filenames> is a text file listing <image filename>s and corresponding class labels\n");
	printf("      separated by a <TAB> character (a tab delimited file, or .tsv). Lines beginning with '#' are ignored\n");       
	printf("  <image directory> is a directory of image files. The class cannot be specified so these can only be used as a <test set>.\n");    
	printf("  <image filename> is the full path to an image file. The classes cannot be specified so these can only be used as a <test set>.\n");    
	printf("  <train set> is anything that qualifies as a <dataset>, but must contain at least two (2) defined classes.\n");
	printf("      An <image filename> or <image directory> cannot define classes.\n");
	printf("  <test set> is anything that qualifies as a <dataset>.  The <train set> will be used to classify the <test set>.\n");
	printf("      This parameter is required for 'classify' and is optional for 'test' (when doing internal tests of the <train set>\n");
	printf("  <report_file> is a report of the test/classify results in html format (must end in .htm or .html).\n");
	
	printf("\nImage sampling options (require re-computing features):\n========================================================\n");
	printf("m - Allow running multiple instances of this program concurrently, save (and re-use) pre-calculated .sig files.\n");
	printf("    This will save and re-use .sig files, making this option useful for single instances/processors as well\n");
	printf("R - Add rotations to training images (0,90,180,270 degrees).\n");
	printf("t[#][^]C[xR] - split the image into CxC or CxR tiles if R is specified. The default is 1.\n");
	printf("    If '#' is specified, each tile location is used as a separate dataset (for testing only!). \n");
	printf("      - If both '#' and '^' are specified only the closest tile is used. \n");
	printf("    If only C is specified (e.g. -t2), tiling will be CxC (e.g. 2 columns by 2 rows). \n");
	printf("      - If both C and R are specified (e.g. -t2x3), tiling will be CxR (e.g. 2 columns by 3 rows). \n");
	printf("dN - Downsample the images (N percents, where N is 1 to 100)\n");
	printf("Sx[:y] - normalize the images such that the mean is set to x and (optinally) the stddev is set to y.\n");   
	printf("Bx,y,w,h - compute features only from the (x,y,w,h) block of the image.\n");      
	
	printf("\nImage Feature options:\n======================\n");
	printf("l - Use a large image feature set.\n");
	printf("c - Compute color features.\n");
	printf("o - force overwriting pre-computed .sig files.\n");   
	printf("O - if there are pre-computed .sig files accompanying images that have the old-style naming pattern,\n" );
	printf("    skip the check to see that they were calculated with the same wndchrm parameters as the current experiment.\n");   
	
	printf("\nFeature reduction options:\n==========================\n");
	printf("fN[:M] - maximum number of features out of the dataset (0,1) . The default is 0.15. \n");
	printf("v[r|w|+|-][path] - read/write/add/subtract the feature weights from a file.\n");   
	printf("A - assess the contribution of each group of image features independently.\n");
	
	printf("\nClassifier options:\n===================\n");
	printf("w - Classify with wnn instead of wnd. \n");
	printf("qN - the number of first closest classes among which the presence of the right class is considered a match.\n");
	printf("r[#]N - Fraction of images/samples to be used for training (0,1). The default is 0.75 of\n");
	printf("        the smallest class. if '#' is specified, force unbalanced training\n");
	printf("i[#]N - Set a maximal number of training images (for each class). If the '#' is specified then\n");
	printf("        the class is ignored if it doesn't have at least N samples.\n");
	printf("jN - Set a maximal number of test images (for each class). \n");
	printf("nN - Number of repeated random splits. The default is 1.\n");
	printf("Nx - set the maximum number of classes (use only the first x classes).\n");
	//printf("C - *highly experimental* perform interpolation on a continuous scale rather than discrete classes\n");
	//printf("    All class labels must be interpretable as numbers.\n");
	
	printf("\nOutput options:\n===============\n");
	printf("s - silent mode. Optionally followed by a verbosity level (higher = more verbose)\n");
	printf("p[+][k][#][path] - Report options.\n");
	printf("   'path' is an optional path to a PHYLIP installation root directory for generating dendrograms.\n");
	printf("   The optinal '+' creates a 'tsv' directory and exports report data into tsv files.\n");
	printf("   'k' is an optional digit (1..3) of the specific phylip algorithm to be used.\n");
	printf("   '#' generates a similarity map of the test images\n");
	printf("P[N] - pair-wise distance algorithms for comparing classes\n");
	printf("   The class probability matrix is the average of marginal probabilities for the images in each class \n");
	printf("   The similarity matrix is the class probability matrix, where each row is normalized to make the class identity column equal to 1.0\n");
	printf("   The dis-similarity (i.e. 1.0 - similarity) between two classes can be interpreted as a \"morphological distance\".\n");
	printf("   There are two entries in the similarity matrix for each comparison: Class 1 classified as Class 2, and Class 2 classified as Class 1.\n");
	printf("   N = 1: Use the maximum of the two dis-similarities.\n");
	printf("   N = 2: Use the average of the two dis-similarities.\n");
	printf("   N = 3: Use the top triangle only (i.e. Class 1 classified as Class 2)\n");
	printf("   N = 4: Use the bottom triangle only (i.e. Class 2 classified as Class 1)\n");
	printf("   N = 5: Use the class probability matrix as a set of coordinates for each class centroid in a \"Marginal Probability Space\". Use Euclidean distance.\n");
	printf("   The default method is 5. Method 2 was described in ref [1], and method 5 was described in ref [2].\n");
	printf("D[path] - feature file name (.fit file) to save the <dataset> or <train set>.\n");
	printf("T[path] - feature file name (.fit file) to save the <test set>.\n");
	printf("h - show this note.\n");
	
	printf("\nExamples:\n=========\n");
	printf("train:\n");
	printf("  wndchrm train /path/to/dataset/ dataset.fit\n");
	printf("  wndchrm train -mcl /path/to/dataset/ testset.fit\n");
	printf("test:\n");
	printf("  wndchrm test -f0.1 dataset.fit\n");
	printf("  wndchrm test -f0.1 -r0.9 -n5 dataset.fit testset.fit\n");
	printf("  wndchrm test -f0.2 -i50 -j20 -n5 -p/path/to/phylip3.65 dataset.fit testset.fit report.html\n");
	printf("  N.B.: By default, the -r or -i parameters will be used to make a balanced training set (equal number of images per class).\n");
	printf("       -r#N can be used to override this default, so that the N fraction of each class will be used for training.\n");
	printf("       If a <test set> is specified, it will be used as the test set for each 'split', but training images will\n");
	printf("       still be randomly chosen from <train set>)\n");
	printf("classify:\n");
	printf("   wndchrm classify dataset.fit /path/to/image.tiff\n");
	printf("   wndchrm classify -f0.2 -cl /path/to/root/dir /path/to/image/directory/\n");
	printf("   wndchrm classify -f0.2 -cl -Ttestset.fit dataset.fit /path/to/image/file_of_filenames.tsv\n");
	printf("   N.B.: classify will use -r or -i to train with fewer than all of the images in <dataset>\n");
	printf("       Unlike 'test', 'classify' will chose the training images in order rather than randomly.\n");
	printf("       classify will ignore the -n parameter because the result will be the same for each run or split.\n");
	printf("       The default -r for 'classify' is 1.0 rather than the 0.75 used in 'test'.\n");
	printf("\nAdditional help:\n================\n");
	printf("A detailed description can be found in: Shamir, L., Orlov, N., Eckley, D.M., Macura, T., Johnston, J., Goldberg, I.\n");
	printf("  [1] \"Wndchrm - an open source utility for biological image analysis\", BMC Source Code for Biology and Medicine, 3:13, 2008.\n");   
	printf("An application of pattern recognition for a quantitative biological assay based on morphology can be found in:\n");
	printf("  [2] Johnston, J., Iser W. B., Chow, D. K., Goldberg, I. G., Wolkow, C. A. \"Quantitative Image Analysis Reveals\n");
	printf("  Distinct Structural Transitions during Aging in Caenorhabditis elegans Tissues\", PLoS ONE, 3:7:e2821, 2008.\n");
	
	printf("\nIf you have questions or problems with this software, please visit our Google code page <http://code.google.com/p/wnd-charm/> \n\n");
	return;
}


int main(int argc, char *argv[])
{   char *dataset_path=NULL, *testset_path=NULL;
    int multi_processor=0;
    int arg_index=1;
    int tile_areas=0;
    int method=1;
    int report=0;
    int splits_num=1;
    double split_ratio=0.75;
    double max_features=0.15;
    double used_mrmr=0.0;
    int max_training_images=0;
    int max_test_images=0;
    int train=0;
    int test=0;
    int classify=0;
    char phylib_path_buffer[256];
    char *phylib_path=NULL;
    char report_file_buffer[256];
    char *report_file=NULL;
    int export_tsv=0;
    int phylip_algorithm=0;
    int distance_method=5;
    int exact_training_images=0;
    long first_n=1;
    char weight_file_buffer[256];
    char weight_vector_action='\0';
    int N=0;                         /* use only the first N classes                               */
    int assess_features=0;           /* assess the contribution of each feature to the performance */
    int image_similarities=0;        /* generate a dendrogram showing the similarity of the images */
    int max_tile=0;                  /* use only the closest tile                                  */
	int overwrite=0;                 /* force overwriting of pre-computed .sig files               */
	char *dataset_save_fit=NULL;     /* path to save the dataset/train set                         */
	char *testset_save_fit=NULL;     /* path to save the test set                                  */
	char *char_p,*char_p2;
	int balanced_splits=1;           /* when 1, use balanced training.  Override with -r#          */
	int random_splits=1;             /* when 1 randomly chose training images, when 0 add them in read order */
	int do_continuous=0;
	int save_sigs=1;
	int skip_sig_check = 0;

	assert (ComputationTaskInstances::initialized() && "Failed to initialize computation tasks");

	featureset_t featureset;         /* for recording the sampling params for images               */
	memset (&featureset,0,sizeof(featureset));
	preproc_opts_t *preproc_opts = &(featureset.preproc_opts);
	sampling_opts_t *sampling_opts = &(featureset.sampling_opts);
	feature_opts_t *feature_opts = &(featureset.feature_opts);

// initialize the param structs
	strcpy (preproc_opts->bounding_rect_base,"B");
	preproc_opts->bounding_rect.x = -1;            /* a bounding rect from which features should be computed     */
	preproc_opts->bounding_rect.y = -1;
	preproc_opts->bounding_rect.w = -1;
	preproc_opts->bounding_rect.h = -1;
	strcpy (preproc_opts->downsample_base,"d");
	preproc_opts->downsample = 100;
	strcpy (preproc_opts->normalize_base,"S");
	preproc_opts->mean = -1;                      /* normalize all image to a sepcified mean                    */
	preproc_opts->stddev = -1;                  /* normalize all image to a sepcified standard deviation      */

	strcpy (sampling_opts->rot_base,"R");
	sampling_opts->rotations = 1;
	strcpy (sampling_opts->tile_base,"t");
	sampling_opts->tiles_x = 1;
	sampling_opts->tiles_y = 1;

	strcpy (feature_opts->compute_colors_base,"c");
	feature_opts->compute_colors = 0;
	strcpy (feature_opts->large_set_base,"l");
	feature_opts->large_set = 0;


    /* read parameters */
    if (argc<2)
    {  ShowHelp();
       return(1);
    }

    if (strcmp(argv[arg_index],"train")==0) train=1;
    if (strcmp(argv[arg_index],"test")==0) {
    	test=1;
    	split_ratio = 0.75;
    	random_splits = 1;
    }
    if (strcmp(argv[arg_index],"classify")==0) {
    	classify=1;
    	split_ratio = 1.0;
    	random_splits = 0; // use order in the input file
    }
	if (!train && !test && !classify) {
		ShowHelp();
		showError(1,"Either 'train', 'test' or 'classify' must be specified.\n");
		return(1);
	}
    arg_index++;

	/* read the switches */
    while (argv[arg_index][0]=='-')
    {   char *p,arg[32];
	    if (argv[arg_index][1]=='p')
        {  report=1;
		   if ((strchr(argv[arg_index],'p')[1])=='+') export_tsv=1;
		   if (isdigit(strchr(argv[arg_index],'p')[1+export_tsv])) phylip_algorithm=(strchr(argv[arg_index],'p')[1+export_tsv])-'0';
		   image_similarities=((strchr(argv[arg_index],'p')[1+export_tsv+(phylip_algorithm>0)])=='#');
           if ((strchr(argv[arg_index],'p')[1+export_tsv+image_similarities+(phylip_algorithm>0)])=='/' || (strchr(argv[arg_index],'p')[1+export_tsv+image_similarities+(phylip_algorithm>0)])=='.')
		   {   strcpy(phylib_path_buffer,&(strchr(argv[arg_index],'p')[1+export_tsv+image_similarities+(phylip_algorithm>0)]));
               phylib_path=phylib_path_buffer;
		   }
		   if (phylip_algorithm<=0) phylip_algorithm=1;   /* set the default */
		   arg_index++;
		   continue;	/* so that the path will not trigger other switches */
        }
        if (strchr(argv[arg_index],'P')) distance_method=atoi(&(strchr(argv[arg_index],'P')[1]));
		if (argv[arg_index][1]=='v' && strlen(argv[arg_index])>3)
		{  weight_vector_action=argv[arg_index][2];
		   if (weight_vector_action!='r' && weight_vector_action!='w' && weight_vector_action!='+' && weight_vector_action!='-')
		     showError(1,"Unspecified weight vector action (-v switch)\n");
		   strcpy(weight_file_buffer,&(strchr(argv[arg_index],'v')[2]));
		   arg_index++;
		   continue;   /* so that the path will not trigger other switches */
		}
		if (argv[arg_index][1]=='D') {
			dataset_save_fit = argv[arg_index]+1;
	    	arg_index++;
			continue;	/* so that the path will not trigger other switches */
		}
		if (argv[arg_index][1]=='T') {
			testset_save_fit = argv[arg_index]+1;
	    	arg_index++;
			continue;	/* so that the path will not trigger other switches */
		}
        /* a block for computing features */
        if ( (char_p = strchr(argv[arg_index],'B'))  && isdigit (*(char_p+1)) ) {
			strcpy(arg,char_p+1);
			p=strtok(arg," ,;");
			preproc_opts->bounding_rect.x=atoi(p);
			p=strtok(NULL," ,;");
			preproc_opts->bounding_rect.y=atoi(p);
			p=strtok(NULL," ,;");
			preproc_opts->bounding_rect.w=atoi(p);
			p=strtok(NULL," ,;");
			preproc_opts->bounding_rect.h=atoi(p);
		}
        /* mean and stabndard deviation for normalizing the images */
        if (strchr(argv[arg_index],'S'))   
        {  strcpy(arg,argv[arg_index]);
		   p=strchr(arg,':');
		   if (p)                          /* standard deviation is specified */
		   {  preproc_opts->stddev=atoi(&(p[1]));
              *p='\0';
		   }
		   preproc_opts->mean=atoi(&(strchr(arg,'S')[1]));   /* mean */
        }
	    if (strchr(argv[arg_index],'m')) multi_processor=1;
        if (strchr(argv[arg_index],'n')) splits_num=atoi(&(strchr(argv[arg_index],'n')[1]));
        if( (char_p = strchr( argv[arg_index],'s') ) ) {
			if( isdigit( *(char_p+1) ) ) {
				verbosity = atoi( char_p+1 );
			} else {
				verbosity = 0;
			}
		}
        if (strchr(argv[arg_index],'o')) overwrite=1;
        if (strchr(argv[arg_index],'O')) skip_sig_check=1;
        if (strchr(argv[arg_index],'l')) feature_opts->large_set=1;
        if (strchr(argv[arg_index],'c')) feature_opts->compute_colors=1;
        if (strchr(argv[arg_index],'C')) do_continuous=1;
        if (strchr(argv[arg_index],'d')) preproc_opts->downsample=atoi(&(strchr(argv[arg_index],'d')[1]));
        if ( (char_p = strchr(argv[arg_index],'f')) ) {
            if ( (char_p2=strchr(char_p,':')) ) used_mrmr=atof(char_p2+1);
		    max_features=atof(char_p+1);
		}
        if ( (char_p = strchr(argv[arg_index],'r')) ) {
			if (*(char_p+1)=='#') {
				balanced_splits = 0;
				char_p++;
			}
           split_ratio=atof(char_p+1);
        }
        if (strchr(argv[arg_index],'q')) first_n=atoi(&(strchr(argv[arg_index],'q')[1]));
        if (strchr(argv[arg_index],'N')) N=atoi(&(strchr(argv[arg_index],'N')[1]));
        if (strchr(argv[arg_index],'A')) assess_features=200; 
        if (strchr(argv[arg_index],'R')) sampling_opts->rotations=4;
        if ( (char_p = strchr(argv[arg_index],'t')) ) {
			if (*(char_p+1)=='#') {
				tile_areas = 1;
				char_p++;
				if (*(char_p+1)=='^') {
					max_tile = 1;
					char_p++;
				}
			}
			sampling_opts->tiles_x = strtol(char_p+1, &char_p2, 0);
			if (*char_p2 == 'x' || *char_p2 == 'X') {
				char *char_p3;
				sampling_opts->tiles_y = strtol(char_p2+1, &char_p3, 0);
				if (char_p3 <= char_p2+1) sampling_opts->tiles_y = sampling_opts->tiles_x;
			} else if (char_p2 > char_p+1) { // no 'x', but got a number after -t: y same as x
				sampling_opts->tiles_y = sampling_opts->tiles_x;
			} else if (char_p2 <= char_p+1) { // no t parameter is an error
				showError(1,"Unspecified tiling scheme (-t switch)\n");
			}
			if (! (sampling_opts->tiles_x > 0 && sampling_opts->tiles_y > 0) ) {
				showError(1,"Badly defined tiling scheme (-t switch). Tiles in both directions must be > 0\n");
			}
       }
        if ( (char_p = strchr(argv[arg_index],'i')) ) {
			if (*(char_p+1)=='#') {
				exact_training_images = 1;
				char_p++;
			}
           max_training_images=atoi(char_p+1);
        }
        if (strchr(argv[arg_index],'j')) max_test_images=atoi(&(strchr(argv[arg_index],'j')[1]));
        if (strchr(argv[arg_index],'w')) method=0;
        if (strchr(argv[arg_index],'h'))
        {  ShowHelp();
           return(1);
        }
        if (strchr(argv[arg_index],'P')) distance_method=atoi(&(strchr(argv[arg_index],'P')[1]));


        arg_index++;
     }

/* check that the values in the switches are correct */
	if (test && splits_num<=0) showError(1,"splits number (n) must be an integer greater than 0");
	if (test && max_training_images<0) showError(1,"Maximal number of training images (i) must be an integer greater than 0");
	if (test && max_test_images<0) showError(1,"maximal number of test images (j) must be an integer greater than 0");
	if (test && report && arg_index==argc-1) showError(1,"a report html file must be specified");
	if (sampling_opts->tiles_x<=0 || sampling_opts->tiles_y <=0) showError(1,"number of tiles (t) must be an integer greater than 0");
	if (preproc_opts->downsample<1 || preproc_opts->downsample>100) showError(1,"downsample size (d) must be an integer between 1 to 100");
	if (split_ratio<0 || split_ratio>1) showError(1,"training fraction (r) must be > 0 and < 1");
	if (splits_num<1 || splits_num>MAX_SPLITS) showError(1,"splits num out of range");
	if (weight_vector_action!='\0' && weight_vector_action!='r' && weight_vector_action!='w' && weight_vector_action!='-' && weight_vector_action!='+') showError(1,"-v must be followed with either 'w' (write) or 'r' (read) ");
	if (distance_method < 1 || distance_method > 5) showError(1,"Unrecognized distance method %d.  Must be between 1 and 5.",distance_method);

	 /* run */
	randomize();   /* random numbers are used for selecting random samples for testing and training */
	setup_featureset (&featureset);
	if (arg_index<argc) {
		int res;
		dataset_path=argv[arg_index++];
		TrainingSet *dataset=new TrainingSet(MAX_SAMPLES,MAX_CLASS_NUM);

		if (train) {
			if (!dataset_save_fit && arg_index < argc && argv[arg_index] && *(argv[arg_index])) dataset_save_fit = argv[arg_index];
			else if (!dataset_save_fit) showError (1,"No output file specified");
				
		// Make sure we can write to the output file before calculating anything.
			FILE *out_file;
			if (!(out_file=fopen(dataset_save_fit,"a"))) {// don't truncate if exists.
				showError (1,"Couldn't open '%s' for writing\n",dataset_save_fit);
				return(0);
			}
			fclose (out_file);
			res=dataset->LoadFromPath(dataset_path, save_sigs, &featureset, do_continuous, skip_sig_check );
			if (res < 1) showError(1,"Errors reading from '%s'\n",dataset_path);
			res = dataset->SaveToFile (dataset_save_fit);
			if (res < 1) showError (1,"Could not save dataset to '%s'.\n",dataset_save_fit);
			if (verbosity>=2) printf ("Saved dataset to '%s'.\n",dataset_save_fit);
	
			// report any warnings
			showError (0,NULL);

       } else if (test || classify) {
			int ignore_group=0;
			TrainingSet *testset=NULL;
		// Make sure we can write to the dataset output file if we got one before calculating anything.
			FILE *out_file;
			if (dataset_save_fit && !(out_file=fopen(dataset_save_fit,"a"))) {// don't truncate if exists.
				showError (1,"Couldn't open '%s' for writing\n",dataset_save_fit);
				return(0);
			} else if (dataset_save_fit) fclose (out_file);
	
			/* check if there is a test set feature file */
			if (arg_index<argc && strstr(argv[arg_index],".htm")==NULL) testset_path=argv[arg_index++];
			if (classify && !testset_path) {
				showError (1,"The classify command must have a test set to work on.\n");
			}


			if (verbosity>=2) printf ("Processing training set '%s'.\n",dataset_path);
			res=dataset->LoadFromPath(dataset_path, save_sigs, &featureset, do_continuous, skip_sig_check);
			if (res < 1) showError(1,"Errors reading from '%s'\n",dataset_path);
			if (dataset_save_fit) {
				res = dataset->SaveToFile (dataset_save_fit);
				if (res < 1) showError (1,"Could not save dataset to '%s'.\n",dataset_save_fit);
				if (verbosity>=2) printf ("Saved dataset to '%s'.\n",dataset_save_fit);
			}

			/* check if there is a report file name */
			if (arg_index<argc) {
				strcpy(report_file_buffer,argv[arg_index]);
				report_file=report_file_buffer;
				report=1;   /* assume that the user wanted a report if a report file was specified */
			}
			
			// Load the test set if there is one
			if (testset_path) {
				if (testset_save_fit && !(out_file=fopen(testset_save_fit,"a"))) {// don't truncate if exists.
					showError (1,"Couldn't open '%s' for writing\n",testset_save_fit);
					return(0);
				} else if (testset_save_fit) fclose (out_file);
				if (verbosity>=2) printf ("Processing test set '%s'.\n",testset_path);
				testset=new TrainingSet(MAX_SAMPLES,MAX_CLASS_NUM);
				res=testset->LoadFromPath(testset_path, save_sigs, &featureset, do_continuous, skip_sig_check);
				if (res < 1) showError(1,"Errors reading from '%s'\n",testset_path);
				if (testset_save_fit) {
					res = testset->SaveToFile (testset_save_fit);
					if (res < 1) showError (1,"Could not save testset to '%s'.\n",testset_save_fit);
					if (verbosity>=2) printf ("Saved testset to '%s'.\n",testset_save_fit);
				}
			}
			if (classify) {
				if (splits_num > 1) catError ("WARNING: -n option is ignored for 'classify'.  Results are based on a single test because there is no randomization.\n");
				splits_num = 1;
				random_splits = 0;
				// defaults are different (-r = 1.0 and random_splits = 0.  Set above, though -r can still be modified).
			} else if (test) {
			// Warn about change from old behavior when specifying a testset.
			// We could make the -r default 1 above if there's a testset, but:
			//  this would be an arbitrary change in defaults for someone new to wndchrm - would they expect this change?
			//  for the sake of maintaining consistency for someone familiar with wndchrm's previous undocumented behavior.
			//  leaving us wondering:  Which is worse?
			// For now, the compromise is to issue a warning, while keeping the defaults as documented.
			// Additionally, previously all the training samples were used balanced or not, and now its balanced except with -r#
			// This is also inconsistent with previous behavior.
// 				if (testset && (split_ratio < 1.0 || balanced_splits)) {
// 					catError (
// 						"WARNING: Change from previous versions when specifying a testset with the test command.\n"
// 						"  Previously the entire dataset was used for training (balanced or not) if a testset was specified,\n"
// 						"  making it impossible to test random training splits on a testset. This was undocumented behavior.\n"
// 						"  If the old behavior is desired, use the classify command, or specify -r1 (or -r#1).\n"
// 					);
// 				}
			}

			for (ignore_group=0;ignore_group<=assess_features;ignore_group++) {
				split_and_test(dataset, report_file, argc, argv, MAX_CLASS_NUM, method, &featureset, split_ratio, balanced_splits, max_features, used_mrmr,splits_num,report,max_training_images,
					exact_training_images,max_test_images,phylib_path,distance_method,phylip_algorithm,export_tsv,first_n,weight_file_buffer,weight_vector_action,N,
					testset,ignore_group,tile_areas,max_tile,image_similarities, random_splits);
			}
	
			// report any warnings
			showError (0,NULL);
       } // test or classify.
 
 
       
     } // no params left for dataset / test set.
     else ShowHelp();

     return(1);
}



