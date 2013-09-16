/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2007 Open Microscopy Environment                             */
/*         Massachusetts Institue of Technology,                                 */
/*         National Institutes of Health,                                        */
/*         University of Dundee                                                  */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
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


#ifndef TrainingSetH
#define TrainingSetH
//---------------------------------------------------------------------------

// STL stuff
#include <vector>
#include <list>
#include <string>
#include "signatures.h"
#include "config.h" // for version info

#define MAX_CLASS_NUM 1024
#define MAX_CLASS_NAME_LENGTH 50
#define MAX_SAMPLES_PER_IMAGE 4096

#define WNN 0
#define WND 1

// N.B.: There is almost certainly code to fix if this is other than 1
#define CONTINUOUS_CLASS_LABEL ""
#define CONTINUOUS_CLASS_INDEX 1

#define UNKNOWN_CLASS_LABEL ""
// N.B.: There is almost certainly code to fix if this is other than 0
#define UNKNOWN_CLASS_INDEX 0

#define CANT_OPEN_DIRECTORY                -1
#define CANT_OPEN_FIT                      -2
#define CANT_LOAD_ALL_SIGS                 -3
#define ADDING_CLASS_TO_CONTINUOUS_DATASET -4
#define CANT_ADD_UNORDERED_CLASS           -5
#define TOO_MANY_CLASSES                   -6
#define CONTINUOUS_DATASET_WITH_CLASSES    -7
#define ADDING_SAMPLE_TO_UNDEFINED_CLASS   -8
#define INCONSISTENT_FEATURE_COUNT         -9
#define INCONSISTENT_FEATURE_VECTORS       -10

typedef struct {
	char bounding_rect_base[16];
	rect bounding_rect;
	char downsample_base[16];
	int downsample;
	char normalize_base[16];
	int mean;
	int stddev;
} preproc_opts_t;

typedef struct {
	char rot_base[16]; // CLI option+params
	int rotations;
	char tile_base[16]; // CLI option+params
	int tiles_x;
	int tiles_y;
} sampling_opts_t;

typedef struct {
	char compute_colors_base[16]; // CLI option+params
	int compute_colors;
	char large_set_base[16]; // CLI option+params
	int large_set;
} feature_opts_t;

typedef struct {
	preproc_opts_t preproc_opts;
	sampling_opts_t sampling_opts;
	feature_opts_t feature_opts;
	int n_features;
	int n_samples;
	struct {
		char sample_name[SAMPLE_NAME_LENGTH];
		int rot_index;
		int tile_index_x;
		int tile_index_y;
	} samples[MAX_SAMPLES_PER_IMAGE];	
} featureset_t;


// Set up our struct for keeping track of per-feature-group statistics.
class FeatureGroup;
class FeatureInfo;
typedef struct {
	std::string name;
	const FeatureGroup *featuregroup_info;
	double min;
	double max;
	double sum_weight;
	double sum_weight2;
	double mean;
	double stddev;
	int n_features;
} featuregroup_stats_t;

typedef struct sort_by_mean_weight_t {
	bool operator() (featuregroup_stats_t i,featuregroup_stats_t j) { return (j.mean < i.mean);}
} sort_by_mean_weight_func;

typedef std::vector<featuregroup_stats_t> featuregroups_t;

// Set up our struct for keeping track of per-feature stuff.
typedef struct {
	std::string name; // N.B.: name as read from file
	const FeatureInfo *feature_info;
	double weight;
} feature_stats_t;

typedef struct sort_by_weight_t {
	bool operator() (feature_stats_t i,feature_stats_t j) { return (j.weight < i.weight);}
} sort_by_weight_func;

typedef struct {
	signatures *sig;
	int rot_index;
	int tile_index_x;
	int tile_index_y;
	bool saved;
	bool added;
} feature_vec_info_t;

typedef std::vector<feature_stats_t> features_t;


typedef struct {
	double accuracy;
	double avg_class_accuracies;
	long known_images;
	long accurate_predictions;
	double plus_minus;
	double classification_p_value;

   double *tile_area_accuracy;            /* used for the different accuracies of the different tile areas     */
   unsigned short *confusion_matrix;      
   unsigned short *training_images;       /* array to store number of training images per class */
   unsigned short *testing_images;        /* array to store number of testing images per class */
   double *class_accuracies;              /* array to store accuracies per class */
   double *similarity_matrix;             /* matrix - used for the similarities between the classes            */
   double *class_probability_matrix;      /* matrix - average class probabilities (un-normalized similarities) */
   double *image_similarities;            /* matrix - used for the similarity values between all test images   */
   features_t feature_stats;
   featuregroups_t featuregroups_stats;
   double feature_weight_distance;
   char *individual_images;                /* a string of the individual image predictions. used for the report */
   unsigned short method; 
   double pearson_coefficient;             /* pearson correlation between the predicted and actual value        */
   double avg_abs_dif;                     /* average absolute difference between the actual and the predicted values */
   double pearson_p_value;
	std::vector< std::list< float > > marginal_probabilities; //keep track of marginal probs
} data_split;

class TrainingSet
{
public:
/* properties */
	char name[256];                       /* Name of dataset - source_path from last '/' to last '.'    */
	char source_path[256];                       /* Path we read this set from     */
	int feature_vec_version, feature_vec_type; // from the signatures class
   signatures **samples;                                           /* samples data                              */
   char SignatureNames[MAX_SIGNATURE_NUM][SIGNATURE_NAME_LENGTH];  /* names of the signatures (e.g. "MultiScale Histogram bin 3) */
   double SignatureWeights[MAX_SIGNATURE_NUM];                     /* weights of the samples                    */
   double SignatureMins[MAX_SIGNATURE_NUM];                        /* minimum value of each signature           */
   double SignatureMaxes[MAX_SIGNATURE_NUM];                       /* maximum value of each signature           */
	  featuregroups_t *aggregated_feature_stats;                    // only the top level TrainingSet has this allocated by wndcharm.cpp::split_and_test
   long class_num;                                                 /* number of known/defined classes (may be 0 if all samples are unknown, may be 1 when is_continuous, or for 1 known discrete class */
   char **class_labels;                                            /* labels of the classes                     */
   long *class_nsamples;                                           /* sample counts in each class               */
   int *train_class;                                               /* class indexes into training set           */
   int  is_continuous;                                             /* A numeric/continuous dataset.  sample_class = 0 or 1 for all samples. class_nsamples is valid, but class_labels is not. */
   int  is_numeric;                                                /* All class labels can be interpreted as numeric values (is_continuous can be false when is_numeric is true) */
   int  is_pure_numeric;                                           /* All class labels are numerical (no characters other than those than can be part of a valid double - INF, NAN, etc are technically valid, but not in our case)  */
   long count;                                                     /* the number of samples in the training set */
   long signature_count;                                           /* the number of signatures (< MAX_SIGNATURE_NUM) */
   long color_features;                                            /* color signatures are used                 */
/* methods */
   TrainingSet(long samples_num, long class_num);                  /* constructor                               */
   ~TrainingSet();                                                 /* destructor                                */
   int AddAllSignatures();                                         /* load the sample feature values from corresponding files */
	int AddImageFile(char *filename, unsigned short sample_class, double sample_value, int save_sigs, featureset_t *featureset, int skip_sig_comparison_check = 0);
	int LoadFromFilesDir(char *path, unsigned short sample_class, double sample_value, int save_sigs, featureset_t *featureset, int skip_sig_comparison_check = 0);
	int LoadFromPath(char *path, int save_sigs, featureset_t *featureset, int make_continuous, int skip_sig_comparison_check = 0);
   double ClassifyImage(TrainingSet *TestSet, int test_sample_index,int method, int tiles, int tile_areas, TrainingSet *TilesTrainingSets[], int max_tile,int rank, data_split *split, double *similarities);  /* classify one or more images */
   double Test(TrainingSet *TestSet, int method, int tiles, int tile_areas, TrainingSet *TilesTrainingSets[], int max_tile,long rank, data_split *split);     /* test      */
   int SaveToFile(char *filename);                                 /* save the training set values to a file    */
	bool IsFitFile(char *filename);                                /* checks if its a proper fit file by making sure the first three lines are pure numeric */
   int ReadFromFile(char *filename);                               /* read the training set values from a file  */
   int SaveWeightVector(char *filename);                           /* save the weights of the features into a file */
   double LoadWeightVector(char *filename, double factor);         /* load the weights of the features from a file and assign them to the features of the training set */
   void SetAttrib(TrainingSet *set);                               /* copy the attributes from one training set to another */   
   int split(int randomize,double ratio,TrainingSet *TrainSet,TrainingSet *TestSet, unsigned short tiles, int train_samples, int test_samples, data_split *split); /* random split to train and test */
   int SplitAreas(long tiles_num, TrainingSet **TrainingSets);    /* split a tiled dataset into several datasets such that each dataset is one tile location */
   void RemoveClass(long class_index);                             /* remove a class                            */
	void MakeContinuous(char *label);                              /* Make an existing dataset with defined classes continuous - label is the "units" */
	void MarkUnknown(long class_index);                            /* mark the given class as unknown (move + reorder class labels, reassign sample classes to class 0 */
   int AddClass(char *label);                                      /* add a discrete class    */
   int AddContinuousClass (char *label);                           /* add a continuous class - not that only one can be added */
   int AddSample(signatures *new_sample);                          /* add signatures computed from one image    */
   void normalize();                                               /* normalize the values of the signatures to [0,100] */
   void SetmRMRScores(double used_signatures,double used_mrmr);                     /* set mRMR scores to the features           */
   void SetFisherScores(double used_signatures, double used_mrmr, data_split *split);/* compute the fisher scores for the signatures  */
   int IgnoreFeatureGroup(long index,char *group_name);            /* set the Fisher Score of a group of image features to zero */
   double distance(signatures *sample1, signatures *sample2,double power);  /* Find the weighted Euclidean distance between two samples  */
   long WNNclassify(signatures *test_sample, double *probabilities, double *normalization_factor, signatures **closest_sample);/* classify a sample using weighted nearest neighbor */
   long classify2(char* name, int test_sample_index, signatures *test_sample, double *probabilities,double *normalization_factor); /* classify using -5                         */
   double InterpolateValue(signatures *test_sample, int method, int N, signatures **closest_sample, double *closest_dist);  /* interpolate a value */
   long classify3(signatures *test_sample, double *probabilities,double *normalization_factor);
   double pearson(int tiles,double *avg_abs_dif,double *p_value);                  /* a pearson correlation of the interpolated and the class labels (if all labels are numeric) */
   long PrintConfusion(FILE *output_file, unsigned short *confusion_matrix, double *similarity_matrix);//, unsigned short dend_file, unsigned short method);  /* print a confusion or similarity matrix */
   long dendrogram(FILE *output_file, char *data_set_name, char *phylib_path, int nodes_num,double *similarity_matrix, char **labels,unsigned short sim_method,unsigned short phylip_algorithm);  /* create a dendrogram */
	 long report(FILE *output_file, int argc, char **argv, char *output_file_name, std::vector<data_split> &splits, unsigned short split_num, featureset_t *featureset, int max_train_images,char *phylib_path, int distance_method, int phylip_algorithm, int export_tsv, TrainingSet *testset, int image_similarities, bool use_err_bars);
   void Summarize(featureset_t *featureset);  /* class summary to stdout */
};

int check_numeric (char *s, double *samp_val);
void chomp (char *line);


#endif
