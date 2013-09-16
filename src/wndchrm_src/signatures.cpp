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

#include <cmath>
#include <cfloat> // Has definition of DBL_EPSILON, FLT_EPSILON
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h> // for locking stuff
#include <errno.h>
#include <time.h>
#include <unistd.h> // apparently, for close() only?
#define OUR_EPSILON FLT_EPSILON*6
#define FLOAT_EQ(x,v) (((v - FLT_EPSILON) < x) && (x <( v + FLT_EPSILON)))
#define OUR_EQ(x,v) (((v - OUR_EPSILON) < x) && (x <( v + OUR_EPSILON)))
#include "cmatrix.h"
#include "TrainingSet.h"
#include "colors/FuzzyCalc.h"
#include "signatures.h"
#include "WORMfile.h"
#include "FeatureNames.h"
#include "FeatureAlgorithms.h"


/* global variable */
extern int verbosity;

// static signatures::max_sigs
long signatures::max_sigs = NUM_DEF_FEATURES;

//---------------------------------------------------------------------------
/*  signatures (constructor)
*/
signatures::signatures() {
	data.clear();
	version = 0;
	feature_vec_type = StdFeatureComputationPlans::fv_unknown;
	count=0;
	allocated = 0;
	sample_class=0;
	full_path[0]='\0';
	sample_name[0]='\0';
	NamesTrainingSet=NULL;   
	ScoresTrainingSet=NULL;
	wf = NULL;
}
//---------------------------------------------------------------------------

signatures::~signatures() {
	if (wf) delete wf;
	wf = NULL;
}
/* duplicate
*/
signatures *signatures::duplicate() {
	signatures *new_samp;
	new_samp=new signatures();
	new_samp->sample_class=sample_class;
	new_samp->sample_value=sample_value;   
	new_samp->interpolated_value=interpolated_value;
	new_samp->count=count;
	new_samp->NamesTrainingSet=NamesTrainingSet;   
	new_samp->ScoresTrainingSet=ScoresTrainingSet;
	strcpy(new_samp->full_path,full_path);

	new_samp->Resize (count);
	new_samp->data = data;
	wf = NULL;
	new_samp->version = version;
	new_samp->feature_vec_type = feature_vec_type;
	return(new_samp);
}

/* Resize
   Allocate memory for specified number of signatures
   nsigs -size_t - number of signatures to preallocate
*/
void signatures::Resize(size_t nsigs) {
	allocated = nsigs;
	data.resize (nsigs);
}

/* Add
   add a signature
   name -char *- the name of the signature (e.g. Multiscale Histogram bin 3)
   value -double- the value to add
*/
void signatures::Add(const char *name,double value) {
	if (name && *name && NamesTrainingSet) {
		char *char_p = ((TrainingSet *)(NamesTrainingSet))->SignatureNames[count];
		if (! *char_p) strcpy(char_p,name);
	}
	
	if (count == 0 && allocated == 0) Resize (max_sigs);
	else if (count >= allocated) Resize (count + 1024);
	data[count]=value;
	count++;
	if (count > max_sigs) max_sigs = count;
}


void signatures::SetFeatureVectorType () {
	if (feature_vec_type == StdFeatureComputationPlans::fv_unknown) {
		switch (count) {
			case NUM_LC_FEATURES:
				feature_vec_type = StdFeatureComputationPlans::fv_long_color;
			break;
			case NUM_L_FEATURES:
				feature_vec_type = StdFeatureComputationPlans::fv_long;
			break;
			case NUM_C_FEATURES:
				feature_vec_type = StdFeatureComputationPlans::fv_short_color;
			break;
			case NUM_DEF_FEATURES:
				feature_vec_type = StdFeatureComputationPlans::fv_short;
			break;
			default:
			break;
		}
	}
}

/* Clear
   clear all signature values
*/
void signatures::Clear() {
	data.clear();
	allocated = 0;
	count = 0;
	feature_vec_type = StdFeatureComputationPlans::fv_unknown;
}

int signatures::IsNeeded(long start_index, long group_length)
{  int sig_index;
   if (!ScoresTrainingSet) return(1);
   for (sig_index=start_index;sig_index<start_index+group_length;sig_index++)
     if (((TrainingSet *)(ScoresTrainingSet))->SignatureWeights[sig_index]>0) return(1);
   return(0);
}

void signatures::compute_plan (const ImageMatrix &matrix, const FeatureComputationPlan *plan) {
	FeatureComputationPlanExecutor executor (plan);
	
	version = CURRENT_FEATURE_VERSION;
	feature_vec_type = plan->feature_vec_type;
	
	Resize (plan->n_features);
	executor.run(&matrix, data, 0);
	
	// update the feature count and the max_count;
	count = plan->n_features;
	if (count > max_sigs) max_sigs = count;
	
	// If we have an attached NamesTrainingSet, copy the feature names over, but only the first time.
	if (NamesTrainingSet) {
		for (int i = 0; i < count; i++) {
			char *char_p = ((TrainingSet *)(NamesTrainingSet))->SignatureNames[i];
			if (! *char_p) {
				strcpy(char_p,plan->getFeatureNameByIndex(i).c_str());
			}
		}
	}
}


/* normalize
   normalize the signature values using the maximum and minimum values of the training set
   ts -TrainingSet *- the training set according which the signature values should be normalized
*/
void signatures::normalize(void *TrainSet)
{
	int sig_index;
	TrainingSet *ts;
	ts=(TrainingSet *)TrainSet;
	double sig_val, sig_min, sig_max;
	for( sig_index = 0; sig_index < count; sig_index++ ) {
		sig_val = data[ sig_index ];
		sig_min = ts->SignatureMins[ sig_index ];
		sig_max = ts->SignatureMaxes[ sig_index ];

		if (std::isnan(sig_val) || sig_val < sig_min || (sig_max - sig_min) < DBL_EPSILON) sig_val = 0; 
		else if( sig_val > sig_max )
			sig_val = 100;
		else
			sig_val = 100 * ( (sig_val - sig_min) / (sig_max - sig_min) );

		data[ sig_index ] = sig_val;
	}
}


/* FileClose
   Closes a value file.  This is the closing command for files opened with ReadFromFile.
   This closes the stream as well as filedescriptor
*/
void signatures::FileClose()
{
	if (wf) {
		wf->finish();
	}
}

int signatures::SaveToFile (int save_feature_names) {
	int sig_index;
	char buffer[IMAGE_PATH_LENGTH+SAMPLE_NAME_LENGTH+1];

	if (!wf) {
		if (strlen (full_path) > 0)
			wf = new WORMfile (GetFileName (buffer));
		else {
			fprintf(stderr, "Cannot write to .sig file - full_path not set.\n");
			return(0);
		}
	}

	if (!wf || !(wf->status == WORMfile::WORM_WR) ) {
		printf("Cannot write to .sig file: cannot open sigfile for writing\n");
		return(0);
	}
	FILE *wf_fp = wf->fp();

	if ( NamesTrainingSet && ((TrainingSet *)(NamesTrainingSet))->is_continuous ) {
		fprintf(wf_fp,"%f\t%d.%d\n",sample_value,version,feature_vec_type);  /* save the continouos value */
	} else {
		fprintf(wf_fp,"%d\t%d.%d\n",sample_class,version,feature_vec_type);  /* save the class index */
	}
	fprintf(wf_fp,"%s\n",full_path);
	for (sig_index=0; sig_index < count; sig_index++) {
		if (save_feature_names && NamesTrainingSet)
			fprintf(wf_fp,"%f\t%s\n",data[sig_index],((TrainingSet *)NamesTrainingSet)->SignatureNames[sig_index]);
		else
			fprintf(wf_fp,"%f\n",data[sig_index]);
	}
   return(1);
}


int signatures::LoadFromFile(char *filename) {
	char buffer[IMAGE_PATH_LENGTH+SAMPLE_NAME_LENGTH+1];
	WORMfile *wf_temp = NULL;
	int ret = 0;

	if (!filename || *filename == '\0')
		GetFileName (buffer);
	else strncpy (buffer,filename,sizeof(buffer));

	wf_temp = new WORMfile (buffer, true); // readonly
	if (wf_temp->status == WORMfile::WORM_RD) {
		LoadFromFilep (wf_temp->fp());
		ret = 1;
	}
	delete wf_temp;  // closes readonly, unlinks write-locked.
	return (ret);
}

void signatures::LoadFromFilep (FILE *value_file) {
	char buffer[IMAGE_PATH_LENGTH+SAMPLE_NAME_LENGTH+1],*p_buffer, name[SIGNATURE_NAME_LENGTH];
	int version_maj = 0, version_min = 0;
	double val;

	/* read the class or value and version */
	fgets(buffer,sizeof(buffer),value_file);
	if (NamesTrainingSet && ((TrainingSet *)(NamesTrainingSet))->is_continuous) {
		sscanf (buffer, "%lf%*[\t ]%d.%d", &sample_value, &version_maj, &version_min);
		sample_class = 1;
	} else {
		sscanf (buffer, "%hu%*[\t ]%d.%d", &sample_class, &version_maj, &version_min);
	}
	// If we did not read a version, then it is 1.0
	if (version_maj == 0) {
		version = 1;
		feature_vec_type = StdFeatureComputationPlans::fv_unknown;
	} else {
		version = version_maj;
		feature_vec_type = version_min;
	}
	/* read the path */
	fgets(buffer,sizeof(buffer),value_file);
	chomp (buffer);
	strcpy(full_path,buffer);
	
	/* read the feature values */
	p_buffer=fgets(buffer,sizeof(buffer),value_file);
	chomp (p_buffer);
	while (p_buffer) {
		name[0] = '\0';
		sscanf (buffer, "%lf%*[\t ]%[^\t\r\n]", &val, name);
		Add(name,val);
		p_buffer=fgets(buffer,sizeof(buffer),value_file);
		chomp (p_buffer);
	}

	// FIXME: There is opportunity here to check for inconsistent number of features if minor version is specified.
	SetFeatureVectorType();
}



/*
  Yet another variant of reading from a file.
  The filename is computed from full_path and sample_name using GetFileName
  If the file is successfully opened and write-locked, return 0 (wf.status = WORMfile::WORM_WR).
  If another process has a lock, return 0 (wf.status = WORMfile::WORM_BUSY).
  If the file exists, and is not locked, the sigs will be loaded from it, and no lock will be issued. (return 1, (wf.status = WORMfile::WORM_FINISHED))
  If an error occurs in obtaining the lock (if necessary) or creating the file (if necessary) or reading it (if possible), return -1.
*/
int signatures::ReadFromFile (bool wait) {
	char buffer[IMAGE_PATH_LENGTH+SAMPLE_NAME_LENGTH+1];

	if (!wf) wf = new WORMfile (GetFileName (buffer), wait, wait);
	else wf->reopen(wait, wait);

	if (!wf) return (-1);
	if (wf->status == WORMfile::WORM_BUSY) {
		return (0);
	} else if (wf->status == WORMfile::WORM_WR) {
		return (0);
	} else if (wf->status == WORMfile::WORM_RD) {
		Clear(); // reset sample count
		LoadFromFilep (wf->fp());
		wf->finish(); // this unlocks, closes, etc.
		// Of course, if it was empty afterall, its an error.
		if (count < 1) {
			return (NO_SIGS_IN_FILE);
		} else {
			return (1);
		}
	} else {
	// I/O error
		return (-1);
	}
}

/*
  get the filename for storing the signature.
  The filename is generated from the full path of the image, plus a sample name.
  The sample name (i.e. _0_0 for tile 0,0) is set externally depending on what sample of the image (referred to by full_path) the signature pertains to.
  It is stored internally so that a signature object can get to its own file without additional information.
*/
char *signatures::GetFileName (char *buffer) {
	char *char_p;

	if (wf) {
		strcpy (buffer, wf->path.c_str());
		return buffer;
	}

	strcpy(buffer,full_path);
	char_p = strrchr(buffer,'.');
	if (!char_p) char_p=buffer+strlen(buffer);

	sprintf(char_p,"%s.sig",sample_name);
	return (buffer);
}


// Based on
// Usable AlmostEqual function
// By Bruce Dawson
// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
// returns the ulps (floating point representations) b/w the two inputs
// uses a union for punning to avoid strict aliasing rules
/*
int diffUlps(float A, float B)
{
	union fi_union {
	int32_t i;
	float f;
	};
	fi_union fiA,fiB;
	fiA.f = A;
	fiB.f = B;

    int32_t aInt = fiA.i;
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;
    // Make bInt lexicographically ordered as a twos-complement int
    int32_t bInt = fiB.i;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;
    int intDiff = abs(aInt - bInt);
    return (intDiff);
}
*/

/*
  This function is used to determine (somewhat quickly) if the features stored in a file match those
  that would be calculated from the passed-in matrix.
  A partial signature calculation is done to determine the match.  An exact match (within FLT_EPSILON) of every feature is required to return 1.
  If the file can't be opened, or if the match is inexact, 0 is returned.
*/
int signatures::CompareToFile (const ImageMatrix &matrix, char *filename, int compute_colors, int large_set) {
// 	signatures file_sigs;
// 	double vec[72];
// 	int i,file_index;
// 
// 	if (! file_sigs.LoadFromFile (filename) ) {
// 		// Any errors are ignored
// 		errno = 0;
// 		return (0);
// 	}
// 	if (verbosity > 1) printf ("compare %s to computed \n",filename);
// 
// 	// 20 features long: 323-342, standard: N/A
// 	if (large_set) {
// 		matrix.fractal2D(20,vec);
// 		file_index = 323;
// // for (i = 0; i< 20; i++) printf ("fractal2D computed %15.10f\tfrom file: %15.10f\tdiff: %f\tulps: %d\n",vec[i],file_sigs.data[file_index+i],
// // (file_sigs.data[file_index+i] - vec[i])/FLT_EPSILON
// // ,diffUlps(file_sigs.data[file_index+i],vec[i])
// // );
// 		for (i = 0; i< 20; i++) if (!OUR_EQ(file_sigs.data[file_index+i],vec[i])) {
// 			if (verbosity > 1) printf ("fractal2D mismatch computed %15.10f\tfrom file: %15.10f\n",vec[i],file_sigs.data[file_index+i]);
// 			return (0);
// 		}
// 	}
// 	if (verbosity > 1) printf ("fractal2D match\n");
// 	
// 	// 28 features long: 253-280, standard: 485-512
// 	matrix.HaralickTexture2D(0,vec);
// 	if (large_set) file_index = 253;
// 	else file_index = 485;
// 	for (i = 0; i < 28; i++) if (!OUR_EQ(file_sigs.data[file_index+i],vec[i])) {
// 		if (verbosity > 1) printf ("HaralickTexture2D mismatch computed %15.10f\tfrom file: %15.10f\n",vec[i],file_sigs.data[file_index+i]);
// 		return (0);
// 	}
// 	if (verbosity > 1) printf ("HaralickTexture2D match\n");
// 	
// 	// 72 features long: 133-204, standard: 881-952
// // 	long output_size;   /* output size is normally 72 */
// // 	matrix->zernike2D(vec,&output_size);
// // 	if (large_set) file_index = 133;
// // 	else file_index = 881;
// // 	for (i = 0; i < 72; i++) if (!OUR_EQ(file_sigs.data[file_index+i],vec[i])) {
// // 		if (verbosity > 1) printf ("zernike2D mismatch computed %15.10f\tfrom file: %15.10f\n",vec[i],file_sigs.data[file_index+i]);
// // 		return (0);
// // 	}
// // 	if (verbosity > 1) printf ("zernike2D match.\n");
// 
// 	if (verbosity > 1) printf ("Match found.\n");
// 	return (1);


return (0);
}


