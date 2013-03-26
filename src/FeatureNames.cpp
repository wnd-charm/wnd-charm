/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2007 Open Microscopy Environment                                */
/*       Massachusetts Institue of Technology,                                   */
/*       National Institutes of Health,                                          */
/*       University of Dundee                                                    */
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
/* Written by:                                                                   */
/*      Ilya G. Goldberg <goldbergil [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

// This defines a big C string constant (oldFeatureNamesFileStr) with tab-delimited old_feature new_feature lines.
// Its defined this way because doing map declarations on the stack takes forever to compile,
// blows up memory during compilation with optimization, and results in a monstrously huge object file 5x bigger than the rest of the library
#include "OldFeatureNamesFileStr.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "FeatureAlgorithms.h"
#include "ImageTransforms.h"
#include "FeatureNames.h"

/* global variable */
extern int verbosity;

// hopeful pre-main initialization of the featurename lookup map
const bool FeatureNames::initialized() {
	static cnm_t  &channels = channels_map ();
	static tnm_t  &transforms = transforms_map ();
	static fam_t  &feature_algorithms = feature_algorithms_map ();
	static fgnm_t &feature_groups = feature_groups_map ();
	static fnm_t  &feature_names = feature_names_map ();
	static ofnm_t &old_features = old_features_map ();

	return ( ! (
		channels.empty() ||
		transforms.empty() ||
		feature_algorithms.empty() ||
		feature_groups.empty() ||
		feature_names.empty() ||
		old_features.empty()
	) );
}
const bool FeatureNames::init_maps_ = FeatureNames::initialized ();


// testing only
// #define CYCLES 100
// 
// typedef struct {
// 	std::string name;
// 	const FeatureNames::FeatureInfo *feature_info;
// 	int n_hits;
// } featuregroup_stats_t;
// struct sort_by_n_features_t {
// 	bool operator() (featuregroup_stats_t i,featuregroup_stats_t j) { return (i.feature_info->group->algorithm->n_features<j.feature_info->group->algorithm->n_features);}
// } sort_by_n_features;
// // testing only
// 
// int main () {
// 	int i,j,n_featnames;
// 	FeatureNames::ofnm_t::const_iterator ofnm_it;
// 	std::vector<std::string> featurenames;
//	static ofnm_t &old_features = old_features_map();
// 
// 	featurenames.reserve(FeatureNames::old_features.size());
// 	for(ofnm_it = FeatureNames::old_features.begin(); ofnm_it != FeatureNames::old_features.end(); ++ofnm_it ) {
// 		featurenames.push_back( ofnm_it->first );
// 	}
// 	n_featnames = featurenames.size();
// 
// 	std::set<std::string> fgs;
// 	std::set<std::string>::iterator fgs_it;
// 	std::vector<featuregroup_stats_t> fgsv;
// 	std::vector<featuregroup_stats_t>::iterator fgsv_it;
// 	featuregroup_stats_t featuregroup_stats;
// 	FeatureNames::FeatureInfo const *featureinfo;
// 	for (i = 0; i< n_featnames; i++) {
// 		featureinfo = FeatureNames::getFeatureInfoByName ( featurenames[i].c_str() );
// 		fgs_it = fgs.find(featureinfo->group->name);
// 		if (fgs_it == fgs.end()) {
// 			featuregroup_stats.name = featureinfo->group->name;
// 			featuregroup_stats.feature_info = featureinfo;
// 			featuregroup_stats.n_hits = 0;
// 			fgs.insert(featureinfo->group->name);
// 			fgsv.push_back (featuregroup_stats);
// 		}
// 	}
// 
// 	sort (fgsv.begin(), fgsv.end(), sort_by_n_features);
// 	for(fgsv_it = fgsv.begin(); fgsv_it != fgsv.end(); ++fgsv_it ) {
// 		printf ("%s: %d\n",fgsv_it->name.c_str(), fgsv_it->feature_info->group->algorithm->n_features);
// 	}
// 
// 
// 
// // 	fgnm_t *fgnm = featureGroups();
// // 	fgnm_it_t fgnm_it;
// // 	for(fgnm_it = (*fgnm).begin(); fgnm_it != (*fgnm).end(); ++fgnm_it ) {
// // 		printf ("%s: %d\n",fgnm_it->second.name.c_str(), fgnm_it->second.n_features);
// // 	}
// 
// 
// 
// 	int found=0;
// 	int missed=0;
// 	double cpu_secs;
// 	clock_t start, end;
// 
// // time parse feature - found
// 	found=0;
// 	missed=0;
// 	start = clock();
// 	for (j = 0; j < CYCLES; j++) {
// 		for (i = 0; i< n_featnames; i++) {
// 			if (! FeatureNames::getFeatureInfoByName ( featurenames[i].c_str() ) ) missed++;
// 			else found++;
// 		}
// 	}
// 	end = clock();
// 	cpu_secs = ((double) (end - start)) / CLOCKS_PER_SEC;
// 	printf ("       +parse %9d lookups %9d misses %9.2g secs, %9.0f lookups/sec\n",found+missed,missed,cpu_secs,(double)(found+missed)/cpu_secs);
// 
// 
// // time parse feature - nonexistant
// 	found=0;
// 	missed=0;
// 	start = clock();
// 	for (j = 0; j < CYCLES; j++) {
// 		for (i = 0; i< n_featnames; i++) {
// 		// worst- case can't find algorithm
// 		//	if ( !FeatureNames::getFeatureInfoByName ( "CoOcMat_MaximalCorrelationCoefficientDif_ChebyshevFFT MaxCorrCoef () [123]" ) ) missed++;
// 		// second-worst - malformed invalid algorithm: no '()'
// 			if ( !FeatureNames::getFeatureInfoByName ( "CoOcMat_MaximalCorrelationCoefficientDif_ChebyshevFFT MaxCorrCoef [123]" ) ) missed++;
// 		// best-case invalid name - no []
// 		//	if ( !FeatureNames::getFeatureInfoByName ( "CoOcMat_MaximalCorrelationCoefficientDif_ChebyshevFFT MaxCorrCoef" ) ) missed++;
// 			else found++;
// 		}
// 	}
// 	end = clock();
// 	cpu_secs = ((double) (end - start)) / CLOCKS_PER_SEC;
// 	printf ("       -parse %9d lookups %9d misses %9.2g secs, %9.0f lookups/sec\n",found+missed,missed,cpu_secs,(double)(found+missed)/cpu_secs);
// 
// }

/*
New-style feature names follow this style:
  Zernike Coefficients (Wavelet (Edge ())) [21]
  Zernike Coefficients (Wavelet ()) [21]
  Zernike Coefficients () [21]
Algorithm name followed by a set of nested transforms in parentheses, then feature index within the group in square brackets.
The inner-most parentheses contain an optional channel label
White-space is not used as a delimiter - only '(',')','[' and ']'. Leading and trailing whitespace is eliminated from transforms and algorithm names.
The parentheses are *required* to parse transforms. In their absence, the entire feature name (other than square brackets) is taken as the algorithm name.
The square brackets are required for indexes.  In their absence, an index of 0 is assumed.
The index is 0-based.
The transforms are applied in reverse order to what is listed, as expected for the nested representation.
In the example above, Edge transform first, then Wavelet, then Zernike coefficients on the Wavelet.
The feature and group names reported in .name fields are normalized for whitespace as in the example above.
*/

const FeatureInfo *FeatureNames::getFeatureInfoByName (const std::string &featurename_in) {
	static fnm_t &feature_names = FeatureNames::feature_names_map();

	fnm_t::const_iterator fnm_it = feature_names.find(featurename_in);
	if (fnm_it != feature_names.end()) return (fnm_it->second);

	std::string featurename;
	const std::string &featurename_old = oldFeatureNameLookup (featurename_in);
	if (!featurename_old.empty()) featurename = featurename_old;
	else featurename = featurename_in;

// parse out the group index
	int index = -1;
	size_t found_left = featurename.find_last_of ('[');
	size_t found_right = featurename.find_last_of (']');
	if (found_left == std::string::npos || found_right == std::string::npos || found_right-found_left-1 < 1)
		index = -1;
	else
		index = atoi ( (featurename.substr (found_left+1,found_right-found_left-1)).c_str() );
// parse out the group name
	size_t found;
	std::string groupname;
	groupname = featurename.substr (0,found_left);
	// clean trailing whitespace
	found=groupname.find_last_not_of(" \t\f\v\n\r");
	if (found != std::string::npos)
		groupname.erase(found+1);

	const FeatureGroup *featuregroup = getGroupByName (groupname);

	// For now, if we got an invalid (unknown) index, assume that its 0.
	// note that the feature name is normalized for whitespace
	if (index < 0) index = 0;
	else if (index < featuregroup->algorithm->n_features) {
		featurename = featuregroup->labels[index];
	} else {
		char buf[64];
		sprintf(buf," [%d]",index);
		featurename = featuregroup->name + buf;
	}

	// For unknown featuregroups, make sure that n_features is always big enough to accomodate the index (if any)
	if (index > featuregroup->algorithm->n_features) {
	// bypass static constraint
		int *index_p = (int *)&(featuregroup->algorithm->n_features);
		*index_p = index;
	}

	return (feature_names[featurename_in] = new FeatureInfo (featurename, featuregroup, index));
}

// This returns an iterator to the algorithm map by string lookup
const FeatureAlgorithm *FeatureNames::getFeatureAlgorithmByName (const std::string &name) {
	static fam_t  &feature_algorithms = feature_algorithms_map ();

	fam_t::const_iterator fam_it = feature_algorithms.find(name);

	if (fam_it == feature_algorithms.end()) {
		return (NULL);
	} else {
		return (fam_it->second);
	}
}

bool FeatureNames::registerFeatureAlgorithm (const FeatureAlgorithm *algorithm) {
	static fam_t  &feature_algorithms = feature_algorithms_map ();
	fam_t::const_iterator fam_it = feature_algorithms.find(algorithm->name);

	if (fam_it != feature_algorithms.end())
		return true;

	feature_algorithms[algorithm->name] = algorithm;
	return true;
}

// This should return a channel object by string lookup
const Channel *FeatureNames::getChannelByName (const std::string &name) {
	static cnm_t  &channels = channels_map ();

	cnm_t::const_iterator cnm_it = channels.find(name);

	if (cnm_it == channels.end()) {
		return (channels[name] = new Channel (name));
	} else {
		return (cnm_it->second);
	}
}

// This should return a transform object by string lookup
const ImageTransform *FeatureNames::getTransformByName (const std::string &name) {
	static tnm_t &image_transforms = transforms_map ();
	tnm_t::const_iterator tnm_it = image_transforms.find(name);

	if (tnm_it == image_transforms.end()) {
		return (image_transforms[name] = new EmptyTransform (name));
	} else {
		return (tnm_it->second);
	}
}

bool FeatureNames::registerImageTransform (const ImageTransform *algorithm) {
	static tnm_t &image_transforms = transforms_map ();
	tnm_t::const_iterator tnm_it = image_transforms.find(algorithm->name);

	if (tnm_it != image_transforms.end())
		return true;

	image_transforms[algorithm->name] = algorithm;
	return true;
}

FeatureGroup::FeatureGroup (const std::string &s, const FeatureAlgorithm *f, const Channel *c, std::vector<ImageTransform const *> &t) {
	name = s; algorithm = f; channel = c; transforms = t;
	if (algorithm)
		for (int i = 0; i < algorithm->n_features; i++) {
			char buf[16];
			sprintf (buf, " [%d]",i);
			labels.insert (labels.end(), name + buf);
		}
}


const FeatureGroup *FeatureNames::getGroupByName (const std::string &name) {
	static fgnm_t &feature_groups = feature_groups_map();

	fgnm_t::const_iterator fgnm_it = feature_groups.find(name);
	if (fgnm_it != feature_groups.end()) return (fgnm_it->second);

	size_t found;
	
//printf ("groupname cache miss\n");

// parse out algorithm name: everything up to the first '(' without trailing whitespace
	std::string algorithmname;
	FeatureAlgorithm *algorithm;
	size_t found_parens = name.find_first_of ('(');
	algorithmname.assign (name,0,found_parens);
	// clean trailing whitespace
	found=algorithmname.find_last_not_of(" \t\f\v\n\r");
	if (found != std::string::npos)
		algorithmname.erase(found+1);

	algorithm = (FeatureAlgorithm *) getFeatureAlgorithmByName (algorithmname);
	assert (algorithm != NULL && "algorithm not found when calling FeatureNames::getGroupByName");

// parse out the transforms - separated by '('
	size_t found_trans_s;
	size_t found_trans_e;
	std::string transform_name;
	std::vector<ImageTransform const *> transforms;
	const ImageTransform *transform;

	found_trans_s = name.find_first_not_of(" ()",found_parens+1);
	if (found_trans_s != std::string::npos) found_trans_e = name.find_first_of('(',found_trans_s+1);
	else (found_trans_e = std::string::npos);
	while (found_trans_e != std::string::npos) {
		transform_name.assign (name.substr (found_trans_s,found_trans_e-found_trans_s));
	// clean trailing whitespace
		found=transform_name.find_last_not_of(" \t\f\v\n\r");
		if (found != std::string::npos) {
			transform_name.erase(found+1);
			transform = getTransformByName(transform_name);
			transforms.push_back ( transform );
		}
		found_trans_s = name.find_first_not_of(" ()",found_trans_e+1);
		if (found_trans_s != std::string::npos) found_trans_e = name.find_first_of('(',found_trans_s+1);
		else (found_trans_e = std::string::npos);
	}
// The vector holds the transforms in application order, which opposite of left-to-right read order.
	std::reverse(transforms.begin(),transforms.end());

// Parse the channel
// If there is a channel specified, its in the inner parens.
// found_trans_s points at its first char, and found_trans_e points at npos
	std::string channel_name;
	Channel *channel = NULL;
	if (found_trans_s != std::string::npos && (found_trans_e = name.find_first_of(')',found_trans_s+1)) != std::string::npos ) {
		channel_name.assign (name.substr (found_trans_s,found_trans_e-found_trans_s));
		found=channel_name.find_last_not_of(" \t\f\v\n\r");
		if (found != std::string::npos)	{
			channel_name.erase(found+1);
			channel = (Channel *)getChannelByName(channel_name);
	// Empty parens - or not closed
		} else {
			channel = NULL;
		}
// Empty parens
	} else {
		channel = NULL;
	}

// Normalize the whitespace in the group name
	std::string name_norm = algorithm->name;
	size_t i;
	name_norm += " (";
	for (i = transforms.size(); i > 0  ; i--) {
		name_norm += transforms[i-1]->name + " (";
	}
	if (channel) name_norm += channel->name;
	name_norm += ")";
	for (i = 0; i < transforms.size(); i++) name_norm += ")";

// end of validation checks
// 	printf ("%s: [%s]",name.c._str(),featuregroup->algorithm.c_str());
// 	if (featuregroup->transforms.size()) printf (" [%s]",featuregroup->transforms[featuregroup->transforms.size()-1].c_str());
// 	for (int i= featuregroup->transforms.size()-2; i >= 0 ; i--) printf ("->[%s]",featuregroup->transforms[i].c_str());
// 	printf ("[%d]\n",index);

	return (feature_groups[name] = new FeatureGroup (name_norm, algorithm, channel, transforms));
}


const std::string &FeatureNames::oldFeatureNameLookup (const std::string &featurename_in) {
	static ofnm_t &old_features = old_features_map();
	const static std::string emptyString = "";

	ofnm_t::const_iterator ofnm_it = old_features.find(featurename_in);

	if (ofnm_it == old_features.end()) return (emptyString);
	else return (ofnm_it->second);
}



// Storage for class statics
// This is to avoid the "static initialization order fiasco" (see http://www.parashift.com/c++-faq/static-init-order-on-first-use-members.html)
FeatureNames::cnm_t  &FeatureNames::channels_map () {
	static cnm_t* channels_ = new cnm_t();
	return *channels_;
}
FeatureNames::tnm_t  &FeatureNames::transforms_map () {
	static tnm_t* transforms_ = new tnm_t();
	return *transforms_;
}
FeatureNames::fam_t  &FeatureNames::feature_algorithms_map () {
	static fam_t* feature_algorithms_ = new fam_t();
	return *feature_algorithms_;
}
FeatureNames::fgnm_t &FeatureNames::feature_groups_map () {
	static fgnm_t* feature_groups_ = new fgnm_t();
	return *feature_groups_;
}
FeatureNames::fnm_t  &FeatureNames::feature_names_map () {
	static fnm_t* feature_names_ = new fnm_t();
	return *feature_names_;
}
FeatureNames::ofnm_t &FeatureNames::old_features_map () {
	static ofnm_t* old_features_ = new ofnm_t();
	if (!old_features_->empty()) return (*old_features_);

	char *p = oldFeatureNamesFileStr;
	enum {st_eol, st_leading_space, st_ignore_line, st_key, st_val} state = st_eol, new_state = st_eol;
	std::string key, val;

	while (*p) {
		key = val = "";
		while (*p && state != st_key) {
			if (*p == '\n') new_state = st_eol;
			else if ( (state == st_eol || state == st_leading_space) && isspace (*p)) new_state = st_leading_space;
			else if ( (state == st_eol || state == st_leading_space) && !isalpha (*p)) new_state = st_ignore_line;
			else if (state == st_eol || state == st_leading_space) new_state = st_key;
		
			state = new_state;
			if (state != st_key) p++;
		}

		// key is everything up until the first tab.
		while (*p && state != st_val) {
			if (*p == '\t') new_state = st_val;
			else key += *p;

			state = new_state;
			if (state != st_val) p++;
		}
		
		// consume whitespace
		while (*p && isspace (*p)) p++;

		// value is everything until '\n'
		while (*p && state != st_eol) {
			if (*p == '\n') new_state = st_eol;
			else val += *p;

			state = new_state;
			if (state != st_eol) p++;
		}
		if (key.length()) {
			(*old_features_)[key] = val;
		}
	}

	return (*old_features_);
}
