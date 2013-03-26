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
#include <assert.h>
#include <string>
#include <iostream>
#include "Tasks.h"
#include "FeatureNames.h"
#include "ImageTransforms.h"
#include "FeatureAlgorithms.h"
#include "cmatrix.h"


/* global variable */
extern int verbosity;

const char *ComputationTask::typeLabels (size_t type_idx) {
	const char *type_labels[] = {
		"Unknown Task",
		"Image Transform Task",
		"Feature Algorithm Task"
	};
	assert (type_idx > 0 && type_idx < (sizeof(type_labels) / sizeof (type_labels[0])) && "Attempt to get out of range TaskType label");
	return type_labels[type_idx];
}

void ComputationTask::print_info () const {
	std::cout << typeLabel () << " '" << name << "'" << std::endl;
}


bool ComputationTaskInstances::initialized () {
	static std::vector<const ComputationTask *> &instances = getInstances();
	return (! instances.empty());
}

// storage for static vector of ComputationTask instances in a singleton class
// Done in a static member function holding a static to avoid "static initialization order fiasco"
std::vector<const ComputationTask *> &ComputationTaskInstances::getInstances () {
	static std::vector<const ComputationTask *> *ComputationTasks = new std::vector<const ComputationTask *>;
	return (*ComputationTasks);
}

bool ComputationTaskInstances::add (const ComputationTask *task) {
	static std::vector<const ComputationTask *> &instances = getInstances();

	if (verbosity > 4) std::cout << "Registering " << task->typeLabel () << " '" << task->name << "'" << std::endl;
	instances.push_back (task);
	return (task->register_task());
}

ComputationTaskInstances::~ComputationTaskInstances() {
	static std::vector<const ComputationTask *> &instances = getInstances();

	for (size_t idx = 0; idx < instances.size(); idx++) {
		delete (instances[idx]);
	}
	instances.clear();
}

void ComputationTaskNode::print_info () const {
	for (size_t level = 0; level < depth; level++) std::cout << "- ";
	std::cout << name << " dependents: " << num_dependent_nodes;

	if (task)
		std::cout << ": " << task->typeLabel () << " '" << task->name << "'";
	std::cout << std::endl;

	for (size_t i = 0; i < dependent_tasks.size(); i++) {
		dependent_tasks[i]->print_info();
	}
}

size_t ComputationTaskNode::get_num_dependent_nodes () {
	num_dependent_nodes = 0;
	
	// We cast away the constness of the dependent tasks so we can update them.
	// Don't want to modify the declaration of the class to make these non-const.
	// This is trickier to do before nodes are added as dependent nodes
	for (size_t i = 0; i < dependent_tasks.size(); i++) {
		num_dependent_nodes += const_cast<ComputationTaskNode *>(dependent_tasks[i])->get_num_dependent_nodes();
		num_dependent_nodes++; // don't forget the dependent task itself
	}
	return (num_dependent_nodes);
}


void FeatureComputationPlan::add (const FeatureGroup *fg) {
	const ComputationTaskNode *source_node = root;
	nodemap_t::iterator nodemap_it;
	std::string node_key;

	assert(!isFinal && "Attempt to call FeatureComputationPlan::add() to a finalized plan");
	assert(fg && "Attempt to call FeatureComputationPlan::add() with NULL feature group");
	assert(FG_offset_map.find (fg->name) == FG_offset_map.end() && "Calling FeatureComputationPlan::add() with duplicate FeatureGroup");

	// Add the feature group to the list in the order they are added.
	feature_groups.push_back (fg);

	// Add nodes as necessary for transform dependencies
	// N.B.: add_get_node stores an internal reference to the node in the base class nodemap
	std::string trans_node_name;
	for (size_t i = 0; i < fg->transforms.size(); i++) {
		node_key += fg->transforms[i]->name;
		if (!trans_node_name.empty()) trans_node_name += "->";
		trans_node_name += fg->transforms[i]->name;
		source_node = add_get_node (node_key, source_node, fg->transforms[i], trans_node_name);
	}

	// Add the feature algorithm
	node_key += fg->algorithm->name;
	FG_node_map[fg->name] = add_get_node (node_key, source_node, fg->algorithm, fg->name);

	// Determine the column where to put this FG's results, and update the feature count.
	size_t start_idx = n_features;
	n_features += fg->algorithm->n_features;
	FG_offset_map[fg->name] = start_idx;
	
	// Add a lookup by index, mapping to the FG labels as well as the FG itself
	for (size_t idx = 0; idx < fg->labels.size(); idx++) {
		offset_FG_map[start_idx + idx] = fg;
		offset_FN_map[start_idx + idx] = &(fg->labels[idx]);
	}
}

void FeatureComputationPlan::add (const std::string &fg_name) {
	add ( FeatureNames::getGroupByName (fg_name) );
};

const std::string &FeatureComputationPlan::getFeatureName (size_t offset) const {
	static const std::string null_string = "";
	offset_FN_map_t::const_iterator it = offset_FN_map.find (offset);
	if (it != offset_FN_map.end()) {
		return *(it->second);
	} else {
		return null_string;
	}
}


bool compare_dependencies (const ComputationTaskNode *first, const ComputationTaskNode *second) {
	return first->num_dependent_nodes < second->num_dependent_nodes;
}
void ComputationPlanExecutor::make_dependencies_executable (const ComputationTaskNode *exec_node) {
	// add dependencies to executable_nodes
	executable_nodes.insert (executable_nodes.end(), exec_node->dependent_tasks.begin(), exec_node->dependent_tasks.end());
	make_heap (executable_nodes.begin(), executable_nodes.end(), compare_dependencies);
};

const ComputationTaskNode *ComputationPlanExecutor::get_next_executable_node () {
// get the first node off the heap and pop it off the vector

	pop_heap (executable_nodes.begin(), executable_nodes.end(), compare_dependencies);
	const ComputationTaskNode *exec_node = executable_nodes.back();
	executable_nodes.pop_back();
	return (exec_node);
}

void FeatureComputationPlanExecutor::execute_node (const ComputationTaskNode *exec_node) {
	// Put it in the executing nodes set
	ComputationPlanExecutor::execute_node (exec_node);

	const ComputationTask *task = exec_node->task;
	const ImageMatrix *IM_in = IM_map[exec_node->source_task->node_key];
	assert (IM_in != NULL && "Attempt to execute a FeatureComputationPlan node with a NULL source ImageMatrix");

	if (verbosity > 5) std::cout << "** executing node '" << exec_node->name << "' with " << exec_node->num_dependent_nodes << " total dependents. IM_in=" << IM_in;
	switch (task->type) {
		case ComputationTask::ImageTransformTask: {
			const ImageTransform *IT_task = dynamic_cast<const ImageTransform *>(exec_node->task);
			assert (IT_task && "Attempt to cast task as a (const ImageTransform *) failed.");
			// The ImageMatrix cache is keyed by node_key
			assert (IM_map.find(exec_node->node_key) == IM_map.end() && "Attempt to execute a transform which is already cached.");
			
			ImageMatrix *IM_out = new ImageMatrix;
			if (verbosity > 5) std::cout << " ImageTransform task '" << IT_task->name << "'" << std::endl;
			IT_task->execute (*IM_in, *IM_out);
			IM_map[exec_node->node_key] = IM_out;
		} break;
		
		case ComputationTask::FeatureAlgorithmTask: {
			const FeatureAlgorithm *FA_task = dynamic_cast<const FeatureAlgorithm *>(exec_node->task);
			assert (FA_task && "Attempt to cast task as a (const FeatureAlgorithm *) failed.");
			// The name of a node with a FeatureAlgorithmTask is the FeatureGroup name
			size_t offset = plan->getFGoffset (exec_node->name);
			assert (offset != (size_t)-1 && "execute_node() called with an unknown offset for the Feature Group");

			if (verbosity > 5) std::cout << " FeatureAlgorithm task '" << FA_task->name << "' (@ " << offset << ", " << FA_task->n_features << " features)" << std::endl;
			size_t mat_offset = (plan->n_features * current_feature_mat_row) + offset;
			// FIXME: we need to pass the result matrix chunck as a parameter instead of on the stack
			// construct a vector for the result
			std::vector<double> res_vec = FA_task->execute (*IM_in);
			for (int idx = 0; idx < FA_task->n_features; idx++) feature_mat[mat_offset+idx] = res_vec[idx];
		} break;
		
		default:
			std::cout << std::endl;
			assert (false && "Attempt to execute a node with an undefined task type");
		break;
	}
}

// FIXME: this can go into the base class (?) if its not specialized for task types
void FeatureComputationPlanExecutor::finish_node_execution (const ComputationTaskNode *exec_node) {
	if (verbosity > 6) std::cout << "finished '" << exec_node->name << "'" << std::endl;
	// remove it from executing_nodes
	ComputationPlanExecutor::finish_node_execution (exec_node);

	// N.B.:  The execute_node() method is responsible for storing the execution results
	make_dependencies_executable (exec_node);

	if (verbosity > 6) {
		if (exec_node->dependent_tasks.size()) {
		std::cout << "inserting into executable_nodes:" << std::endl;
		for (size_t i = 0; i < exec_node->dependent_tasks.size(); i++)
		std::cout << "  *" << exec_node->dependent_tasks[i]->name << std::endl;
		}
	}
}


void FeatureComputationPlanExecutor::run (const ImageMatrix *source_mat, std::vector<double> &feature_mat_in, size_t dest_row) {

	reset();

	feature_mat = &feature_mat_in[0];
	current_feature_mat_row = dest_row;
	// put the source_mat into the cache
	IM_map["root"] = source_mat;

	finish_node_execution(plan->root);

	const ComputationTaskNode *exec_node;
	
	while (! executable_nodes.empty() ) {
		exec_node = get_next_executable_node();
		execute_node (exec_node);
		finish_node_execution(exec_node);
	}
	// The caches get cleaned up in reset() above, or in the destructor
	if (verbosity > 5) std::cout << "Finished running execution plan '" << plan->name << "'" << std::endl;
}

void FeatureComputationPlanExecutor::reset () {
	// The ImageMatrixes in IM_map were all created within the execution, so they must all be deleted.
	// EXCEPT the root node, which was a parameter to run().
	// The root node must be called 'root'
	IM_map.erase ("root");
	IM_map_t::iterator IM_map_it;
	for(IM_map_it = IM_map.begin(); IM_map_it != IM_map.end(); IM_map_it++) {
		if (verbosity > 7) std::cout << "FeatureComputationPlanExecutor::reset(): deleting IM for node_key=" << IM_map_it->first << std::endl;
		delete (IM_map_it->second);
	}
	IM_map.clear();
	feature_mat = NULL;
	current_feature_mat_row = size_t(-1);
	// note that the plan stays.
}

const FeatureComputationPlan *StdFeatureComputationPlans::getFeatureSet () {
	static FeatureComputationPlan *the_plan = new FeatureComputationPlan ("Standard Feature Set");
	if ( the_plan->isFinalized() ) return the_plan;

	addStdFeatures (the_plan);
	the_plan->feature_vec_type = fv_short;
	the_plan->finalize();
	return (the_plan);
}

const FeatureComputationPlan *StdFeatureComputationPlans::getFeatureSetColor () {
	static FeatureComputationPlan *the_plan = new FeatureComputationPlan ("Color Feature Set");
	if ( the_plan->isFinalized() ) return the_plan;

	addStdFeatures (the_plan);
	addColorFeatures (the_plan);
	the_plan->feature_vec_type = fv_short_color;
	the_plan->finalize();
	return (the_plan);
}

const FeatureComputationPlan *StdFeatureComputationPlans::getFeatureSetLong () {
	static FeatureComputationPlan *the_plan = new FeatureComputationPlan ("Long Feature Set");
	if ( the_plan->isFinalized() ) return the_plan;

	addLongFeatures (the_plan, false);
	the_plan->feature_vec_type = fv_long;
	the_plan->finalize();
	return (the_plan);
}

const FeatureComputationPlan *StdFeatureComputationPlans::getFeatureSetLongColor () {
	static FeatureComputationPlan *the_plan = new FeatureComputationPlan ("Long Color Feature Set");
	if ( the_plan->isFinalized() ) return the_plan;

	addLongFeatures (the_plan, true);
	the_plan->feature_vec_type = fv_long_color;
	the_plan->finalize();
	return (the_plan);
}

void StdFeatureComputationPlans::addLongFeatures (FeatureComputationPlan *the_plan, bool color) {

	addGroupAFeatures (the_plan, "()");
	addGroupBFeatures (the_plan, "()");
	addGroupCFeatures (the_plan, "()");

	if (color)
		addColorFeatures (the_plan);

	addGroupBFeatures (the_plan, "(Fourier ())");
	addGroupCFeatures (the_plan, "(Fourier ())");

	addGroupBFeatures (the_plan, "(Wavelet ())");
	addGroupCFeatures (the_plan, "(Wavelet ())");

	addGroupBFeatures (the_plan, "(Chebyshev ())");
	addGroupCFeatures (the_plan, "(Chebyshev ())");

	// Fourier, then Chebyshev
	addGroupCFeatures (the_plan, "(Chebyshev (Fourier ()))");

	// Fourier, then Wavelet
	addGroupCFeatures (the_plan, "(Wavelet (Fourier ()))");

	// Wavelet, then Fourier
	addGroupBFeatures (the_plan, "(Fourier (Wavelet ()))");
	addGroupCFeatures (the_plan, "(Fourier (Wavelet ()))");

	// Chebyshev, then Fourier
	addGroupCFeatures (the_plan, "(Fourier (Chebyshev ()))");

	// Wavelet, then Chebyshev
	addGroupCFeatures (the_plan, "(Chebyshev (Wavelet ()))");

	// Edge
	addGroupBFeatures (the_plan, "(Edge ())");
	addGroupCFeatures (the_plan, "(Edge ())");

	// Edge, then Fourier
	addGroupBFeatures (the_plan, "(Fourier (Edge ()))");
	addGroupCFeatures (the_plan, "(Fourier (Edge ()))");

	// Edge, then wavelet
	addGroupBFeatures (the_plan, "(Wavelet (Edge ()))");
	addGroupCFeatures (the_plan, "(Wavelet (Edge ()))");

}

void StdFeatureComputationPlans::addGroupAFeatures (FeatureComputationPlan *the_plan, std::string transform) {
	const char* the_fs[] = {
	// Group A:
		"Edge Features ",
		"Otsu Object Features ",
		"Inverse-Otsu Object Features ",
		"Gabor Textures "
	};
	for (size_t i = 0; i < ( sizeof(the_fs) / sizeof(the_fs[0]) ); i++) {
		the_plan->add(std::string(the_fs[i]) + transform);
	}
}

void StdFeatureComputationPlans::addGroupBFeatures (FeatureComputationPlan *the_plan, std::string transform) {
	const char* the_fs[] = {
	// Group B:
		"Chebyshev-Fourier Coefficients ",
		"Chebyshev Coefficients ",
		"Zernike Coefficients "
	};
	for (size_t i = 0; i < ( sizeof(the_fs) / sizeof(the_fs[0]) ); i++) {
		the_plan->add(std::string(the_fs[i]) + transform);
	}
}

void StdFeatureComputationPlans::addGroupCFeatures (FeatureComputationPlan *the_plan, std::string transform) {
	const char* the_fs[] = {
	// Group C:
		"Comb Moments ",
		"Haralick Textures ",
		"Multiscale Histograms ",
		"Tamura Textures ",
		"Radon Coefficients ",
		"Fractal Features ",   
		"Pixel Intensity Statistics ",   
		"Gini Coefficient "
	};
	for (size_t i = 0; i < ( sizeof(the_fs) / sizeof(the_fs[0]) ); i++) {
		the_plan->add(std::string(the_fs[i]) + transform);
	}
}


void StdFeatureComputationPlans::addColorFeatures (FeatureComputationPlan *the_plan) {
	the_plan->add("Color Histogram ()");

	addGroupBFeatures (the_plan, "(Color Transform ())");
	addGroupCFeatures (the_plan, "(Color Transform ())");

	addGroupBFeatures (the_plan, "(Hue ())");
	addGroupCFeatures (the_plan, "(Hue ())");

	addGroupBFeatures (the_plan, "(Fourier (Hue ()))");
	addGroupCFeatures (the_plan, "(Fourier (Hue ()))");

	addGroupBFeatures (the_plan, "(Chebyshev (Hue ()))");
	addGroupCFeatures (the_plan, "(Chebyshev (Hue ()))");
}


void StdFeatureComputationPlans::addStdFeatures (FeatureComputationPlan *the_plan) {

	const char* the_fs[] = {
		"Chebyshev-Fourier Coefficients ()",
		"Chebyshev-Fourier Coefficients (Fourier ())",
		"Chebyshev Coefficients ()",
		"Chebyshev Coefficients (Fourier ())",
		"Comb Moments ()",
		"Comb Moments (Chebyshev ())",
		"Comb Moments (Chebyshev (Fourier ()))",
		"Comb Moments (Fourier ())",
		"Comb Moments (Wavelet ())",
		"Comb Moments (Wavelet (Fourier ()))",
		"Edge Features ()",
		"Otsu Object Features ()",
		"Inverse-Otsu Object Features ()",
		"Gabor Textures ()",
		"Haralick Textures ()",
		"Haralick Textures (Chebyshev ())",
		"Haralick Textures (Chebyshev (Fourier ()))",
		"Haralick Textures (Fourier ())",
		"Haralick Textures (Wavelet ())",
		"Haralick Textures (Wavelet (Fourier ()))",
		"Multiscale Histograms ()",
		"Multiscale Histograms (Chebyshev ())",
		"Multiscale Histograms (Chebyshev (Fourier ()))",
		"Multiscale Histograms (Fourier ())",
		"Multiscale Histograms (Wavelet ())",
		"Multiscale Histograms (Wavelet (Fourier ()))",
		"Radon Coefficients ()",
		"Radon Coefficients (Chebyshev ())",
		"Radon Coefficients (Chebyshev (Fourier ()))",
		"Radon Coefficients (Fourier ())",
		"Tamura Textures ()",
		"Tamura Textures (Chebyshev ())",
		"Tamura Textures (Chebyshev (Fourier ()))",
		"Tamura Textures (Fourier ())",
		"Tamura Textures (Wavelet ())",
		"Tamura Textures (Wavelet (Fourier ()))",
		"Zernike Coefficients ()",
		"Zernike Coefficients (Fourier ())"
	};

	for (size_t i = 0; i < ( sizeof(the_fs) / sizeof(the_fs[0]) ); i++) {
		the_plan->add(the_fs[i]);
	}
}


