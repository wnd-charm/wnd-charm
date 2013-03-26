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
#ifndef __TASKS_H_
#define __TASKS_H_

#include <assert.h>
#include <vector>
#include <string>
#include <set>
// defines OUR_UNORDERED_MAP based on what's available
#include "unordered_map_dfn.h"

// ComputationTask
// 	hasa registration functionality
// 	hasa name
// FeatureAlgorithm
// 	isa ComputationTask
// ImageTransform
// 	isa ComputationTask
// 
// ComputationTaskNode
// 	hasa ComputationTask source_task
// 	hasmany ComputationTask dependent_tasks
// ComputationPlan
// 	hasmany ComputationTaskNode
// 	hasa vector of FeatureGroup
// 	add (FeatureGroup)
// 	start (src_mat, dest_mat, dest_row) // calculates column offsets
// 	run()
// 
// const static char[] ShortFeatureGroups
// const static char[] LongFeatureGroups
// const static char[] ColorFeatureGroups
// "static" ComputationPlan LongFeatureVector
// "static" ComputationPlan LongColorFeatureVector
// "static" ComputationPlan ShortFeatureVector
// "static" ComputationPlan ShortColorFeatureVector

class FeatureComputationPlan;
class ComputationTaskNode;

// This is the parent class for classes that implement algorithms that do stuff.
// The parent class is virtual. Inherited classes must implement execute() and register_task()
class ComputationTask {
	public:
		enum TaskType {
			UnknownTaskType = 0,
			ImageTransformTask = 1,
			FeatureAlgorithmTask = 2,
		};
		std::string name;
		TaskType type;

		virtual bool register_task() const = 0;

		virtual void print_info() const;
		static const char *typeLabels (size_t type_idx);
		const char *typeLabel () const {return (typeLabels(type));}
	// Protected so inherited classes can still call the parent constructor
	protected:
		ComputationTask (const std::string &s, TaskType t) { name = s; type = t;}
	private:
        ComputationTask(ComputationTask const&);              // Don't Implement
        void operator=(ComputationTask const&); // Don't implement
};

// This singleton class owns all of the instances of ComputationTask singletons.
//  (i.e. the actual FeatureAlgorithms. ImageTransforms, etc. that "do" stuff)
// so that we can get references (pointers) to them later for looking them up by name using a map.
// Done in a static member function holding a static to avoid "static initialization order fiasco"
class ComputationTaskInstances {
	public:
		static bool initialized ();
		static bool add (const ComputationTask *task);
		static std::vector<const ComputationTask *> &getInstances ();
		~ComputationTaskInstances();
	private:
		ComputationTaskInstances(); // private constructor makes this a static class
        ComputationTaskInstances(ComputationTaskInstances const&); // Don't Implement
        void operator=(ComputationTaskInstances const&);           // Don't implement
};

// A node in the ComputationPlan
// The root node may or may not have a task defined.
// All other nodes have a defined source_task and a defined task.
class ComputationTaskNode {
	public:
		const ComputationTaskNode *source_task; // root = NULL
		const ComputationTask *task;
		std::string name;
		std::string node_key;
		size_t num_dependent_nodes; // full traversal, valid only after the ComputationPlan is finalized
		size_t depth; // root = 0
		// The nodes are owned by a ComputationPlan
		std::vector <const ComputationTaskNode *> dependent_tasks;
		void print_info() const;
		size_t get_num_dependent_nodes ();
		ComputationTaskNode(ComputationTaskNode *source_node_in, const ComputationTask *task_in, const std::string &name_in = "", const std::string &key_in = "") {
			assert (! (task_in && !source_node_in) && "Attempt to create a ComputationTaskNode with an undefined source but defined task");
			assert (! (!task_in && source_node_in) && "Attempt to create a ComputationTaskNode with an undefined task but defined source");

			source_task = source_node_in; // NULL for root nodes
			task = task_in;               // Maybe NULL for root nodes
			num_dependent_nodes = 0;


			if (task_in && name_in.empty())
				name = task_in->name;
			else if (name_in.empty())
				name = "root"; // the root node must be named 'root'
			else
				name = name_in;
			// attach the node to the source
			if (source_node_in) { // not the root node
				source_node_in->dependent_tasks.push_back (this);
				depth = source_node_in->depth + 1;
				if (key_in.empty()) node_key = source_node_in->node_key + name;
				else node_key = key_in;
			} else { // the root node
				depth = 0;
				if (key_in.empty()) node_key = name; // the root node's key must be 'root'
				else node_key = key_in;
				// Hopefully, the caller also does something with a root node 
			}
		}
};


// This class maintains a tree of ComputationTaskNode objects.
// The idea is to keep the plan's nodes undisturbed because the plan can be applied many times over.
// The ComputationPlanExecutor class is used to run a specific plan.  This class adds members needed at run-time.
// ComputationPlan
//     has name, root, nodemap, isFinal
//     has virtual add (string)
//     has add_get_node
//     finalize()
// FeatureComputationPlan : ComputationPlan
//     has feature_groups vector
//     has FG_offset_map
//     has n_features
class ComputationPlan {
	public:
		std::string name;
		ComputationTaskNode *root;
		virtual void add (const std::string &node_name) = 0;

		// just a helper to check for a node_key and either return an existing node, or add a new one.
		// N.B.: the source_node parameter is passed in as const since we want nodes to be const outside of the API.
		// If we end up adding this node, we have to cast away this constness to add the new node as a dependency to the source_node.
		const ComputationTaskNode *add_get_node (const std::string &node_key, const ComputationTaskNode *source_node,
			const ComputationTask *task, const std::string &name_in = "") {
				nodemap_t::iterator nodemap_it;

				nodemap_it = nodemap.find(node_key);
				if (nodemap_it != nodemap.end()) {
					return (nodemap_it->second);
				} else {
					// Note the constness of the source node is cast away here
					return (nodemap[node_key] = new ComputationTaskNode(const_cast<ComputationTaskNode *>(source_node), task, name_in, node_key));
				}
		}

		// finalize() does whatever's necessary between the nodes being complete and the first execution.
		void finalize () {
			root->get_num_dependent_nodes();
			isFinal = true;
		}
		bool isFinalized() { return isFinal; }

		ComputationPlan (const std::string &name_in) {
			name = name_in;
			// the root isn't in the maps - this is its one reference.
			root = new ComputationTaskNode(NULL,NULL);
			isFinal = false;
		}
		~ComputationPlan() {
			nodemap_t::iterator nodemap_it;
			for(nodemap_it = nodemap.begin(); nodemap_it != nodemap.end(); nodemap_it++) {
				delete (nodemap_it->second);
			}
			nodemap.clear();
			delete root;
		}
	protected:
		bool isFinal;
		// node keys are composed of catenating dependent task names
		typedef OUR_UNORDERED_MAP<std::string, const ComputationTaskNode *> nodemap_t;
		nodemap_t nodemap;
	private:
		ComputationPlan();                       // Don't implement
        ComputationPlan(ComputationPlan const&); // Don't Implement
        void operator=(ComputationPlan const&);  // Don't implement
};

// ComputationPlanExecutor
//     has plan
//     has executable_nodes
//     has executing_nodes
//     virtual run()
//     virtual finish_node_execution()
//     virtual execute_node()
// FeatureComputationPlanExecutor : ComputationPlanExecutor
//     has IM_cache (IM_map)
// FeatureComputationPlanConcurrentExecutor : FeatureComputationPlanExecutor
//     has concurrency specialization (execute_node(), finish_node_execution())
//
// ComputationPlanExecutor is a virtual class for running specific types of plans
// The run() method is pure-virtual, which must be defined in derived classes.
// The run() method in derived classes usually needs parameters, so the non-parameter run() method
//   must also be overridden with a noop in this case.
class ComputationPlanExecutor {
	public:
		const ComputationPlan *plan;
		
		const ComputationTaskNode *get_next_executable_node ();
		// Call once per iteration, supplying iteration parameters.
		// Pure virtual - must have override in inherited class
		virtual void run () = 0;

		virtual void finish_node_execution (const ComputationTaskNode *exec_node) {
			executing_nodes.erase (exec_node->node_key);
		}

		// This should be called at the end of an over-ridden finish_node_execution()
		void make_dependencies_executable (const ComputationTaskNode *exec_node);

		// Sub-classes must assign their specific plan pointer in their constructor
		// Relying on the parent constructor to set the plan member doesn't result in the proper setting
		// of the sub-class plan member if it has a different pointer type than a plain ComputationPlan *
		// One of those C++ things...  something to do with shadowing members.  blech.
		// N.B.: It is very wrong to end up with an executor with no plan.  Not sure how to best enforce this in C++
		ComputationPlanExecutor(const ComputationPlan *plan_in) {
			plan = plan_in;
		}
		~ComputationPlanExecutor () {
			reset();
		}
	protected:
		// executable_nodes is a vector that gets heap-ified to act as a priority queue.
		// The only reason its not an actual std::priority_queue is that we want to add nodes to it in a block, then heapify.
		typedef std::vector<const ComputationTaskNode *> executable_nodes_t;
		executable_nodes_t executable_nodes;
		// The executing_nodes is an std::set of nodes (no order, just a list we can add/remove from)
		typedef OUR_UNORDERED_MAP<std::string, const ComputationTaskNode *> executing_nodes_t;
		executing_nodes_t executing_nodes;

		virtual void execute_node (const ComputationTaskNode *exec_node) {
			assert (executing_nodes.find(exec_node->node_key) == executing_nodes.end() && "Attempt to execute a node which is already executing.");
			executing_nodes.insert ( std::pair<std::string, const ComputationTaskNode *>(exec_node->node_key, exec_node) );
		}
		// This resets the object for the next call to run()
		// run() calls reset(), so this is protected.  The destructor also calls reset()
		// The base class doesn't create things during a run, so there's nothing to destroy here between runs.
		virtual void reset () {};
	private:
		ComputationPlanExecutor();                                // Don't implement
        ComputationPlanExecutor(ComputationPlanExecutor const&);  // Don't Implement
        void operator=(ComputationPlanExecutor const&);           // Don't implement
};


// 
// forward declarations
class ImageMatrix;
class FeatureGroup;
// This class has additional members and methods specific for computing features
#define CURRENT_FEATURE_VERSION 2
class FeatureComputationPlan : public ComputationPlan {
	public:
		size_t n_features;
		int feature_vec_type;              // stores the integer value of the feature_vec_types enum.

		virtual void add (const std::string &FGname);
		void add (const FeatureGroup *fg);

		size_t getFGoffset (const std::string &FGname) const {
			FG_offset_map_t::const_iterator it = FG_offset_map.find (FGname);
			if (it != FG_offset_map.end())
				return it->second;
			else
				return -1;
		}

		const std::string &getFeatureName (size_t offset) const;
		FeatureComputationPlan (const std::string &name_in) : ComputationPlan (name_in) {
			n_features = 0;
			feature_vec_type = 0;
		}
		// parent destructor takes care of CalculationTask objects
		// This plan doesn't own any of the objects it has references to
		~FeatureComputationPlan() {}
	private:
		std::vector<const FeatureGroup *> feature_groups;

		// FG_offset_map keys are feature group names. The value is the column where the FG vector starts.
		typedef OUR_UNORDERED_MAP<std::string, size_t> FG_offset_map_t;
		FG_offset_map_t FG_offset_map;
		// offset_FG_map keys are feature numbers/offsets, values are FG pointers
		typedef OUR_UNORDERED_MAP<size_t, const FeatureGroup *> offset_FG_map_t;
		offset_FG_map_t offset_FG_map;
		// offset_FN_map keys are feature numbers/offsets, values are individual feature names, which are references to FG.labels
		// can't store references in a map in c++, so they have to be pointers.
		// The strings pointed to are owned by the FG object they point to.
		typedef OUR_UNORDERED_MAP<size_t, const std::string *> offset_FN_map_t;
		offset_FN_map_t offset_FN_map;
		// FG_node_map keys are feature group names. The value is a pointer to the corresponding node object
		// FIXME: We're not really using this for anything - is it necessary?
		typedef OUR_UNORDERED_MAP<std::string, const ComputationTaskNode *> FG_node_map_t;
		FG_node_map_t FG_node_map;

};

class FeatureComputationPlanExecutor : public ComputationPlanExecutor {
	public:
		const FeatureComputationPlan *plan;
		double *feature_mat;
		size_t current_feature_mat_row;

		virtual void finish_node_execution (const ComputationTaskNode *exec_node);
		virtual void run (const ImageMatrix *source_mat, std::vector<double> &feature_mat_in, size_t dest_row);
		// in the parent, the run method signature has no parameters and is pure virtual
		// this class has to have run parameters, so we override the paren't virtual run() with a noop
		virtual void run () {}
		~FeatureComputationPlanExecutor () {
			reset();
		}
		// sub-classes must set their own plan.  Relying on the parent class to do this doesn't work.
		FeatureComputationPlanExecutor (const FeatureComputationPlan *plan_in) : ComputationPlanExecutor (plan_in) {
			plan = plan_in;
			feature_mat = NULL;
			current_feature_mat_row = size_t(-1);
		}
	protected:
		// ImageMatrix cache
		// IM_map keys are node_keys for transform nodes (source->node_key)
		typedef OUR_UNORDERED_MAP<std::string, const ImageMatrix *> IM_map_t;
		IM_map_t IM_map;

		virtual void execute_node (const ComputationTaskNode *exec_node);
		// This resets the object for the next call to run() (run() calls reset)
		virtual void reset ();

};
// two options for concurrent executors:  Derive from ConcurrentExecutor, or derive from plan-specific executor.
// class FeatureComputationPlanConcurrentExecutor : public ConcurrentPlanExecutor
//   -- and --  class ConcurrentPlanExecutor : public ComputationPlanExecutor
// -- or --
// class FeatureCalculationConcurrentExecutor : public FeatureComputationPlanExecutor

class StdFeatureComputationPlans {
	private:
		StdFeatureComputationPlans(); // private constructor: static class
	public:
		enum feature_vec_types {
			fv_unknown = 0,  // this must evaluate to false
			fv_short = 1,
			fv_long = 2,
			fv_short_color = 3,
			fv_long_color = 4
		};
		static const FeatureComputationPlan *getFeatureSet();
		static const FeatureComputationPlan *getFeatureSetColor();
		static const FeatureComputationPlan *getFeatureSetLong();
		static const FeatureComputationPlan *getFeatureSetLongColor();
		static void addLongFeatures (FeatureComputationPlan *the_plan, bool color);
		static void addGroupAFeatures (FeatureComputationPlan *the_plan, std::string transform);
		static void addGroupBFeatures (FeatureComputationPlan *the_plan, std::string transform);
		static void addGroupCFeatures (FeatureComputationPlan *the_plan, std::string transform);
		static void addColorFeatures (FeatureComputationPlan *the_plan);
		static void addStdFeatures (FeatureComputationPlan *the_plan);
};

#endif //__TASKS_H_
