#!/usr/bin/env python
from pychrm.FeatureSet import *
import multiprocessing as mp
import logging

#================================================================
def CleanShmemCache():
	"""read the lock files in /tmp and acquire their shmem blocks so they can be destroyed"""
	#from tempfile import gettempdir
	from glob import glob
	import os

	#globcmd = gettempdir() + sep + 'wndchrm*'
	globcmd = '/tmp/wndchrm*'
	lock_files = glob( globcmd )
	print len(lock_files)
	for lock_file in lock_files:
		tmp_path, tmp_file = os.path.split( lock_file )
		print tmp_file
		shimmat = pychrm.SharedImageMatrix()
		shimmat.fromCache( '/'+tmp_file )	


#================================================================
def worker( input_task_queue, output_result_queue ):
	for func, args in iter( input_task_queue.get, 'STOP' ):
		result = func( *args )
		output_result_queue.put( result )

#================================================================
def ConcurrentTransformFunc( tform_name, input_shmem_name ):
	"""returns the mmap path to the transformed pixel plane """

	mp.log_to_stderr()
	logger = mp.get_logger()
	logger.setLevel( logging.INFO )
	import os

	required_tform = None
	if isinstance( tform_name, tuple ):
		required_tform = tform_name[-1]
	else:
		required_tform = tform_name

	retval = None
	try:
		print "<<<<<<<<<<<{} (pid {}): attempting transform {} on input {}".format(
				mp.current_process(), os.getpid(), required_tform, input_shmem_name )

		original = pychrm.SharedImageMatrix()
		original.fromCache( input_shmem_name )
		if original.csREAD != original.Status(): 
			raise ValueError( 'Could not build an SharedImageMatrix from {0}, check the file.'.\
												 format( input_shmem_name ) )

		ret_px_plane = pychrm.SharedImageMatrix()
		ret_px_plane.transform( original, Transforms[ required_tform ] )

		original.DisableDestructorCacheCleanup( True )
		ret_px_plane.DisableDestructorCacheCleanup( True )
		ret_px_plane_shmem_name = ret_px_plane.GetShmemName()

		print ">>>>>>>>>>>>{} (pid {}): input {}, transformed \"{}\" = {}".format(
				mp.current_process(), os.getpid(), input_shmem_name, tform_name, ret_px_plane_shmem_name )

		retval = ( tform_name, ret_px_plane_shmem_name )

	except KeyboardInterrupt:
		print "ConcurrentTransformFunc: Caught Keyboard interrupt inside child process {}".format(
				os.getpid() )
	finally:
		return retval

#================================================================
def ConcurrentFeatCalcFunc( input_shmem_name, alg_name, feature_array, offset):
	"""Perform the feature calculation and put into the feature matrix"""

	mp.log_to_stderr()
	logger = mp.get_logger()
	logger.setLevel( logging.INFO )

	this_fns_retval = False

	try:
		print "<<<<<<<<<<<{} (pid {}: attempting FEATURE CALCULATION {} on input {}".format(
				mp.current_process(), os.getpid(), alg_name, input_shmem_name )
		
		px_plane = pychrm.SharedImageMatrix()
		px_plane.fromCache( input_shmem_name )
		if px_plane.csREAD != px_plane.Status(): 
			raise ValueError( 'Could not build an SharedImageMatrix from {0}, check the file.'.\
												 format( input_shmem_name ) )

		px_plane.DisableDestructorCacheCleanup(True)
		alg = Algorithms[ alg_name ]
		num_features = alg.n_features
		feature_array[ offset : offset + num_features ]  = alg.calculate( px_plane )

		print ">>>>>>>>>>>>{} (pid {}: input {}, calculated features {}".format(
				mp.current_process(), os.getpid(), input_shmem_name, alg_name  )

		this_fns_retval = True

	except KeyboardInterrupt:
		print "ConcurrentFeatCalcFunc: Caught Keyboard interrupt inside child process {}".format(
				os.getpid() )
	finally:
		return this_fns_retval
		

#================================================================
def GenerateWorkPlan( featuregroup_strings ):
	"""identify the required transforms and the order in which they need to occur"""
	# FIXME: Two levels of transforms hardcoded for now
	# ASSUMPTION: featuregroup string list coming in contains only unique strings

	first_round_tforms = set()
	second_round_tforms = set()
	parsed_algorithms = []
	feature_names = []

	for fg_string in featuregroup_strings:
		if fg_string.startswith( '#' ):
			continue
		string_rep = fg_string.rstrip( ")" )
		parsed = string_rep.split( ' (' )
		alg_name = parsed[0]

		for i in range( Algorithms[ alg_name ].n_features ):
			feature_names.append( fg_string + ' [' + str(i) + ']' )
		tform_list = parsed[1:]
		try:
			tform_list.remove( "" )
		except ValueError:
			pass

		# Make sure all transforms mentioned are known to pychrm
		if len(tform_list) != 0:
			for tform in tform_list:
				if tform not in Transforms:
					raise KeyError( "Don't know about a transform named {0}".format( tform ) )

		if len( tform_list ) > 1:
			# reverse() does an in-place list reversal that returns None 
			tform_list.reverse()
			second_round_tforms.add( tuple( tform_list ) )
			parsed_algorithms.append( ( alg_name, ' '.join( tform_list ) )  )
		elif len( tform_list ) == 1:
			first_round_tforms.add( tform_list[0] )
			parsed_algorithms.append( (alg_name, tform_list[0]) )
		else:
			parsed_algorithms.append( (alg_name, "") )

	return first_round_tforms, second_round_tforms, parsed_algorithms, feature_names

#================================================================
def CalcFeatures():

	import sys

	input_tif = sys.argv[1]
	num_cpus = mp.cpu_count()

	# Create queues
	task_queue = mp.Queue()
	done_queue = mp.Queue()



	original = pychrm.SharedImageMatrix()
	if 1 != original.OpenImage( input_tif , 0, None, 0, 0 ):
		raise ValueError( 'Could not build an SharedImageMatrix from {0}, check the file.'.\
							 format( path_to_image ) )
	else:
		print "successfully loaded image"

	pixel_planes = {}
	pp_shmem_names = {}
	original_shmem_name = original.GetShmemName()
	pp_shmem_names[ "" ] = original_shmem_name

	work_order = large_featureset_featuregroup_strings.split( '\n' )
	first_round_tforms, second_round_tforms, algs, feat_names = GenerateWorkPlan( work_order )

	# Enqueue first round tasks
	for tform_name in first_round_tforms:
		fn_args = ( tform_name, original_shmem_name )
		task_queue.put( ( ConcurrentTransformFunc, fn_args ) )

	# Start worker processes
	for i in range( num_cpus ):
		mp.Process( target = worker, args = (task_queue, done_queue) ).start()

	# Deal with the returned vals as they come off the done_queue
	for i in range( len( first_round_tforms ) ):
		tform_name, tformed_shmem_name = done_queue.get()
		new_ShImMat = pychrm.SharedImageMatrix()
		new_ShImMat.fromCache( original_shmem_name, tform_name )
		status = new_ShImMat.Status()
		if new_ShImMat.csREAD != status: 
			raise RuntimeError( 'Mother process got a csstatus of {} which should be csREAD'.\
												 format( status ) )
		assert( tformed_shmem_name == new_ShImMat.GetShmemName() )
		pixel_planes[ tform_name ] = new_ShImMat
		pp_shmem_names[ tform_name ] = tformed_shmem_name

	print "\n\n\n************************************ROUND 1 COMPLETE*************************"
	print "Round 1 Summary:"
	for tform_name in pp_shmem_names:
		print tform_name, "\t", pp_shmem_names[ tform_name ]
	print "\n\n\n"

	# Enqueue second round tasks
	for first_tform_name, second_tform_name in second_round_tforms:
		print "performing compound tform: ", first_tform_name, second_tform_name
		first_tform_shmem_addr = pp_shmem_names[ first_tform_name ]
		fn_args = ( (first_tform_name, second_tform_name), first_tform_shmem_addr )
		task_queue.put( ( ConcurrentTransformFunc, fn_args ) )

	# Deal with the returned vals as they come off the done_queue
	for i in range( len( second_round_tforms ) ):
		tform_pair, tformed_shmem_name = done_queue.get()
		new_ShImMat = pychrm.SharedImageMatrix()
		new_ShImMat.fromCache( pp_shmem_names[ tform_pair[0] ] , tform_pair[1] )
		status = new_ShImMat.Status()
		if new_ShImMat.csREAD != status: 
			raise RuntimeError( 'Mother process got a csstatus of {} which should be csREAD'.\
												 format( status ) )
		assert( new_ShImMat.GetShmemName() ==  tformed_shmem_name )
		compund_tform_name = ' '.join( tform_pair )
		pixel_planes[ compound_tform_name ] = new_ShImMat
		pp_shmem_names[ compound_tform_name ] = tformed_shmem_name 

	print "\n\n\n************************************ROUND 2 COMPLETE*************************"
	print "Round 2 Summary:"
	for tform_name in pp_shmem_names:
		print tform_name, "\t", pp_shmem_names[ tform_name ]
	print "\n\n\n"

	print "*****************TOTAL NUM TRANSFORMS {}****************".format( len( pixel_planes ) )
	# After all transforms are completed, all dependencies have been removed.
	# Fire off all feature calculation

	image_index = 0 # will be incremented when we start looping over images
	num_images = 1
	num_features = len( feat_names )
	import ctypes 

	shared_array_base = mp.Array( ctypes.c_double, num_images * num_features )

	column_offset = 0
	for algname, required_tform in algs:
		offset = image_index * num_features + column_offset
		required_pp_shmem_name = pp_shmem_names[ required_tform ]
		fn_args = ( required_pp_shmem_name, algname, shared_array_base, offset )
		task_queue.put( ( ConcurrentFeatCalcFunc, fn_args ) )
		column_offset += Algorithms[ algname ].n_features

	for i in range( len( algs ) ):
		done_queue.get()

	for shimmat in pixel_planes:
		pixel_planes[ shimmat ].DisableDestructorCacheCleanup(False)
	print dir( shared_array_base )

#================================================================
if __name__ == '__main__':
	CleanShmemCache()
	CalcFeatures()
