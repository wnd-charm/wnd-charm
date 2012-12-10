#!/usr/bin/env python
# -------- preamble to get the test data --------------------
from pychrm.TrainingSet import *
import os
import sys
# This script is in googlecode/pychrm/trunk/tests/
# test data in googlecode/wndchrm/tests
# get the directory three levels up from this script, then wndchrm/tests
my_dir = os.path.dirname(os.path.realpath(__file__))
if (len (sys.argv) > 1):
	test_dir = sys.argv[1]
	if not os.path.isdir (test_dir):
		print "The supplied path to the tests directory '{0}' is not a directory".format (test_dir)
		test_dir = None
else:
	test_dir = os.path.join (os.path.dirname(os.path.dirname(os.path.dirname(my_dir))),'wndchrm','tests')
	if not os.path.isdir (test_dir):
		print "The path to the tests directory relative to this script '{0}' is not a directory".format (test_dir)
		test_dir = None
if not (test_dir):
	print "The tests directory can be checked out from svn and the test re-run using the following commands:"
	print "svn checkout http://wnd-charm.googlecode.com/svn/wndchrm/tests tests"
	print "{0} tests".format (sys.argv[0])
	sys.exit(0)
from pychrm import __version__ as pychrm_version
# -------- END preamble to get the test data --------------------

test_name = "WND5 Classification"
max_diff_pass = 0.000001
max_mean_pass = 0.000001
test_sig = os.path.join (test_dir,'t1_s01_c05_ij-l_precalculated.sig')
test_fit = os.path.join (test_dir,'test-l.fit')
test_fit_wght = os.path.join (test_dir,'test_fit-l.weights')
test_tif = os.path.join (test_dir,'t1_s01_c05_ij.tif')

test_result_2873 = [0.033,0.967]
test_result_431 = [0.030,0.970]
test_result_200 = [0.017,0.983]


test_sigs = Signatures.NewFromSigFile( test_sig, test_tif )

ts = FeatureSet_Discrete.NewFromFitFile( test_fit )
ts.featurenames_list = FeatureNameMap.TranslateToNewStyle( ts.featurenames_list )
ts.Normalize()
test_wghts = FisherFeatureWeights.NewFromFile (test_fit_wght)
test_wghts.names = FeatureNameMap.TranslateToNewStyle( test_wghts.names )
ts = ts.FeatureReduce( test_wghts.names )
test_sigs.Normalize( ts )

result_2873 = DiscreteImageClassificationResult.NewWND5( ts, test_wghts, test_sigs )

test_wghts_431 = test_wghts.Threshold( 431 )
ts_431 = ts.FeatureReduce( test_wghts_431.names )
test_sigs_431 = test_sigs.FeatureReduce ( test_wghts_431.names )
result_431 = DiscreteImageClassificationResult.NewWND5( ts_431, test_wghts_431, test_sigs_431 )

test_wghts_200 = test_wghts.Threshold( 200 )
ts_200 = ts.FeatureReduce( test_wghts_200.names )
test_sigs_200 = test_sigs.FeatureReduce ( test_wghts_200.names )
result_200 = DiscreteImageClassificationResult.NewWND5( ts_200, test_wghts_200, test_sigs_200 )

epsilon = 0.00001
max_diff = 0.
sum_diff = 0.
num_diffs = 0.

for idx in range( ts.num_classes ):
	test_val = test_result_2873[idx]
	calc_val = round (result_2873.marginal_probabilities[idx],3)
	diff = abs(calc_val - test_val)
	sum_diff += diff
	num_diffs += 1.0
	if diff > max_diff:
		max_diff = diff
	if ( diff > epsilon):
		print "2873 features: class '{0}' calculated probability ({1}) differs from expected ({2}) by {3}".format (
			ts.classnames_list [idx], calc_val, test_val, diff )

	test_val = test_result_431[idx]
	calc_val = round (result_431.marginal_probabilities[idx],3)
	diff = abs(calc_val - test_val)
	sum_diff += diff
	num_diffs += 1.0
	if diff > max_diff:
		max_diff = diff
	if ( diff > epsilon):
		print "431 features: class '{0}' calculated probability ({1}) differs from expected ({2}) by {3}".format (
			ts.classnames_list [idx], calc_val, test_val, diff )

	test_val = test_result_200[idx]
	calc_val = round (result_200.marginal_probabilities[idx],3)
	diff = abs(calc_val - test_val)
	sum_diff += diff
	num_diffs += 1.0
	if diff > max_diff:
		max_diff = diff
	if ( diff > epsilon):
		print "200 features: class '{0}' calculated probability ({1}) differs from expected ({2}) by {3}".format (
			ts.classnames_list [idx], calc_val, test_val, diff )

print "{0} comparissons, maximum diff = {1}, mean = {2}".format ( int (num_diffs), max_diff, sum_diff / num_diffs )
if (max_diff > max_diff_pass or (sum_diff / num_diffs) > max_mean_pass):
	print "pychrm {0} {1} test: {2}".format (pychrm_version, test_name, 'FAIL')
else:
	print "pychrm {0} {1} test: {2}".format (pychrm_version, test_name, 'PASS')
