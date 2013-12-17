#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2013 University of Dundee & Open Microscopy Environment.
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#
#

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

import re
import numpy as np

from pychrm.FeatureSet import DiscreteBatchClassificationResult, \
    FeatureSet_Discrete, FisherFeatureWeights, Signatures


class TestPychrm(unittest.TestCase):
    """
    Test the pychrm module to check it still works with OmeroPychrm
    """

    def setUp(self):
        self.image1 = 'test-0032-0008-0008.tif'
        self.image2 = 'test-0032-0016-0016.tif'
        self.feature_vector_version = '12345.45689'
        self.SmallFSexpectedLen = 1059
        self.LargeFSexpectedLen = 2919

    def createSignature(self, a, b):
        s = Signatures()
        s.names = ['ft [0]', 'ft [1]']
        s.values = [a, b]
        s.source_file = '%s %s' % (a, b)
        s.version = self.feature_vector_version
        return s

    def createSignatures(self):
        s1 = self.createSignature(1, 10)
        s2 = self.createSignature(2, 10)
        s3 = self.createSignature(3, 20)
        s4 = self.createSignature(4, 20)
        return s1, s2, s3, s4


    def test_featuresInvalidImagePath(self):
        nonExistentImage = 'non-existent-image.tif'
        self.assertRaises(
            ValueError, Signatures.SmallFeatureSet, nonExistentImage)
        self.assertRaises(
            ValueError, Signatures.LargeFeatureSet, nonExistentImage)

    def test_calculateSmallFeatureSet(self):

        ft1 = Signatures.SmallFeatureSet(self.image1)
        self.assertEqual(len(ft1.names), self.SmallFSexpectedLen)
        self.assertTrue(all([re.match('^.+ \[\d+\]$', n) for n in ft1.names]))
        self.assertEqual(len(set(ft1.names)), self.SmallFSexpectedLen)
        self.assertEqual(len(ft1.values), self.SmallFSexpectedLen)
        self.assertFalse(any(np.isinf(ft1.values)))
        self.assertFalse(any(np.isnan(ft1.values)))

        ft2 = Signatures.SmallFeatureSet(self.image2)
        self.assertEqual(ft1.names, ft2.names)
        self.assertEqual(len(ft2.values), self.SmallFSexpectedLen)
        self.assertFalse(any(np.isinf(ft2.values)))
        self.assertFalse(any(np.isnan(ft2.values)))

    def test_calculateLargeFeatureSet(self):

        ft1 = Signatures.LargeFeatureSet(self.image1)
        self.assertEqual(len(ft1.names), self.LargeFSexpectedLen)
        self.assertTrue(all([re.match('^.+ \[\d+\]$', n) for n in ft1.names]))
        self.assertEqual(len(set(ft1.names)), self.LargeFSexpectedLen)
        self.assertEqual(len(ft1.values), self.LargeFSexpectedLen)
        self.assertFalse(any(np.isinf(ft1.values)))
        self.assertFalse(any(np.isnan(ft1.values)))

        ft2 = Signatures.LargeFeatureSet(self.image2)
        self.assertEqual(ft1.names, ft2.names)
        self.assertEqual(len(ft2.values), self.LargeFSexpectedLen)
        self.assertFalse(any(np.isinf(ft2.values)))
        self.assertFalse(any(np.isnan(ft2.values)))

    def test_createFeatureSet(self):
        sig1, sig2, sig3, sig4 = self.createSignatures()
        fts = FeatureSet_Discrete()

        # Add classes out of order
        fts.AddSignature(sig3, 1)

        self.assertEqual(fts.num_classes, 2)
        self.assertEqual(fts.num_features, 2)
        self.assertEqual(fts.num_images, 1)
        self.assertEqual(len(fts.data_list), 2)
        self.assertIsNone(fts.data_list[0])
        #self.assertSequenceEqual
        np.testing.assert_almost_equal(fts.data_list[1], sig3.values)

        self.assertEqual(fts.classsizes_list, [0, 1])
        self.assertEqual(fts.classnames_list, ['UNKNOWN1', 'UNKNOWN2'])

        fts.AddSignature(sig1, 0)

        self.assertEqual(fts.num_classes, 2)
        self.assertEqual(fts.num_features, 2)
        self.assertEqual(fts.num_images, 2)
        self.assertEqual(len(fts.data_list), 2)
        np.testing.assert_almost_equal(fts.data_list[0], sig1.values)
        np.testing.assert_almost_equal(fts.data_list[1], sig3.values)

        self.assertEqual(fts.classsizes_list, [1, 1])
        self.assertEqual(fts.classnames_list, ['UNKNOWN1', 'UNKNOWN2'])

        # fts.ContiguousDataMatrix() fails unless there are at least two images
        # per class, is this really necessary?
        #tmp = fts.ContiguousDataMatrix()
        fts.AddSignature(sig2, 0)
        fts.AddSignature(sig4, 1)
        self.assertEqual(fts.classsizes_list, [2, 2])
        self.assertEqual(fts.num_images, 4)

        tmp = fts.ContiguousDataMatrix()
        self.assertEqual(fts.data_matrix.shape, (4, 2))
        np.testing.assert_almost_equal(fts.data_matrix[0], sig1.values)
        np.testing.assert_almost_equal(fts.data_matrix[1], sig2.values)
        np.testing.assert_almost_equal(fts.data_matrix[2], sig3.values)
        np.testing.assert_almost_equal(fts.data_matrix[3], sig4.values)

    def test_incompatibleFeatureVersion(self):
        s = self.createSignature(1, 10)
        fts = FeatureSet_Discrete()
        fts.feature_vector_version = '0.0'
        self.assertRaises(ValueError, fts.AddSignature, s, 1)

    def test_fisherFeatureWeights(self):
        sig1, sig2, sig3, sig4 = self.createSignatures()

        fts = FeatureSet_Discrete()
        fts.AddSignature(sig1, 0)
        fts.AddSignature(sig2, 0)
        fts.AddSignature(sig3, 1)
        fts.AddSignature(sig4, 1)
        tmp = fts.ContiguousDataMatrix()

        # TODO: weight[1]==0, presumably because the intra-class variance=0,
        # even though feature[1] is a perfect discriminator?
        fts.Normalize()

        wts = FisherFeatureWeights.NewFromFeatureSet(fts)

        np.testing.assert_almost_equal(wts.values, [4.0, 0.0])
        self.assertEqual(wts.names, ['ft [0]', 'ft [1]'])

    def test_thresholdWeights(self):
        w = FisherFeatureWeights()
        w.names = ['a', 'b', 'c']
        w.values = [1.0, 2.0, 4.0]

        w1 = w.Threshold(2)
        self.assertEqual(w1.names, ['c', 'b'])
        self.assertAlmostEqual(w1.values, [4.0, 2.0])

    def test_reduceFeatures(self):
        ftnames = ['c', 'b']

        s = Signatures()
        s.names = ['a', 'b', 'c']
        s.values = [1.0, 2.0, 3.0]
        s.source_file = ''

        s1 = s.FeatureReduce(ftnames)
        self.assertEqual(s1.names, ftnames)
        self.assertAlmostEqual(s1.values, [3.0, 2.0])


    def test_predictDiscrete(self):
        s1 = self.createSignature(1, -1)
        s2 = self.createSignature(2, -2)
        s3 = self.createSignature(3, -3)
        s4 = self.createSignature(4, -4)

        trainFts = FeatureSet_Discrete()
        trainFts.AddSignature(s1, 0)
        trainFts.AddSignature(s2, 0)
        trainFts.AddSignature(s3, 1)
        trainFts.AddSignature(s4, 1)
        tmp = trainFts.ContiguousDataMatrix()

        # Values are chosen so that different weights give different predictions
        s5 = self.createSignature(0, -5)
        s6 = self.createSignature(5, 0)

        testFts = FeatureSet_Discrete()
        testFts.AddSignature(s5, 0)
        testFts.AddSignature(s6, 0)
        tmp = testFts.ContiguousDataMatrix()

        weights = FisherFeatureWeights()
        weights.names = ['ft [0]', 'ft [1]']
        weights.values = [2.0, 1.0]

        pred = DiscreteBatchClassificationResult.New(trainFts, testFts, weights)
        self.assertEqual(len(pred.individual_results), 2)
        r1, r2 = pred.individual_results
        np.testing.assert_almost_equal(
            r1.marginal_probabilities, [0.975, 0.025], decimal=3)
        np.testing.assert_almost_equal(
            r2.marginal_probabilities, [0.025, 0.975], decimal=3)

        weights = FisherFeatureWeights()
        weights.names = ['ft [0]', 'ft [1]']
        weights.values = [1.0, 2.0]

        pred = DiscreteBatchClassificationResult.New(trainFts, testFts, weights)
        self.assertEqual(len(pred.individual_results), 2)
        r1, r2 = pred.individual_results
        np.testing.assert_almost_equal(
            r1.marginal_probabilities, [0.025, 0.975], decimal=3)
        np.testing.assert_almost_equal(
            r2.marginal_probabilities, [0.975, 0.025], decimal=3)



if __name__ == '__main__':
    unittest.main()
