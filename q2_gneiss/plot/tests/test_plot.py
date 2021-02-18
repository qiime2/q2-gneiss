# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import unittest

import numpy as np
import pandas as pd
from qiime2 import CategoricalMetadataColumn
from scipy.cluster.hierarchy import ward
from skbio import TreeNode, DistanceMatrix

from q2_gneiss.plot._plot import dendrogram_heatmap


class TestHeatmap(unittest.TestCase):

    def setUp(self):
        self.results = "results"
        if not os.path.exists(self.results):
            os.mkdir(self.results)

    def tearDown(self):
        shutil.rmtree(self.results)

    def test_visualization(self):
        np.random.seed(0)
        num_otus = 500  # otus
        index = pd.Index(np.arange(5).astype(np.str), name='id')
        table = pd.DataFrame(np.random.random((len(index), num_otus)),
                             index=index,
                             columns=np.arange(num_otus).astype(np.str))

        x = np.random.rand(num_otus)
        dm = DistanceMatrix.from_iterable(x, lambda x, y: np.abs(x-y))
        lm = ward(dm.condensed_form())
        t = TreeNode.from_linkage_matrix(lm, np.arange(len(x)).astype(np.str))

        for i, n in enumerate(t.postorder()):
            if not n.is_tip():
                n.name = "y%d" % i
            n.length = np.random.rand()*3

        md = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a', 'b', 'b'], index=index,
                      name='column-name'))

        dendrogram_heatmap(self.results, table, t, md)

        index_fp = os.path.join(self.results, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Dendrogram heatmap</h1>',
                          html)

    def test_visualization_small(self):
        # tests the scenario where ndim > number of tips
        np.random.seed(0)
        num_otus = 11  # otus
        index = pd.Index(np.arange(5).astype(np.str), name='id')
        table = pd.DataFrame(np.random.random((len(index), num_otus)),
                             index=index,
                             columns=np.arange(num_otus).astype(np.str))

        x = np.random.rand(num_otus)
        dm = DistanceMatrix.from_iterable(x, lambda x, y: np.abs(x-y))
        lm = ward(dm.condensed_form())
        t = TreeNode.from_linkage_matrix(lm, np.arange(len(x)).astype(np.str))

        for i, n in enumerate(t.postorder()):
            if not n.is_tip():
                n.name = "y%d" % i
            n.length = np.random.rand()*3

        md = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a', 'b', 'b'], index=index,
                      name='column-name'))

        dendrogram_heatmap(self.results, table, t, md)

        index_fp = os.path.join(self.results, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Dendrogram heatmap</h1>',
                          html)

    def test_visualization_garbage_metadata(self):
        # tests the scenario where ndim > number of tips
        np.random.seed(0)
        num_otus = 10  # otus
        num_samples = 5
        table = pd.DataFrame(np.random.random((num_samples, num_otus)),
                             index=np.arange(num_samples).astype(np.str),
                             columns=np.arange(num_otus).astype(np.str))

        x = np.random.rand(num_otus)
        dm = DistanceMatrix.from_iterable(x, lambda x, y: np.abs(x-y))
        lm = ward(dm.condensed_form())
        t = TreeNode.from_linkage_matrix(lm, np.arange(len(x)).astype(np.str))

        for i, n in enumerate(t.postorder()):
            if not n.is_tip():
                n.name = "y%d" % i
            n.length = np.random.rand()*3

        md = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a', 'b', 'b', 'foo', 'foo'],
                      index=pd.Index(np.arange(7).astype(np.str), name='id'),
                      name='column-name'))

        dendrogram_heatmap(self.results, table, t, md)

        index_fp = os.path.join(self.results, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Dendrogram heatmap</h1>',
                          html)

    def test_heatmap_extra_tips(self):
        # Adds in test scenario where there more tips than features
        # in the table
        np.random.seed(0)
        num_otus = 11  # otus
        index = pd.Index(np.arange(5).astype(np.str), name='id')
        table = pd.DataFrame(np.random.random((len(index), num_otus)),
                             index=index,
                             columns=np.arange(num_otus).astype(np.str))

        x = np.random.rand(num_otus*2)
        dm = DistanceMatrix.from_iterable(x, lambda x, y: np.abs(x-y))
        lm = ward(dm.condensed_form())
        t = TreeNode.from_linkage_matrix(lm, np.arange(len(x)).astype(np.str))

        for i, n in enumerate(t.postorder()):
            if not n.is_tip():
                n.name = "y%d" % i
            n.length = np.random.rand()*3

        md = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a', 'b', 'b'], index=index,
                      name='column-name'))

        dendrogram_heatmap(self.results, table, t, md)

        index_fp = os.path.join(self.results, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Dendrogram heatmap</h1>',
                          html)
