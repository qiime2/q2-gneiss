# ----------------------------------------------------------------------------
# Copyright (c) 2017-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import os
import shutil

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
from scipy.cluster.hierarchy import ward

from skbio import TreeNode, DistanceMatrix
from skbio.stats.composition import ilr_inv
from gneiss.balances import balance_basis
from q2_gneiss.plot._plot import dendrogram_heatmap, balance_taxonomy
from qiime2 import CategoricalMetadataColumn, NumericMetadataColumn


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


class TestBalanceTaxonomy(unittest.TestCase):

    def setUp(self):
        self.results = "results"
        if not os.path.exists(self.results):
            os.mkdir(self.results)
        self.balances = pd.DataFrame(
            {'a': [-2, -1, 0, 1, 2],
             'b': [-2, 0, 0, 0, 0]},
            index=['a1', 'a2', 'a3', 'a4', 'a5']
        )
        self.tree = TreeNode.read([r'((k, q)d, ((x, y)a, z)b)c;'])

        self.taxonomy = pd.DataFrame(
            [['foo;barf;a;b;c;d;e', 1],
             ['foo;bark;f;g;h;i;j', 1],
             ['foo;bark;f;g;h;w;j', 1],
             ['nom;tu;k;l;m;n;o', 0.9],
             ['nom;tu;k;l;m;t;o', 0.9]],
            columns=['Taxon', 'Confidence'],
            index=['x', 'y', 'z', 'k', 'q'])

        self.balances = pd.DataFrame(
            [[1, 2, 3, 4, 5, 6, 7],
             [-3.1, -2.9, -3, 3, 2.9, 3.2, 3.1],
             [1, 1, 1, 1, 1, 1, 1],
             [3, 2, 1, 0, -1, -2, -3]],
            index=['d', 'a', 'b', 'c'],
            columns=['s1', 's2', 's3', 's4', 's5', 's6', 's7']
        ).T
        basis, _ = balance_basis(self.tree)
        self.table = pd.DataFrame(
            ilr_inv(self.balances, basis),
            columns=['x', 'y', 'z', 'k', 'q'],
            index=['s1', 's2', 's3', 's4', 's5', 's6', 's7']
        )

        index = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6', 's7'], name='id')
        self.categorical = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a', 'b', 'b', 'b', 'b'],
                      index=index, name='categorical'))
        self.multi_categorical = CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'c', 'b', 'b', 'b', 'c'],
                      index=index, name='multi_categorical'))
        self.partial_numerical_categorical = CategoricalMetadataColumn(
            pd.Series(['1', '1', '1', '2', '2', '2', 'a'],
                      index=index, name='multi_categorical'))
        self.full_numerical_categorical = CategoricalMetadataColumn(
            pd.Series(['1', '1', '1.0', '2', '2', '2.0', '3'],
                      index=index, name='numerical_categorical'))
        self.continuous = NumericMetadataColumn(
            pd.Series(np.arange(7), index=index, name='continuous'))

    def tearDown(self):
        shutil.rmtree(self.results)

    def test_balance_taxonomy(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='c')
        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Balance Taxonomy</h1>', html)
            self.assertIn('Numerator taxa', html)
            self.assertIn('Denominator taxa', html)

        # extract csv files and test for contents
        exp = pd.DataFrame(
            [['foo', 'barf', 'a', 'b', 'c', 'd', 'e'],
             ['foo', 'bark', 'f', 'g', 'h', 'i', 'j'],
             ['foo', 'bark', 'f', 'g', 'h', 'w', 'j']],
            columns=['0', '1', '2', '3', '4', '5', '6'],
            index=['x', 'y', 'z'])
        res = pd.read_csv(num_fp, index_col=0)
        pdt.assert_frame_equal(exp, res.sort_index())

        exp = pd.DataFrame([['nom', 'tu', 'k', 'l', 'm', 't', 'o'],
                            ['nom', 'tu', 'k', 'l', 'm', 'n', 'o']],
                           columns=['0', '1', '2', '3', '4', '5', '6'],
                           index=['q', 'k']).sort_index()
        res = pd.read_csv(denom_fp, index_col=0)
        pdt.assert_frame_equal(exp, res.sort_index())

    def test_balance_taxonomy_tips(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='a')
        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))

        exp = pd.DataFrame(['foo', 'bark', 'f', 'g', 'h', 'i', 'j'],
                           index=['0', '1', '2', '3', '4',
                                  '5', '6'],
                           columns=['y']).T
        res = pd.read_csv(num_fp, index_col=0)
        pdt.assert_frame_equal(exp, res)

        res = pd.read_csv(denom_fp, index_col=0)
        exp = pd.DataFrame(['foo', 'barf', 'a', 'b', 'c', 'd', 'e'],
                           index=['0', '1', '2', '3', '4', '5', '6'],
                           columns=['x']).T
        pdt.assert_frame_equal(exp, res)

    def test_balance_taxonomy_categorical(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='a',
                         metadata=self.categorical)
        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))
        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))

        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Balance Taxonomy</h1>', html)
            self.assertIn('Numerator taxa', html)
            self.assertIn('Denominator taxa', html)
            self.assertIn('Proportion', html)

    def test_balance_taxonomy_categorical_error(self):
        with self.assertRaises(ValueError):
            balance_taxonomy(self.results, self.table, self.tree,
                             self.taxonomy, balance_name='a',
                             metadata=self.categorical,
                             threshold=100.)

    def test_balance_taxonomy_multi_categorical(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='a',
                         metadata=self.multi_categorical)
        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))
        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))
        prop_fp = os.path.join(self.results, 'proportion_plot.pdf')
        self.assertFalse(os.path.exists(prop_fp))

        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))
        prop_fp = os.path.join(self.results, 'proportion_plot.pdf')
        self.assertFalse(os.path.exists(prop_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Balance Taxonomy</h1>', html)
            self.assertIn('Numerator taxa', html)
            self.assertIn('Denominator taxa', html)
            self.assertNotIn('Proportion', html)

    def test_balance_taxonomy_partial_numerical_categorical(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='a',
                         metadata=self.partial_numerical_categorical)

        for file in ['numerator.csv', 'denominator.csv',
                     'balance_metadata.pdf']:
            self.assertTrue(os.path.exists(os.path.join(self.results, file)))
        self.assertFalse(os.path.exists(os.path.join(self.results,
                                                     'proportion_plot.pdf')))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            for exp in ['<h1>Balance Taxonomy</h1>', 'Numerator taxa',
                        'Denominator taxa']:
                self.assertIn(exp, html)
            self.assertNotIn('Proportion', html)

    def test_balance_taxonomy_full_numerical_categorical(self):
        with self.assertRaisesRegex(ValueError, 'only numerical values'):
            balance_taxonomy(self.results, self.table, self.tree,
                             self.taxonomy, balance_name='a',
                             metadata=self.full_numerical_categorical)

    def test_balance_taxonomy_continuous(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='a',
                         metadata=self.continuous)

        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))
        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))
        prop_fp = os.path.join(self.results, 'proportion_plot.pdf')
        self.assertTrue(os.path.exists(prop_fp))

        box_fp = os.path.join(self.results, 'balance_metadata.pdf')
        self.assertTrue(os.path.exists(box_fp))

        with open(index_fp, 'r') as fh:
            html = fh.read()
            self.assertIn('<h1>Balance Taxonomy</h1>', html)
            self.assertIn('Numerator taxa', html)
            self.assertIn('Denominator taxa', html)
            self.assertIn('Proportion', html)

    def test_balance_taxonomy_genus(self):
        index_fp = os.path.join(self.results, 'index.html')
        balance_taxonomy(self.results, self.table, self.tree,
                         self.taxonomy, balance_name='c',
                         taxa_level=6)
        self.assertTrue(os.path.exists(index_fp))
        # test to make sure that the numerator file is there
        num_fp = os.path.join(self.results, 'numerator.csv')
        self.assertTrue(os.path.exists(num_fp))
        # test to make sure that the denominator file is there
        denom_fp = os.path.join(self.results, 'denominator.csv')
        self.assertTrue(os.path.exists(denom_fp))


if __name__ == "__main__":
    unittest.main()
