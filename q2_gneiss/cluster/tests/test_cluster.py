# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import qiime2
import pandas as pd
from skbio.util import get_data_path
from skbio import TreeNode
import pandas.util.testing as pdt
import numpy.testing as npt


class TestClusteringPlugin(unittest.TestCase):

    def assert_tree_almost_equals(self, a, b):
        for n, m in zip(a.postorder(include_self=True),
                        b.postorder(include_self=True)):
            self.assertEqual(n.name, m.name)
            if n.length is None or m.length is None:
                self.assertEqual(n.length, m.length)
            else:
                npt.assert_almost_equal(n.length, m.length)

    def test_proportional_artifact(self):
        from qiime2.plugins.gneiss.methods import correlation_clustering
        table_f = get_data_path("feature-table.qza")
        in_table = qiime2.Artifact.load(table_f)

        res = correlation_clustering(in_table, pseudocount=0.1)
        res_clust = res.clustering._view(TreeNode)
        exp_str = ('((F4:0.228723591874,(F5:0.074748541601,'
                   '(F1:0.00010428164962,F2:0.00010428164962)'
                   'y4:0.0746442599513)y3:0.153975050273)'
                   'y1:0.70266138894,(F3:0.266841737789,F6:0.266841737789)'
                   'y2:0.664543243026)y0;\n')
        exp_tree = TreeNode.read([exp_str])
        self.assert_tree_almost_equals(exp_tree, res_clust)

    def test_gradient_artifact(self):
        from qiime2.plugins.gneiss.methods import gradient_clustering
        table_f = get_data_path("test_gradient.biom.qza")
        metadata_f = get_data_path("test_metadata.txt")
        in_table = qiime2.Artifact.load(table_f)
        in_metadata = qiime2.Metadata.load(metadata_f)

        res = gradient_clustering(in_table, in_metadata.get_column('x'))
        res_clust = res.clustering._view(TreeNode)
        exp_str = '((o1:0.5,o2:0.5)y1:0.5,(o3:0.5,o4:0.5)y2:0.5)y0;\n'
        self.assertEqual(exp_str, str(res_clust))

    def test_gradient_match(self):
        # there are extra rows in match that need to be filtered out
        from qiime2.plugins.gneiss.methods import gradient_clustering
        table_f = get_data_path("test_gradient.biom.qza")
        metadata_f = get_data_path("test_metadata2.txt")
        in_table = qiime2.Artifact.load(table_f)
        in_metadata = qiime2.Metadata.load(metadata_f)

        res = gradient_clustering(in_table, in_metadata.get_column('x'))
        res_clust = res.clustering._view(TreeNode)
        exp_str = '((o1:0.5,o2:0.5)y1:0.5,(o3:0.5,o4:0.5)y2:0.5)y0;\n'
        self.assertEqual(exp_str, str(res_clust))

    def test_gradient_artifact_weighted(self):
        from qiime2.plugins.gneiss.methods import gradient_clustering
        table_f = get_data_path("weighted.biom.qza")
        metadata_f = get_data_path("test_metadata.txt")
        in_table = qiime2.Artifact.load(table_f)
        in_metadata = qiime2.Metadata.load(metadata_f)

        res_uw = gradient_clustering(in_table, in_metadata.get_column('x'),
                                     weighted=False)
        res_w = gradient_clustering(in_table, in_metadata.get_column('x'),
                                    weighted=True)
        res_clust_uw = res_uw.clustering._view(TreeNode)
        res_clust_w = res_w.clustering._view(TreeNode)

        self.assertNotEqual(str(res_clust_uw), str(res_clust_w))

    def test_gradient_missing_samples(self):
        from qiime2.plugins.gneiss.methods import gradient_clustering
        table = pd.DataFrame({"x": 1, "y": 2}, index=["a", "s1"])
        table = qiime2.Artifact.import_data("FeatureTable[Frequency]", table)
        metadata = qiime2.Metadata.load(get_data_path("test_metadata.txt"))

        with self.assertRaisesRegex(KeyError, "not present.*a"):
            gradient_clustering(table, metadata.get_column("x"))

    def test_gradient_ignore_missing_samples(self):
        from qiime2.plugins.gneiss.methods import gradient_clustering
        table = pd.DataFrame({"x": 1, "y": 2}, index=["a", "s1"])
        table = qiime2.Artifact.import_data("FeatureTable[Frequency]", table)
        metadata = qiime2.Metadata.load(get_data_path("test_metadata.txt"))

        gradient_clustering(table, metadata.get_column("x"),
                            ignore_missing_samples=True)
        # Checkpoint assertion
        self.assertTrue(True)

    def test_assign_ids(self):
        from qiime2.plugins.gneiss.methods import assign_ids
        tree_f = get_data_path("tree.qza")
        table_f = get_data_path("tree_table.qza")
        tree = qiime2.Artifact.load(tree_f)
        table = qiime2.Artifact.load(table_f)
        out_tree = assign_ids(input_tree=tree, input_table=table)
        res_t = out_tree.output_tree._view(TreeNode)
        for n in res_t.levelorder(include_self=True):
            self.assertTrue(n.name is not None)

    def test_assign_ids_intersect(self):
        from qiime2.plugins.gneiss.methods import assign_ids
        tree_f = get_data_path("tree_extra.qza")
        table_f = get_data_path("polytomy_table.qza")
        tree = qiime2.Artifact.load(tree_f)
        table = qiime2.Artifact.load(table_f)
        output = assign_ids(input_tree=tree, input_table=table)
        res_tree = output.output_tree._view(TreeNode)
        res_table = output.output_table._view(pd.DataFrame)
        for n in res_tree.levelorder(include_self=True):
            self.assertTrue(n.name is not None)
        exp = list('abde')
        res = [n.name for n in res_tree.tips()]
        self.assertEqual(exp, res)

        exp = pd.DataFrame({
            's1': [1.0, 2.0, 4.0, 5.0],
            's2': [1.0, 5.0, 6.0, 0.0]
        }, index=['a', 'b', 'd', 'e']).T

        pdt.assert_frame_equal(exp, res_table)

    def test_assign_ids_polytomy(self):
        from qiime2.plugins.gneiss.methods import assign_ids
        tree_f = get_data_path("polytomy.qza")
        table_f = get_data_path("polytomy_table.qza")
        tree = qiime2.Artifact.load(tree_f)
        table = qiime2.Artifact.load(table_f)
        out_tree = assign_ids(input_tree=tree, input_table=table)
        res_t = out_tree.output_tree._view(TreeNode)
        res_nontips = []
        for n in res_t.levelorder(include_self=True):
            self.assertTrue(n.name is not None)
            if not n.is_tip():
                res_nontips.append(n.name)
        self.assertEqual(len(res_nontips), 4)


if __name__ == '__main__':
    unittest.main()
