# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import skbio
from gneiss.composition import ilr_transform
from gneiss.util import rename_internal_nodes


def add_pseudocount(table: pd.DataFrame, pseudocount: float=0.5) -> (
                    pd.DataFrame):
    return table.replace(0, pseudocount)


def ilr_hierarchical(table: pd.DataFrame, tree: skbio.TreeNode,
                     pseudocount: float=0.5) -> pd.DataFrame:
    return ilr_transform(add_pseudocount(table, pseudocount), tree)


def ilr_phylogenetic(table: pd.DataFrame, tree: skbio.TreeNode,
                     pseudocount: float=0.5) -> (
                     pd.DataFrame, skbio.TreeNode):
    t = tree.copy()
    t.bifurcate()
    t = rename_internal_nodes(t)
    return ilr_transform(add_pseudocount(table, pseudocount), t), t
