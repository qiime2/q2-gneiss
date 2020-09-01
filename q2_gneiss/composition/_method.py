# ----------------------------------------------------------------------------
# Copyright (c) 2017-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import skbio
from gneiss.composition import ilr_transform
from gneiss.util import match_tips
from gneiss.util import rename_internal_nodes

from q2_gneiss._util import add_pseudocount
from gneiss.balances import _balance_basis
from skbio import OrdinationResults


def ilr_hierarchical(table: pd.DataFrame, tree: skbio.TreeNode,
                     pseudocount: float = 0.5) -> pd.DataFrame:
    return ilr_transform(add_pseudocount(table, pseudocount), tree)


def ilr_phylogenetic(table: pd.DataFrame, tree: skbio.TreeNode,
                     pseudocount: float = 0.5) -> (
                     pd.DataFrame, skbio.TreeNode):
    t = tree.copy()
    t.bifurcate()
    table, t = match_tips(table, t)
    t = rename_internal_nodes(t)
    return ilr_transform(add_pseudocount(table, pseudocount), t), t


def ilr_phylogenetic_differential(
        differential: pd.DataFrame, tree: skbio.TreeNode) -> (
            pd.DataFrame, skbio.TreeNode):
    t = tree.copy()
    t.bifurcate()
    diff, _tree = match_tips(differential.T, t)
    _tree = rename_internal_nodes(_tree)
    in_nodes = [n.name for n in _tree.levelorder() if not n.is_tip()]
    basis = _balance_basis(_tree)[0]
    basis = pd.DataFrame(basis.T, index=diff.columns, columns=in_nodes)
    diff_balances = diff @ basis
    return diff_balances, t


def ilr_phylogenetic_ordination(table: pd.DataFrame, tree: skbio.TreeNode,
                                pseudocount: float = 0.5,
                                top_k_var: int = 10,
                                clades: str = None) -> (
                                    OrdinationResults,
                                    skbio.TreeNode, pd.DataFrame
                                ):
    t = tree.copy()
    t.bifurcate()
    _table, _tree = match_tips(table, t)
    _tree = rename_internal_nodes(_tree)
    in_nodes = [n.name for n in _tree.levelorder() if not n.is_tip()]
    basis = _balance_basis(_tree)[0]
    _table = add_pseudocount(_table, pseudocount)
    basis = pd.DataFrame(basis.T, index=_table.columns, columns=in_nodes)
    balances = np.log(_table) @ basis
    var = balances.var(axis=0).sort_values(ascending=False)
    if not clades:
        clades = var.index[:top_k_var]
    balances = balances[clades]
    balances.index.name = 'sampleid'
    # feature metadata
    basis = basis[clades]
    eigvals = var[clades]
    prop = var[clades] / var.sum()
    balances = OrdinationResults(
        short_method_name='ILR',
        long_method_name='Phylogenetic Isometric Log Ratio Transform',
        samples=balances,
        eigvals=eigvals,
        proportion_explained=prop
    )
    basis.index.name = 'featureid'
    return balances, _tree, basis
