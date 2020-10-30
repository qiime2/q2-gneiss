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
from gneiss.util import NUMERATOR, DENOMINATOR

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
    diff_balances = (diff @ basis).T
    diff_balances.index.name = 'featureid'
    return diff_balances, t


def get_children(tree, side):
    if tree.children[side].is_tip():
        return [tree.children[side].name]
    else:
        return [n.name for n in tree.children[side].tips()]


def logmean(table, tips, pseudocount):
    if len(tips) == 1:
        return np.array(np.log(table[tips] + pseudocount)).ravel()
    else:
        return np.array(np.log(table[tips] + pseudocount).mean(axis=1))


def _fast_ilr(tree, table, clades, pseudocount=0.5):
    # manually computes the ILR transform on a subset of specified clades
    balances = {}
    basis = {}
    for c in clades:
        subtree = tree.find(c)
        num_tips = get_children(subtree, NUMERATOR)
        denom_tips = get_children(subtree, DENOMINATOR)
        r, s = len(num_tips),  len(denom_tips)
        if r == 0 or s == 0:
            raise ValueError(f'Clade {c} has no children {subtree}')
        Z = np.sqrt(r * s / (r + s))

        num = logmean(table, num_tips, pseudocount)
        denom = logmean(table, denom_tips, pseudocount)
        print(num.shape, denom.shape)
        balances[c] = pd.Series(Z * (num - denom), index=table.index)
        basis[c] = pd.Series(np.zeros(len(table.columns)), index=table.columns)
        basis[c].loc[num_tips] = np.sqrt(s / (r * (r + s)))
        basis[c].loc[denom_tips] = -np.sqrt(s / (s * (r + s)))
    balances = pd.DataFrame(balances)
    basis = pd.DataFrame(basis)
    return balances, basis


def ilr_phylogenetic_ordination(table: pd.DataFrame, tree: skbio.TreeNode,
                                pseudocount: float = 0.5,
                                top_k_var: int = 10,
                                clades: list = None) -> (
                                    OrdinationResults,
                                    skbio.TreeNode, pd.DataFrame
                                ):
    t = tree.copy()
    t.bifurcate()
    _table, _tree = match_tips(table, t)
    _tree = rename_internal_nodes(_tree)
    if not clades:
        clades = var.index[:top_k_var]
        in_nodes = [n.name for n in _tree.levelorder() if not n.is_tip()]
        basis = _balance_basis(_tree)[0]
        _table = add_pseudocount(_table, pseudocount)
        basis = pd.DataFrame(basis.T, index=_table.columns, columns=in_nodes)
        balances = np.log(_table) @ basis
        balances = balances[clades]
        basis = basis[clades]
    else:
        clades = clades[0].split(',')
        balances, basis = _fast_ilr(_tree, _table, clades, pseudocount=0.5)

    var = balances.var(axis=0).sort_values(ascending=False)
    balances.index.name = 'sampleid'
    # feature metadata
    eigvals = var
    prop = var[clades] / var.sum()
    balances = OrdinationResults(
        short_method_name='ILR',
        long_method_name='Phylogenetic Isometric Log Ratio Transform',
        samples=balances,
        features=pd.DataFrame(np.eye(len(clades)), index=clades),
        eigvals=eigvals,
        proportion_explained=prop
    )
    basis.index.name = 'featureid'
    return balances, _tree, basis
