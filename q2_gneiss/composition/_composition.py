# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import skbio
from q2_types.tree import Hierarchy, Phylogeny, Rooted
from q2_gneiss.plugin_setup import plugin
from q2_types.feature_table import (FeatureTable, Frequency, Composition,
                                    Balance)
from q2_gneiss.composition._impute import add_pseudocount
from gneiss.composition import ilr_transform
from gneiss.util import rename_internal_nodes
from qiime2.plugin import Int


def ilr_hierarchical(table: pd.DataFrame, tree: skbio.TreeNode) -> (
                     pd.DataFrame):
    return ilr_transform(table, tree)


def ilr_phylogenetic(table: pd.DataFrame, tree: skbio.TreeNode) -> (
                     pd.DataFrame, skbio.TreeNode):
    t = tree.copy()
    t.bifurcate()
    t = rename_internal_nodes(t)
    return ilr_transform(table, t), t


plugin.methods.register_function(
    function=add_pseudocount,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'pseudocount': Int},
    outputs=[('composition_table', FeatureTable[Composition])],
    input_descriptions={
        'table': 'The feature table to which pseudocounts should be added.'
    },
    parameter_descriptions={
        'pseudocount': 'The value to add to all counts in the feature table.'
    },
    output_descriptions={
        'composition_table': 'The resulting feature table.'
    },
    name='Add pseudocount to table',
    description="Increment all counts in table by pseudocount."
)


plugin.methods.register_function(
    function=ilr_hierarchical,
    inputs={'table': FeatureTable[Composition],
            'tree': Hierarchy},
    outputs=[('balances', FeatureTable[Balance])],
    parameters={},
    name='Isometric Log-ratio Transform',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'the ilr transform will be performed.'),
        'tree': ('A hierarchy of feature identifiers that defines the '
                 'partitions of features.  Each tip in the hierarchy'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.  This assumes that all of the '
                 'internal nodes in the tree have labels.')
    },
    parameter_descriptions={},
    output_descriptions={'balances': ('The resulting balances from the '
                                      'ilr transform.')},
    description="Calculate balances given a hierarchy."
)


plugin.methods.register_function(
    function=ilr_phylogenetic,
    inputs={'table': FeatureTable[Composition],
            'tree': Phylogeny[Rooted]},
    outputs=[('balances', FeatureTable[Balance]),
             ('tree', Phylogeny[Rooted])],
    parameters={},
    name='Isometric Log-ratio Transform',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'the ilr transform will be performed.'),
        'tree': ('A hierarchy of feature identifiers that defines the '
                 'partitions of features.  Each tip in the hierarchy'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.  This assumes that all of the '
                 'internal nodes in the tree have labels.')
    },
    parameter_descriptions={},
    output_descriptions={'balances': ('The resulting balances from the '
                                      'ilr transform.')},
    description="Calculate balances given a hierarchy."
)
