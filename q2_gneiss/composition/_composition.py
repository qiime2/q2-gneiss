# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.tree import Hierarchy, Phylogeny, Rooted
from q2_gneiss.plugin_setup import plugin
from q2_types.feature_table import (FeatureTable, Frequency, Composition,
                                    Balance)
from q2_gneiss.composition._method import ilr_hierarchical, ilr_phylogenetic
from qiime2.plugin import Float


plugin.methods.register_function(
    function=ilr_hierarchical,
    inputs={'table': FeatureTable[Frequency],
            'tree': Hierarchy},
    outputs=[('balances', FeatureTable[Balance])],
    parameters={'pseudocount': Float},
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
    parameter_descriptions={
        'pseudocount': 'The value to add to zero counts in the feature table.'
    },
    output_descriptions={'balances': ('The resulting balances from the '
                                      'ilr transform.')},
    description="Calculate balances given a hierarchy."
)


plugin.methods.register_function(
    function=ilr_phylogenetic,
    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted]},
    outputs=[('balances', FeatureTable[Balance]),
             ('tree', Phylogeny[Rooted])],
    parameters={'pseudocount': Float},
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
    parameter_descriptions={
        'pseudocount': 'The value to add to zero counts in the feature table.'
    },
    output_descriptions={'balances': ('The resulting balances from the '
                                      'ilr transform.')},
    description="Calculate balances given a hierarchy."
)
