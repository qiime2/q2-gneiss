# ----------------------------------------------------------------------------
# Copyright (c) 2017-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.tree import Hierarchy, Phylogeny, Rooted
from q2_gneiss.plugin_setup import plugin
from q2_types.feature_table import (FeatureTable, Frequency,
                                    Balance, Composition)
from q2_types.feature_data import FeatureData, Differential
from q2_types.ordination import PCoAResults

from q2_gneiss.composition._method import (
    ilr_hierarchical, ilr_phylogenetic, ilr_phylogenetic_ordination
)
from qiime2.plugin import Float, Int, List, Str
from ._type import CladeMetadata



plugin.methods.register_function(
    function=ilr_hierarchical,
    inputs={'table': FeatureTable[Frequency | Composition],
            'tree': Hierarchy},
    outputs=[('balances', FeatureTable[Balance])],
    parameters={'pseudocount': Float},
    name='Isometric Log-ratio Transform applied to a hierarchical clustering',
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
    inputs={'table': FeatureTable[Frequency | Composition],
            'tree': Phylogeny[Rooted]},
    outputs=[('balances', FeatureTable[Balance]),
             ('hierarchy', Hierarchy)],
    parameters={'pseudocount': Float},
    name='Isometric Log-ratio Transform applied to a phylogenetic tree',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'the ilr transform will be performed.'),
        'tree': ('A rooted phylogeny of feature identifiers that defines '
                 'the partitions of features.  Each tip in the hierarchy'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.  This assumes that all of the '
                 'internal nodes in the tree have labels. This tree may '
                 'contain polytomic nodes (i.e., nodes with more than '
                 'two children), in which case they will be bifurcated.')
    },
    parameter_descriptions={
        'pseudocount': 'The value to add to zero counts in the feature table.'
    },
    output_descriptions={'balances': ('The resulting balances from the '
                                      'ilr transform.'),
                         'hierarchy': 'Hierarchy from bifurcated phylogeny'},
    description="Calculate balances given a rooted phylogeny."
)


plugin.methods.register_function(
    function=ilr_phylogenetic_ordination,
    inputs={'table': FeatureTable[Frequency | Composition],
            'tree': Phylogeny[Rooted]},
    outputs=[('ordination', PCoAResults),
             ('bifurcated_tree', Phylogeny[Rooted]),
             ('clade_metadata', FeatureData[Differential])],
    parameters={'pseudocount': Float,
                'top_k_var': Int,
                'clades': List[Str]},
    name='Ordination through a phylogenetic Isometric Log Ratio transform.',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'the ilr transform will be performed.'),
        'tree': ('A rooted phylogeny of feature identifiers that defines '
                 'the partitions of features.  Each tip in the hierarchy'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.  This assumes that all of the '
                 'internal nodes in the tree have labels. This tree may '
                 'contain polytomic nodes (i.e., nodes with more than '
                 'two children), in which case they will be bifurcated.')
    },
    parameter_descriptions={
        'pseudocount': 'The value to add to zero counts in the feature table.',
        'top_k_var': 'The top k most variable balances.',
        'clades': 'The names of clades to focus on (overrides top-k-var).',
    },
    output_descriptions={
        'ordination': ('The resulting ordination from the '
                       'ilr transform.'),
        'bifurcated_tree': 'Bifurcating phylogeny',
        'clade_metadata': ('Metadata specifying clade membership.')
    },
    description="Compute an ILR ordination given a rooted phylogeny."
)
