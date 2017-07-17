# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_composition.plugin_setup import Composition
from q2_types.tree import Hierarchy
from q2_gneiss.plugin_setup import plugin
from q2_composition._type import Balance
from q2_types.feature_table import FeatureTable
from gneiss.composition import ilr_transform


plugin.methods.register_function(
    function=ilr_transform,
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
