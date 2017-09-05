# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import skbio

from q2_types.feature_table import FeatureTable, Balance, Frequency
from q2_types.tree import Hierarchy
from qiime2.plugin import Str, Metadata, Float
from q2_gneiss.plugin_setup import plugin
from gneiss.plot._regression_plot import ols_summary, lme_summary
from gneiss.util import match_tips, match
from gneiss.composition import ilr_transform
from gneiss.regression._ols import ols
from gneiss.regression._mixedlm import mixedlm


def ols_regression(output_dir: str,
                   table: pd.DataFrame, tree: skbio.TreeNode,
                   metadata: Metadata, formula: str) -> None:
    res = ols(table=table, metadata=metadata._dataframe,
              formula=formula)
    res.fit()

    ols_summary(output_dir, res, tree)

plugin.visualizers.register_function(
    function=ols_regression,
    inputs={'table': FeatureTable[Balance],
            'tree': Hierarchy},
    parameters={'formula': Str, 'metadata': Metadata},
    name='Simplicial Ordinary Least Squares Regression',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'simplicial regression will be performed.'),
        'tree': ('A hierarchy of feature identifiers where each tip'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.')
    },
    parameter_descriptions={
        'formula': 'Statistical formula specifying the statistical model.',
        'metadata': ('Metadata information that contains the '
                     'covariates of interest.')
    },
    description="Perform linear regression on balances."
)


def core_regression(output_dir: str,
                    table: pd.DataFrame, tree: skbio.TreeNode,
                    metadata: Metadata, formula: str,
                    pseudocount: float=1):

    _table, _tree = match_tips(table, tree)
    balances = ilr_transform(_table+pseudocount, _tree)
    ols_regression(output_dir, balances, _tree,
                   metadata, formula)


plugin.visualizers.register_function(
    function=core_regression,
    inputs={'table': FeatureTable[Frequency],
            'tree': Hierarchy},
    parameters={'formula': Str, 'metadata': Metadata,
                'pseudocount': Float},
    name='Simplicial Ordinary Least Squares Regression',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'simplicial regression will be performed.'),
        'tree': ('A hierarchy of feature identifiers where each tip'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.'),
    },
    parameter_descriptions={
        'formula': 'Statistical formula specifying the statistical model.',
        'metadata': ('Metadata information that contains the '
                     'covariates of interest.'),
        'pseudocount': ('The pseudocount to add to avoid taking the '
                        'logarithm of zero.')
    },
    description="Perform simplicial linear regression."
)


def lme_regression(output_dir: str,
                   table: pd.DataFrame, tree: skbio.TreeNode,
                   metadata: Metadata, formula: str,
                   groups: str) -> None:
    res = mixedlm(table=table, metadata=metadata._dataframe,
                  formula=formula, groups=groups)
    res.fit()
    lme_summary(output_dir, res, tree)


plugin.visualizers.register_function(
    function=lme_regression,
    inputs={'table': FeatureTable[Balance],
            'tree': Hierarchy},
    parameters={'metadata': Metadata, 'formula': Str, 'groups': Str},
    name='Simplicial Linear mixed effects regression',
    input_descriptions={
        'table': ('The feature table containing the samples in which '
                  'simplicial regression with mixed effects will be performed'
                  'will be performed.'),
        'tree': ('A hierarchy of feature identifiers where each tip'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.')
    },
    parameter_descriptions={
        'formula': 'Statistical formula specifying the statistical model.',
        'metadata': ('Metadata information that contains the '
                     'covariates of interest.')
    },
    description="Build and run linear mixed effects model on balances."
)
