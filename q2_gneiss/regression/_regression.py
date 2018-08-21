# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import skbio
from gneiss.regression._ols import ols
from gneiss.regression._mixedlm import mixedlm

from q2_types.feature_table import FeatureTable, Balance
from q2_types.tree import Hierarchy
from qiime2.plugin import Str, Metadata
from q2_gneiss.plugin_setup import plugin
from gneiss.plot._regression_plot import ols_summary, lme_summary
import numpy as np


def ols_regression(output_dir: str,
                   table: pd.DataFrame, tree: skbio.TreeNode,
                   metadata: Metadata, formula: str) -> None:

    if np.any(table.var(axis=0) == 0):
        message = ('Detected zero variance balances - '
                   'double check your table for unobserved features.')
        raise UserWarning(message)

    res = ols(table=table, metadata=metadata.to_dataframe(),
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
        'formula': 'Formula specifying the statistical model. '
                   'In other words, a list of the metadata categories that '
                   'will be used in the regression model, sepearted by "+" '
                   'For more information see https://patsy.readthedocs.io/en/latest/API-reference.html',
        'metadata': ('Metadata information that contains the '
                     'covariates of interest.')
    },
    description="Perform linear regression on balances. This will tell how much variability is explained by metadata categories in your formula."
)


def lme_regression(output_dir: str,
                   table: pd.DataFrame, tree: skbio.TreeNode,
                   metadata: Metadata, formula: str,
                   groups: str) -> None:
    if np.any(table.var(axis=0) == 0):
        message = ('Detected zero variance balances - '
                   'double check your table for unobserved features.')
        raise UserWarning(message)

    res = mixedlm(table=table, metadata=metadata.to_dataframe(),
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
        'formula': 'Statistical formula specifying the statistical model.'
                   'In other words, a list of the metadata categories that '
                   'will be used in the linear mixed effect model, sepearted by "+" '
                   'For more information see https://patsy.readthedocs.io/en/latest/API-reference.html',
        'metadata': ('Metadata information that contains the '
                     'covariates of interest.')
    },
    description="Build and run linear mixed effects model on balances."
)
