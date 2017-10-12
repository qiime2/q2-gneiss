# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._regression import lme_regression, ols_regression, phylogenetic_regression


__all__ = ["lme_regression", "ols_regression", "phylogenetic_regression"]
