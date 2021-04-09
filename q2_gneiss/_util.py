# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def add_pseudocount(table: pd.DataFrame, pseudocount: float = 0.5) -> (
                    pd.DataFrame):
    return table.replace(0, pseudocount)
