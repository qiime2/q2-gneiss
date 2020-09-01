# ----------------------------------------------------------------------------
# Copyright (c) 2017-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._format import CladeMetadataFormat
from ..plugin_setup import plugin
import pandas as pd
import qiime2


# differential types
@plugin.register_transformer
def _322(ff: CladeMetadataFormat) -> pd.DataFrame:
    return qiime2.Metadata.load(str(ff)).to_dataframe()


@plugin.register_transformer
def _323(ff: CladeMetadataFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _324(data: pd.DataFrame) -> CladeMetadataFormat:
    ff = CladeMetadataFormat()
    qiime2.Metadata(data).save(str(ff))
    return ff