import pandas as pd
import qiime2
import . import CladeMetadataFormat


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
