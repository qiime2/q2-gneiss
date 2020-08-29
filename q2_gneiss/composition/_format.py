import qiime2
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class CladeMetadataFormat(model.TextFileFormat):
    def validate(self, *args):
        try:
            md = qiime2.Metadata.load(str(self))
        except qiime2.metadata.MetadataFileError as md_exc:
            raise ValidationError(md_exc) from md_exc

        if md.column_count == 0:
            raise ValidationError('Format must contain at least 1 column')

        filtered_md = md.filter_columns(column_type='numeric')
        if filtered_md.column_count != md.column_count:
            raise ValidationError('Must only contain numeric values.')


CladeMetadataDirectoryFormat = model.SingleFileDirectoryFormat(
    'CladeMetadataDirectoryFormat', 'clade_metadata.tsv', CladeMetadataFormat)
