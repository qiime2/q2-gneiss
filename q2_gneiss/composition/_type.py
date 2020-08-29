# ----------------------------------------------------------------------------
# Copyright (c) 2017-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin import SemanticType
from ..plugin_setup import plugin
from q2_types.feature_data import FeatureData
from ._format import CladeMetadataDirectoryFormat


CladeMetadata = SemanticType('CladeMetadata',
                             variant_of=FeatureData.field['type'])

plugin.register_semantic_types(CladeMetadata)


plugin.register_semantic_type_to_format(
    FeatureData[CladeMetadata], CladeMetadataDirectoryFormat)
