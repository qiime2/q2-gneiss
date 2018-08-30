# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

import qiime2.plugin
import qiime2.sdk
from q2_gneiss import __version__


citations = qiime2.plugin.Citations.load('citations.bib', package='q2_gneiss')

plugin = qiime2.plugin.Plugin(
    name='gneiss',
    version=__version__,
    website='https://biocore.github.io/gneiss/',
    citations=[citations['morton2017balance']],
    short_description=('Plugin for building compositional models.'),
    description=('This is a QIIME 2 plugin supporting statistical models on '
                 'feature tables and metadata using balances.'),
    package='q2_gneiss')

importlib.import_module('q2_gneiss.composition')
importlib.import_module('q2_gneiss.regression')
importlib.import_module('q2_gneiss.plot')
importlib.import_module('q2_gneiss.cluster')
