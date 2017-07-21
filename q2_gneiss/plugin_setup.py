# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

import qiime2.plugin
import qiime2.sdk
from q2_gneiss import __version__


_citation = ('Morton JT, Sanders J, Quinn RA, McDonald D, Gonzalez A, '
             'Vázquez-Baeza Y Navas-Molina JA, Song SJ, Metcalf JL, '
             'Hyde ER, Lladser M, Dorrestein PC, Knight R. 2017. '
             'Balance trees reveal microbial niche differentiation '
             'mSystems 2:e00162-16. '
             'https://doi.org/10.1128/mSystems.00162-16.')

plugin = qiime2.plugin.Plugin(
    name='gneiss',
    version=__version__,
    website='https://biocore.github.io/gneiss/',
    citation_text=_citation,
    short_description=('Plugin for building compositional models.'),
    description=('This is a QIIME 2 plugin supporting statistical models on '
                 'feature tables and metadata using balances.'),
    package='q2_gneiss')

importlib.import_module('q2_gneiss.composition')
importlib.import_module('q2_gneiss.regression')
importlib.import_module('q2_gneiss.plot')
importlib.import_module('q2_gneiss.cluster')
