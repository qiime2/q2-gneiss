# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import qiime2
from gneiss.plot._heatmap import heatmap
from gneiss.util import (match, match_tips)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Hierarchy
from qiime2.plugin import (Int, MetadataColumn, Categorical,
                           Str, Choices, Float)
from skbio import TreeNode
from skbio.stats.composition import clr, centralize

from q2_gneiss._util import add_pseudocount
from q2_gneiss.plugin_setup import plugin

_transform_methods = ['clr', 'log']
_mpl_colormaps = ['viridis', 'inferno', 'plasma', 'magma',
                  'Blues', 'BuGn', 'BuPu',
                  'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                  'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                  'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
                  'afmhot', 'autumn', 'bone', 'cool',
                  'copper', 'gist_heat', 'gray', 'hot',
                  'pink', 'spring', 'summer', 'winter',
                  'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                  'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                  'seismic', 'Accent', 'Dark2', 'Paired', 'Pastel1',
                  'Pastel2', 'Set1', 'Set2', 'Set3', 'Vega10',
                  'Vega20', 'Vega20b', 'Vega20c',
                  'gist_earth', 'terrain', 'ocean', 'gist_stern',
                  'brg', 'CMRmap', 'cubehelix',
                  'gnuplot', 'gnuplot2', 'gist_ncar',
                  'nipy_spectral', 'jet', 'rainbow',
                  'gist_rainbow', 'hsv', 'flag', 'prism']


# Heatmap
def dendrogram_heatmap(output_dir: str, table: pd.DataFrame,
                       tree: TreeNode,
                       metadata: qiime2.CategoricalMetadataColumn,
                       pseudocount: float = 0.5,
                       ndim: int = 10, method: str = 'clr',
                       color_map: str = 'viridis'):

    table, tree = match_tips(add_pseudocount(table, pseudocount), tree)
    nodes = [n.name for n in tree.levelorder() if not n.is_tip()]

    nlen = min(ndim, len(nodes))
    numerator_color, denominator_color = '#fb9a99', '#e31a1c'
    highlights = pd.DataFrame([[numerator_color, denominator_color]] * nlen,
                              index=nodes[:nlen])
    if method == 'clr':
        mat = pd.DataFrame(clr(centralize(table)),
                           index=table.index,
                           columns=table.columns)
    elif method == 'log':
        mat = pd.DataFrame(np.log(table),
                           index=table.index,
                           columns=table.columns)
    c = metadata.to_series()
    table, c = match(table, c)
    # TODO: There are a few hard-coded constants here
    # will need to have some adaptive defaults set in the future
    fig = heatmap(mat, tree, c, highlights, cmap=color_map,
                  highlight_width=0.01, figsize=(12, 8))
    fig.savefig(os.path.join(output_dir, 'heatmap.svg'), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'heatmap.pdf'), bbox_inches='tight')

    css = r"""
        .square {
          float: left;
          width: 100px;
          height: 20px;
          margin: 5px;
          border: 1px solid rgba(0, 0, 0, .2);
        }

        .numerator {
          background: %s;
        }

        .denominator {
          background: %s;
        }
    """ % (numerator_color, denominator_color)

    index_fp = os.path.join(output_dir, 'index.html')
    with open(index_fp, 'w') as index_f:
        index_f.write('<html><body>\n')
        index_f.write('<h1>Dendrogram heatmap</h1>\n')
        index_f.write('<img src="heatmap.svg" alt="heatmap">')
        index_f.write('<a href="heatmap.pdf">')
        index_f.write('Download as PDF</a><br>\n')
        index_f.write('<style>%s</style>' % css)
        index_f.write('<div class="square numerator">'
                      'Numerator<br/></div>')
        index_f.write('<div class="square denominator">'
                      'Denominator<br/></div>')
        index_f.write('</body></html>\n')


plugin.visualizers.register_function(
    function=dendrogram_heatmap,
    inputs={'table': FeatureTable[Frequency],
            'tree': Hierarchy},
    parameters={'metadata': MetadataColumn[Categorical],
                'ndim': Int,
                'pseudocount': Float,
                'method': Str % Choices(_transform_methods),
                'color_map': Str % Choices(_mpl_colormaps)},
    input_descriptions={
        'table': ('The feature table that will be plotted as a heatmap. '
                  'This table is assumed to have strictly positive values.'),
        'tree': ('A hierarchy of feature identifiers where each tip'
                 'corresponds to the feature identifiers in the table. '
                 'This tree can contain tip ids that are not present in '
                 'the table, but all feature ids in the table must be '
                 'present in this tree.')},
    parameter_descriptions={
        'metadata': 'Categorical metadata column to group the samples.',
        'ndim': 'Number of dimensions to highlight.',
        'pseudocount': 'The pseudocount to add to avoid division by zero.',
        'method': ("Specifies how the data should be normalized for display."
                   "Options include 'log' or 'clr' (default='clr')."),
        'color_map': ("Specifies the color map for plotting the heatmap. "
                      "See https://matplotlib.org/examples/color/"
                      "colormaps_reference.html for more details.")
    },
    name='Dendrogram heatmap.',
    description=("Visualize the feature table as a heatmap, "
                 "with samples sorted along a specified categorical metadata "
                 "column and features clustered together specified by the "
                 "tree."),
    deprecated=True,
)
