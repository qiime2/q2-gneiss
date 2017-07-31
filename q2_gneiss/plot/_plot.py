# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skbio import TreeNode
from skbio.stats.composition import clr, centralize

from q2_gneiss.plugin_setup import plugin
from gneiss.plot._heatmap import heatmap
from gneiss.plot._decompose import balance_barplots, balance_boxplot
from gneiss.util import (match, match_tips, NUMERATOR, DENOMINATOR)

from q2_types.tree import Hierarchy
from q2_gneiss.composition._type import Composition, Balance
from q2_types.feature_table import FeatureTable
from q2_types.feature_data import FeatureData, Taxonomy
from qiime2.plugin import Int, MetadataCategory, Str, Choices


def balance_taxonomy(output_dir: str, balances: pd.DataFrame, tree: TreeNode,
                     taxonomy: pd.DataFrame,
                     balance_name: Str,
                     taxa_level: Int = 0,
                     metadata: MetadataCategory = None) -> None:

    # parse out headers for taxonomy
    taxa_data = list(taxonomy['Taxon'].apply(lambda x: x.split(';')).values)
    taxa_df = pd.DataFrame(taxa_data,
                           index=taxonomy.index)

    # fill in NAs
    def f(x):
        y = np.array(list(map(lambda k: k is not None, x)))
        i = max(0, np.where(y)[0][-1])
        x[np.logical_not(y)] = [x[i]] * np.sum(np.logical_not(y))
        return x
    taxa_df = taxa_df.apply(f, axis=1)

    num_clade = tree.find(balance_name).children[NUMERATOR]
    denom_clade = tree.find(balance_name).children[DENOMINATOR]

    if num_clade.is_tip():
        num_features = pd.DataFrame(
            {num_clade.name: taxa_df.loc[num_clade.name]}
            ).T
    else:
        num_features = taxa_df.loc[num_clade.subset()]

    if denom_clade.is_tip():
        denom_features = pd.DataFrame(
            {denom_clade.name: taxa_df.loc[denom_clade.name]}
            ).T
    else:
        denom_features = taxa_df.loc[denom_clade.subset()]

    num_color, denom_color = '#4c72b0', '#4c72b0'

    fig, (ax_num, ax_denom) = plt.subplots(2)
    balance_barplots(tree, balance_name, taxa_level, taxa_df,
                     denom_color=denom_color, num_color=num_color,
                     axes=(ax_num, ax_denom))

    ax_num.set_title(
        r'$%s_{numerator} \; taxa \; (%d \; taxa)$' % (balance_name,
                                                       len(num_features)))
    ax_denom.set_title(
        r'$%s_{denominator} \; taxa \; (%d \; taxa)$' % (balance_name,
                                                         len(denom_features)))
    ax_denom.set_xlabel('Number of unique taxa')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'barplots.svg'))
    fig.savefig(os.path.join(output_dir, 'barplots.pdf'))

    if metadata is not None:
        fig2, ax = plt.subplots()
        c = metadata.to_series()
        data, c = match(balances, c)
        data[c.name] = c
        y = data[balance_name]
        # check if continuous
        try:
            c = c.astype(np.float64)
            ax.scatter(c.values, y)
            ax.set_xlabel(c.name)
        except:
            balance_boxplot(balance_name, data, y=c.name, ax=ax)

        ylabel = (r"$%s = \ln \frac{%s_{numerator}}"
                  "{%s_{denominator}}$") % (balance_name,
                                            balance_name,
                                            balance_name)
        ax.set_title(ylabel, rotation=0)
        ax.set_ylabel('log ratio')
        fig2.savefig(os.path.join(output_dir, 'balance_metadata.svg'))
        fig2.savefig(os.path.join(output_dir, 'balance_metadata.pdf'))

    index_fp = os.path.join(output_dir, 'index.html')
    with open(index_fp, 'w') as index_f:
        index_f.write('<html><body>\n')
        if metadata is not None:
            index_f.write('<h1>Balance vs %s </h1>\n' % c.name)
            index_f.write(('<img src="balance_metadata.svg" '
                           'alt="barplots">\n\n'
                           '<a href="balance_metadata.pdf">'
                           'Download as PDF</a><br>\n'))

        index_f.write(('<h1>Balance Taxonomy</h1>\n'
                       '<img src="barplots.svg" alt="barplots">\n\n'
                       '<a href="barplots.pdf">'
                       'Download as PDF</a><br>\n'
                       '<h3>Numerator taxa</h3>\n'
                       '<a href="numerator.csv">\n'
                       'Download as CSV</a><br>\n'
                       '<h3>Denominator taxa</h3>\n'
                       '<a href="denominator.csv">\n'
                       'Download as CSV</a><br>\n'))

        num_features.to_csv(os.path.join(output_dir, 'numerator.csv'),
                            header=True, index=True)
        denom_features.to_csv(os.path.join(output_dir, 'denominator.csv'),
                              header=True, index=True)
        index_f.write('</body></html>\n')


plugin.visualizers.register_function(
    function=balance_taxonomy,
    inputs={'balances': FeatureTable[Balance], 'tree': Hierarchy,
            'taxonomy': FeatureData[Taxonomy]},
    parameters={'balance_name': Str,
                'taxa_level': Int,
                'metadata': MetadataCategory},
    input_descriptions={
        'balances': 'The table of balances resulting from the ilr transform.',
        'tree': 'The tree used to calculate the balances.',
        'taxonomy': 'Taxonomy information for the OTUs.',
    },
    parameter_descriptions={
        'balance_name': 'Name of the balance to summarize.',
        'taxa_level': 'Level of taxonomy to summarize.',
        'metadata': 'Metadata for plotting the balance (optional).'},
    name='Balance Summary',
    description=("Visualize the distribution of a single balance "
                 "and summarize its numerator and denominator components.")
)


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
                       tree: TreeNode, metadata: MetadataCategory,
                       ndim=10, method='clr', color_map='viridis'):

    table, tree = match_tips(table, tree)
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

    # TODO: There are a few hard-coded constants here
    # will need to have some adaptive defaults set in the future
    fig = heatmap(mat, tree, metadata.to_series(), highlights, cmap=color_map,
                  highlight_width=0.01, figsize=(12, 8))
    fig.savefig(os.path.join(output_dir, 'heatmap.svg'))
    fig.savefig(os.path.join(output_dir, 'heatmap.pdf'))

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
    inputs={'table': FeatureTable[Composition],
            'tree': Hierarchy},
    parameters={'metadata': MetadataCategory, 'ndim': Int,
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
        'metadata': ('Metadata to group the samples. '),
        'ndim': 'Number of dimensions to highlight.',
        'method': ("Specifies how the data should be normalized for display."
                   "Options include 'log' or 'clr' (default='clr')."),
        'color_map': ("Specifies the color map for plotting the heatmap. "
                      "See https://matplotlib.org/examples/color/"
                      "colormaps_reference.html for more details.")
    },
    name='Dendrogram heatmap.',
    description=("Visualize the feature tables as a heatmap. "
                 "with samples sorted along a specified metadata category "
                 "and features clustered together specified by the tree.")
)
