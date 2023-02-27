# ----------------------------------------------------------------------------
# Copyright (c) 2017-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


# These updates allow `gradient_linkage` to work with scipy=1.9.0 which
# prevent euclidean from accepting scalar inputs.
# Instead a lambda which is equivalent to euclidean in the 0d case is used.

# The code below this point originated from gneiss, and is subject to the
# following copyright and terms:

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, gneiss development team.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the names scikit-bio, skbio, or biocore nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

import numpy as np
from scipy.cluster.hierarchy import linkage
from skbio import DistanceMatrix, TreeNode

from gneiss.sort import mean_niche_estimator
from gneiss.util import match, rename_internal_nodes


def gradient_linkage(X, y, method='average'):
    # Taken from https://github.com/biocore/gneiss/blob
    # /5d253d68ef14fa82e26b5f74118ff170f1990585/gneiss/cluster/_pba.py#L132
    _X, _y = match(X, y)
    mean_X = mean_niche_estimator(_X, gradient=_y)
    t = _rank_linkage(mean_X)
    return t


def _rank_linkage(r, method='average'):
    # modified from https://github.com/biocore/gneiss/blob
    # /5d253d68ef14fa82e26b5f74118ff170f1990585/gneiss/cluster/_pba.py#L82
    # START MODIFICATION
    def euclidean_1d(a, b): return np.abs(b-a)
    dm = DistanceMatrix.from_iterable(r, euclidean_1d)
    # END MODIFICATION

    lm = linkage(dm.condensed_form(), method)
    t = TreeNode.from_linkage_matrix(lm, r.index)
    t = rename_internal_nodes(t)
    return t
