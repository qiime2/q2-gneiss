# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages
import re

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')


setup(
    name="q2-gneiss",
    version='0.0.1',
    packages=find_packages(),
    # pandas and q2-dummy-types are only required for the dummy methods and
    # visualizers provided as examples. Remove these dependencies when you're
    # ready to develop your plugin, and add your own dependencies (if there are
    # any).
    install_requires=['pandas', 'numpy', 'gneiss>=0.3.0'],
    author="Jamie Morton",
    author_email="jamietmorton@gmail.com",
    description="Compositional Data Analysis and Visualization Toolbox",
    license='BSD-3-Clause',
    url="https://github.com/biocore/gneiss",
    extra_require={
        'q2': ['qiime2 >= 2017.2.0', 'biom-format', 'seaborn']
    },
    entry_points={
        'qiime2.plugins': ['q2-gneiss=q2_gneiss.plugin_setup:plugin']
    }
)
