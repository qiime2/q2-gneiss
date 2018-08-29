# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-gneiss",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Jamie Morton",
    author_email="jamietmorton@gmail.com",
    description="Compositional Data Analysis and Visualization Toolbox",
    license='BSD-3-Clause',
    url="https://github.com/qiime2/q2-gneiss",
    entry_points={
        'qiime2.plugins': ['q2-gneiss=q2_gneiss.plugin_setup:plugin']
    },
    package_data={
        "q2_gneiss": ['citations.bib'],
        'q2_gneiss.cluster.tests': ['data/*'],
        'q2_gneiss.regression.tests': ['data/*'],
    },
    zip_safe=False,
)
