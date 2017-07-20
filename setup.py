# -*- coding: utf-8 -*-
# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

from setuptools import setup, find_packages
import os

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    requirements = ["numpydoc>=0.5"]
else:
    requirements = [
        'cameo==0.10',
        'sympy>=0.7.5',
        'cobra>=0.6',
        'ipython>=4.0',
        'pyzmq>=14.5',
        'numpy>=1.9.2',
        'bokeh>=0.11.0',
        'lxml>=3.4.4',
        'six>=1.9.0',
        'ipywidgets>=4.0',
        'pandas>=0.16.2'
    ]

setup(
    name='driven',
    version="0.0.3",
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    author='Joao Cardoso',
    author_email='jooaaoo@gmail.com',
    description='driven - data-driven constraint-based analysis',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics high-throughput omics',
    url='http://driven.bio',
    long_description="A package for data-driven modeling and analysis. It implements novel and state-of-the-art methods"
                     " to integrate 'omics' data in genome-scale methods.",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'License :: OSI Approved :: Apache Software License'
    ],
)
