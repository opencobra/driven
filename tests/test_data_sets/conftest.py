# -*- coding: utf-8 -*-

# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University Denmark
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

"""Define module level fixtures."""

from __future__ import absolute_import

import numpy as np
import pytest

from driven.data_sets import ExpressionProfile
from driven.data_sets import FluxConstraints


@pytest.fixture(scope="session")
def mock_expression_profile():
    """Generate mock ExpressionProfile."""
    genes = ["G1", "G2"]
    time_series = ["T1", "T2"]
    expression = np.array([[10, 11],
                           [65, 109]])
    pvalues = np.array([[0.02],
                        [0.048]])
    return ExpressionProfile(expression=expression,
                             identifiers=genes,
                             conditions=time_series,
                             p_values=pvalues)


@pytest.fixture(scope="session")
def mock_flux_constraints():
    """Generate mock FluxConstraints."""
    reaction_ids = ["R1", "R2", "R3"]
    limits = np.array([[0, 10],
                       [0.5, 0.7],
                       [5.1, 5.2]])
    return FluxConstraints(reaction_ids, limits)
