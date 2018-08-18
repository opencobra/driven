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

"""Test functionalities of FluxConstraints."""

from __future__ import absolute_import

import numpy as np
import pandas as pd

from driven.data_sets import FluxConstraints


def test_indexing(mock_flux_constraints):
    """Test indexing of FluxConstraints."""
    assert (mock_flux_constraints["R2"] == np.array([0.5, 0.7])).all()


def test_equality(mock_flux_constraints):
    """Test equality check of FluxConstraints."""
    eq_flux_constraints = FluxConstraints(limits=np.array([[0, 10],
                                                           [0.5, 0.7],
                                                           [5.1, 5.2]]),
                                          reaction_ids=["R1", "R2", "R3"])
    assert mock_flux_constraints == eq_flux_constraints
    del eq_flux_constraints
    assert mock_flux_constraints != ["T", "E", "S", "T"]


def test_from_data_frame(mock_flux_constraints):
    """Test import from pandas.DataFrame."""
    data_frame = pd.DataFrame(np.array([[0, 10],
                                        [0.5, 0.7],
                                        [5.1, 5.2]]),
                              index=["R1", "R2", "R3"],
                              columns=["lower_limit", "upper_limit"])
    converted_data_frame = FluxConstraints.from_data_frame(data_frame,
                                                           data_type="constraints")
    del data_frame
    assert mock_flux_constraints == converted_data_frame


def test_from_csv(mock_flux_constraints, flux_constraints_csv):
    """Test import from CSV file."""
    csv_data = FluxConstraints.from_csv(flux_constraints_csv,
                                        data_type="constraints",
                                        index_col=0)
    assert mock_flux_constraints == csv_data
