# -*- coding: utf-8 -*-

# Copyright (c) 2015 Novo Nordisk Foundation Center for Biosustainability,
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

"""Test functionalities of ExpressionProfile."""

from __future__ import absolute_import

import numpy as np
import pandas as pd

from driven.data_sets import ExpressionProfile


def test_indexing(mock_expression_profile):
    """Test indexing of ExpressionProfile."""
    assert mock_expression_profile["G1", "T1"] == 10
    assert mock_expression_profile[1, 1] == 109


def test_equality(mock_expression_profile):
    """Test equality check of ExpressionProfile."""
    eq_expression_profile = ExpressionProfile(expression=np.array([[10, 11],
                                                                   [65, 109]]),
                                              identifiers=["G1", "G2"],
                                              conditions=["T1", "T2"],
                                              p_values=np.array([[0.02],
                                                                 [0.048]]))
    assert mock_expression_profile == eq_expression_profile
    del eq_expression_profile
    assert mock_expression_profile != ["T", "E", "S", "T"]


def test_from_data_frame(mock_expression_profile):
    """Test import from pandas.DataFrame."""
    data_frame = pd.DataFrame(np.array([[10, 11, 0.02],
                                        [65, 109, 0.048]]),
                              index=["G1", "G2"],
                              columns=["T1", "T2", "T1 T2 p-value"])
    converted_data_frame = ExpressionProfile.from_data_frame(data_frame)
    del data_frame
    assert mock_expression_profile == converted_data_frame


def test_from_csv(mock_expression_profile):
    """Test import from CSV file."""
    csv_data = ExpressionProfile.from_csv("test_expression_profile.csv",
                                          index_col=0)
    assert mock_expression_profile == csv_data


def test_p_value_columns(mock_expression_profile):
    """Test p-value column name lookup."""
    assert mock_expression_profile.p_value_columns == ["T1 T2 p-value"]


def test_p_values(mock_expression_profile):
    """Test p-value lookup."""
    assert (mock_expression_profile.p_values == np.array([[0.02], [0.048]])
            ).all()


def test_to_dict(mock_expression_profile):
    """Test to_dict method of ExpressionProfile."""
    assert mock_expression_profile.to_dict("T1") == {"G1": 10, "G2": 65}


def test_minmax(mock_expression_profile):
    """Test minmax method of ExpressionProfile."""
    assert mock_expression_profile.minmax() == (10, 109)
