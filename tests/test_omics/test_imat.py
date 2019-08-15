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

"""Test the expected functionality of iMAT."""

from __future__ import absolute_import

# import pytest

# from driven.omics import imat


# @pytest.mark.parametrize("cutoff, expected", [
#     ((0.25, 0.75), [("R1", 0.0), ("R2", 0.0)]),
#     ((0.50, 0.75), [("R1", 0.0), ("R2", 0.0)])
# ])
# def test_imat(toy_model, toy_expression_data, cutoff, expected):
#     """Test iMAT."""
#     _, sol = imat(
#         toy_model, toy_expression_data, cutoff, condition="Exp#2")
#     for rxn_id, value in expected:
#         assert sol.fluxes[rxn_id] == pytest.approx(value)
