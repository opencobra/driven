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

"""Test the expected functionality of GIMME."""

from __future__ import absolute_import

# from pytest import approx

# from driven.omics import gimme


# def test_gimme(toy_model, toy_expression_data):
#     """Test GIMME."""
#     _, sol_025 = gimme(toy_model, toy_expression_data, cutoff=0.25)
#     _, sol_050 = gimme(toy_model, toy_expression_data, cutoff=0.5)
#     assert sol_025.objective_value == approx(0.0)
#     assert sol_050.objective_value > 0.0
