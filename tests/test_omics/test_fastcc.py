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

"""Test the expected functionality of FASTCC."""

from __future__ import absolute_import

from driven.omics import fastcc


def test_fastcc(toy_model_fastcore):
    """Test FASTCC."""
    expected_rxns = ['r1', 'r3', 'r4', 'r5', 'r6']
    fastcc_model = fastcc(toy_model_fastcore)
    fastcc_rxns = [rxn.id for rxn in fastcc_model.reactions]
    assert fastcc_rxns == expected_rxns
