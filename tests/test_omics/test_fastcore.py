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

"""Test the expected functionality of FASTCORE."""

from __future__ import absolute_import

import pytest
from cobra import Model, Reaction

from driven.omics import fastcore


@pytest.fixture(scope="module")
def toy_model_fastcore():
    """
    Build a toy model for testing FASTCORE from reference [1]_.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """
    model = Model("FASTCORE")
    r_1 = Reaction("v1")
    r_2 = Reaction("v2")
    r_4 = Reaction("v4")
    r_5 = Reaction("v5")
    r_6 = Reaction("v6")
    r_7 = Reaction("v7")
    r_8 = Reaction("v8")

    model.add_reactions([r_1, r_2, r_4, r_5, r_6, r_7, r_8])

    r_1.reaction = "-> A"
    r_2.reaction = "A -> B"
    r_4.reaction = "B -> D"
    r_5.reaction = "D ->"
    r_6.reaction = "E -> D"
    r_7.reaction = "F -> E"
    r_8.reaction = "A -> F"

    model.objective = r_4
    return model


def test_fastcore(toy_model_fastcore):
    """Test FASTCORE."""
    expected_rxns = ['r1', 'r2', 'r4', 'r5']
    expected_mets = ['A', 'B', 'D']

    fastcore_model = fastcore(toy_model_fastcore)
    fastcore_rxns = [rxn.id for rxn in fastcore_model.reactions]
    fastcore_mets = [met.id for met in fastcore_model.metabolites]

    assert fastcore_rxns == expected_rxns
    assert fastcore_mets == expected_mets
