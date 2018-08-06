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

"""Define package level fixtures."""

from __future__ import absolute_import

import cobra
import numpy as np
import pytest

from driven.data_sets.expression_profile import ExpressionProfile


@pytest.fixture(scope="session")
def toy_model():
    """
    Build a toy metabolic network model according to the reference [1]_.

    References
    ----------
    .. [1] Blazier AS and Papin JA (2012) Integration of expression data in
           genome-scale metabolic network reconstructions.
           Front. Physio. 3:299. doi: 10.3389/fphys.2012.00299
    """
    model = cobra.Model("Blazier et al 2012")
    r_1 = cobra.Reaction("R1")
    r_2 = cobra.Reaction("R2")
    r_3 = cobra.Reaction("R3")
    r_4 = cobra.Reaction("R4")
    r_5 = cobra.Reaction("R5")
    r_6 = cobra.Reaction("R6")
    r_7 = cobra.Reaction("R7")
    r_8 = cobra.Reaction("R8")

    model.add_reactions([r_1, r_2, r_3, r_4, r_5, r_6, r_7, r_8])

    r_1.reaction = "M1_e -> M1"
    r_2.reaction = "M1 -> M2"
    r_3.reaction = "M2 -> M3"
    r_4.reaction = "M3 ->"
    r_5.reaction = "M4_e -> M4"
    r_6.reaction = "M4 -> M5"
    r_7.reaction = "M5 <-> M2"
    r_8.reaction = "M5 ->"

    r_2.gene_reaction_rule = "Gene2"
    r_3.gene_reaction_rule = "Gene3"
    r_6.gene_reaction_rule = "Gene6"
    r_7.gene_reaction_rule = "Gene7"

    EX_M1_e = model.add_boundary(model.metabolites.M1_e, lb=-10.0)
    EX_M4_e = model.add_boundary(model.metabolites.M4_e, lb=-10.0)

    model.objective = r_4
    return model


@pytest.fixture(scope="session")
def toy_expression_data():
    """
    Build a toy expression data set according to the reference [1]_.

    References
    ----------
    .. [1] Blazier AS and Papin JA (2012) Integration of expression data in
           genome-scale metabolic network reconstructions.
           Front. Physio. 3:299. doi: 10.3389/fphys.2012.00299
    """
    expression = np.array(([0.17, 0.20, 0.93],
                           [0.36, 0.83, 0.77],
                           [0.87, 0.65, 0.07],
                           [0.55, 0.49, 0.52]))

    genes = ["Gene2", "Gene3", "Gene6", "Gene7"]
    conditions = ["Exp#1", "Exp#2", "Exp#3"]

    return ExpressionProfile(genes, conditions, expression)
