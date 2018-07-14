# -*- coding: utf-8 -*-
"""Contains test class and functions for the methods."""
# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
import sys
import numpy as np
import cobra

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..')))

from driven import ExpressionProfile
from driven.flux_analysis import gimme, imat


def toy_model():
    """
    Builds a toy model based on the article:
    Blazier AS and Papin JA (2012) Integration of expression data in
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


def toy_expression_data():
    """
    Builds the toy expression data from the article:
    Blazier AS and Papin JA (2012) Integration of expression data in
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


class TestTranscriptomics(object):
    """
    Base class for testing methods for integrating transcriptomics
    data in genome-scale models.
    """

    def test_gimme(self):
        """
        Tests GIMME.
        """
        model = toy_model()
        expression = toy_expression_data()

        _, sol_025 = gimme(model, expression, cutoff=0.25)
        _, sol_050 = gimme(model, expression, cutoff=0.5)
        assert np.isclose(sol_025.objective_value, 0.0)
        assert sol_050.objective_value > 0.0

    def test_imat(self):
        """
        Tests iMAT.
        """
        model = toy_model()
        expression = toy_expression_data()

        _, sol_025_075 = imat(model, expression,
                              (0.25, 0.75), condition="Exp#2")
        assert np.isclose(sol_025_075.fluxes["R1"], 0.0)
        assert np.isclose(sol_025_075.fluxes["R2"], 0.0)

        _, sol_050_075 = imat(model, expression,
                              (0.50, 0.75), condition="Exp#2")

        assert np.isclose(sol_050_075.fluxes["R1"], 0)
        assert np.isclose(sol_050_075.fluxes["R2"], 0)

    # def test_made(self):
    #     raise NotImplementedError

    # def test_eflux(self):
    #     raise NotImplementedError

    # def test_relatch(self):
    #     raise NotImplementedError

    # def test_gx_fba(self):
    #     raise NotImplementedError

    # def test_prom(self):
    #     raise NotImplementedError
