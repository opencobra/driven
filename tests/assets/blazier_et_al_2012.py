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

"""
Model and expression sets used by Anna S. Blazier and Jason A. Papin in their review.

Blazier AS and Papin JA (2012) Integration of expression data in genome-scale metabolic network reconstructions.
Front. Physio. 3:299. doi: 10.3389/fphys.2012.00299
"""


from cobra import Model, Reaction
from cobra import Metabolite
from driven.data_sets.expression_profile import ExpressionProfile
import numpy as np


def build_model():
    m = Model("Blazier et al 2012")
    m1_e = Metabolite("M1_e")
    m1 = Metabolite("M1")
    m2 = Metabolite("M2")
    m3 = Metabolite("M3")
    m4_e = Metabolite("M4_e")
    m4 = Metabolite("M4")
    m5 = Metabolite("M5")

    r1 = Reaction("R1")
    r1.add_metabolites({m1_e: -1, m1: 1})

    r2 = Reaction("R2")
    r2.add_metabolites({m1: -1, m2: 1})
    r2.gene_reaction_rule = "Gene2"

    r3 = Reaction("R3")
    r3.add_metabolites({m2: -1, m3: 1})
    r3.gene_reaction_rule = "Gene3"

    r4 = Reaction("R4")
    r4.add_metabolites({m3: -1})

    r5 = Reaction("R5")
    r5.add_metabolites({m4_e: -1, m4: 1})

    r6 = Reaction("R6")
    r6.add_metabolites({m4: -1, m5: 1})
    r6.gene_reaction_rule = "Gene6"

    r7 = Reaction("R7")
    r7.add_metabolites({m5: -1, m2: 1})
    r7.lower_bound = -r7.upper_bound
    r7.gene_reaction_rule = "Gene7"

    r8 = Reaction("R8")
    r8.add_metabolites({m5: -1})

    m.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8])

    EX_M1_e = m.add_boundary(m1_e)
    EX_M1_e.lower_bound = -10

    EX_M4_e = m.add_boundary(m4_e)
    EX_M4_e.lower_bound = -10

    m.objective = r4

    return m


def build_expression_profile():
    expression = np.zeros((4, 3))
    expression[0, 0] = 0.17
    expression[0, 1] = 0.20
    expression[0, 2] = 0.93
    expression[1, 0] = 0.36
    expression[1, 1] = 0.83
    expression[1, 2] = 0.77
    expression[2, 0] = 0.87
    expression[2, 1] = 0.65
    expression[2, 2] = 0.07
    expression[3, 0] = 0.55
    expression[3, 1] = 0.49
    expression[3, 2] = 0.52

    genes = ["Gene2", "Gene3", "Gene6", "Gene7"]
    conditions = ["Exp#1", "Exp#2", "Exp#3"]

    return ExpressionProfile(genes, conditions, expression)

model = build_model()

expression_profile = build_expression_profile()

