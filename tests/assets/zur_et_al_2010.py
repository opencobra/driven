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
Model and expression sets used by Hadas Zur, Eytan Ruppin and Tomer Shlomi in iMAT publication.

Zur, Hadas; Ruppin, Eytan and Shlomi, Tomer (2010) iMAT: an integrative metabolic analysis tool.
Bioinformatics 26 (24): 3140-3142. doi: 10.1093/bioinformatics/btq602
"""

from __future__ import print_function, absolute_import
from warnings import warn

from cobra import Model, Reaction, Metabolite
from driven.data_sets.expression_profile import ExpressionProfile
import numpy as np


def build_model():
    m = Model("Zur et al 2012")

    m1 = Metabolite("M1")
    m2 = Metabolite("M2")
    m3 = Metabolite("M3")
    m4 = Metabolite("M4")
    m5 = Metabolite("M5")
    m6 = Metabolite("M6")
    m7 = Metabolite("M7")
    m8 = Metabolite("M8")
    m9 = Metabolite("M9")
    m10 = Metabolite("M10")

    r1 = Reaction("R1")
    r1.add_metabolites({m3: 1})

    r2 = Reaction("R2")
    r2.add_metabolites({m1: 1})
    r2.gene_reaction_rule = "G1 or G2"

    r3 = Reaction("R3")
    r3.add_metabolites({m2: 1})
    r3.gene_reaction_rule = "G5"

    r4 = Reaction("R4")
    r4.add_metabolites({m1: -1, m10: 1})
    r4.lower_bound = -r4.upper_bound

    r5 = Reaction("R5")
    r5.add_metabolites({m10: -1, m4: 1})
    r5.lower_bound = -r5.upper_bound

    r6 = Reaction("R6")
    r6.add_metabolites({m1: -1, m4: 1})

    r7 = Reaction("R7")
    r7.add_metabolites({m1: -1, m2: -1, m5: 1, m6: 1})
    r7.gene_reaction_rule = "G6"

    r8 = Reaction("R8")
    r8.add_metabolites({m3: -1, m4: -1, m7: 1, m8: 1})
    r8.gene_reaction_rule = "G3"

    r9 = Reaction("R9")
    r9.add_metabolites({m5: -1})

    r10 = Reaction("R10")
    r10.add_metabolites({m6: -1, m9: 1})
    r10.gene_reaction_rule = "G7"

    r11 = Reaction("R11")
    r11.add_metabolites({m7: -1})

    r12 = Reaction("R12")
    r12.add_metabolites({m8: -1})
    r12.gene_reaction_rule = "G4"

    r13 = Reaction("R13")
    r13.add_metabolites({m9: -1})

    m.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8])

    m.objective = r4

    return m


def build_expression_profile():
    warn("The profile is not correct because the picture is black and white (do not use this profile for testing)")
    expression = np.zeros((7, 1))
    expression[0, 0] = .0
    expression[1, 0] = .0
    expression[2, 0] = .0
    expression[3, 0] = .0
    expression[4, 0] = .0
    expression[5, 0] = .0
    expression[6, 0] = .0

    genes = ["G1", "G2", "G3", "G4", "G5", "G6", "G7"]
    conditions = ["Exp"]

    return ExpressionProfile(genes, conditions, expression)

model = build_model()

expression_profile = build_expression_profile()