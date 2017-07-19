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

from __future__ import absolute_import, print_function

import sympy

from cobra import Reaction
from cobra.flux_analysis import find_essential_genes

__all__ = ["essential_genes_profile_evaluator", "essential_reactions_profile_evaluator"]

One = sympy.singleton.S.One
Add = sympy.Add._from_args
Mul = sympy.Mul._from_args
Real = sympy.RealNumber


def _essential_profile_assertion_score(actual, expected):
    score = 0
    for key in actual:
        exp = expected.get(key, False)
        if exp == actual[key]:
            score += 1
    return float(score)


def essential_genes_profile_evaluator(model, coefficients, candidates, essential, use_reactions=True):
    total = float(len(model.genes))
    with model:
        _set_objective(model, coefficients, candidates, use_reactions)
        predicted_essential = find_essential_genes(model)
        actual = {g.id: g in predicted_essential for g in model.genes}
        return total - _essential_profile_assertion_score(actual, essential)


def essential_reactions_profile_evaluator(model, coefficients, candidates, essential, use_reactions=True):
    total = float(len(model.reactions) - len(model.exchanges))
    with model:
        _set_objective(model, coefficients, candidates, use_reactions)
        predicted_essential = model.essential_reactions()
        actual = {r.id: r in predicted_essential for r in model.reactions if r not in model.exchanges}
        return total - _essential_profile_assertion_score(actual, essential)


def _set_objective(model, coefficients, candidates, use_reactions):
    if use_reactions:
        obj = Add([Mul([Real(coeff), react.flux_expression]) for react, coeff in zip(candidates, coefficients)])
    else:
        obj = Reaction("make_metabolites")
        obj.add_metabolites({met: coeff for met, coeff in zip(candidates, coefficients) if coeff != 0})
        model.add_reactions([obj])
    model.objective = obj
