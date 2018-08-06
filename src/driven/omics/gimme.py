# -*- coding: utf-8 -*-

# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
# Copyright 2018 Synchon Mandal.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Implement GIMME."""

from __future__ import absolute_import

import cobra
from optlang.symbolics import Zero
from six import iteritems


__all__ = ("gimme",)


def gimme(model, expression_profile, cutoff, fraction_of_optimum=0.9,
          condition=0):
    """
    Gene Inactivity Moderated by Metabolism and Expression (GIMME) [1]_.

    Parameters
    ----------
    model: cobra.Model
        A constraint-based model to perform GIMME on.
    expression_profile: ExpressionProfile
        An expression profile to integrate in the model.
    cutoff: float
        The cutoff value to be defined by the user.
    fraction_of_optimum: float
        The fraction of the Required Metabolic Functionalities.
    condition: str or int, optional (default 0)
        The condition from the expression profile.
        If None, the first condition is used.

    Returns
    -------
    context-specific model: cobra.Model
    solution: cobra.Solution

    References
    ----------
    .. [1] Becker, S. a and Palsson, B. O. (2008).
           Context-specific metabolic networks are consistent with experiments.
           PLoS Computational Biology, 4(5), e1000082.
           doi:10.1371/journal.pcbi.1000082
    """
    with model:
        solution = model.slim_optimize()
        prob = model.problem
        rxn_profile = expression_profile.to_reaction_dict(condition, model)

        if model.objective_direction == 'max':
            fix_obj_const = prob.Constraint(model.objective.expression,
                                            lb=fraction_of_optimum * solution,
                                            name="RMF")
        else:
            fix_obj_const = prob.Constraint(model.objective.expression,
                                            ub=fraction_of_optimum * solution,
                                            name="RMF")
        model.add_cons_vars(fix_obj_const)

        coefficients = {rxn_id: cutoff - expression
                        for rxn_id, expression in iteritems(rxn_profile)
                        if cutoff > expression}
        obj_vars = []
        for rxn_id, coefficient in iteritems(coefficients):
            rxn = model.reactions.get_by_id(rxn_id)
            obj_vars.append((rxn.forward_variable, coefficient))
            obj_vars.append((rxn.reverse_variable, coefficient))

        model.objective = prob.Objective(Zero, sloppy=True, direction="min")
        model.objective.set_linear_coefficients({v: c for v, c in obj_vars})
        sol = model.optimize()
        return model, sol
