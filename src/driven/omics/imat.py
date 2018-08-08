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

"""Implement the iMAT tool."""

from __future__ import absolute_import

from optlang.symbolics import add
from six import iteritems


__all__ = ("imat",)


def imat(model, expression_profile, cutoff, epsilon=1, condition=0):
    r"""
    Apply iMAT [1]_ to obtain context-specific metabolic network.

    Integrate Metabolic Analysis Tool (iMAT) does not assume that an
    objective is to be achieved by the cell, thus allowing generation
    of context-specific models of mammalian tissues (does not have a
    trivial objective function) as well. The core concept of iMAT is
    to maximize the number of matches between reaction states (i.e.,
    active or inactive) and corresponding data states (i.e.,
    expressed or not expressed).

    Parameters
    ----------
    model: cobra.Model
        A constraint-based model to perform iMAT on.
    expression_profile: ExpressionProfile
        The expression profile.
    cutoff: 2-tuple of floats (low, high)
        The cut-off value tuple for expression values.
    epsilon: float
        Positive flux threshold.
    condition: str or int, optional (default 0)
        The condition from the expression profile.
        If None, the first condition is used.

    Returns
    -------
    context-specific model: cobra.Model
    solution: cobra.Solution

    Notes
    -----
    The MILP formulation for iMAT is:

    maximize: (\sum_{i\in R_H} (y_i^+ + y_i^-) + \sum_{i\in R_L} (y_i^+))
    s.t.    : Sv = 0
              v_min <= v <= v_max
              v_i + y_i^+ (v_{min,i} - \epsilon) >= v_{min,i}, i \in R_H
              v_i + y_i^- (v_{max,i} + \epsilon) <= v_{max,i}, i \in R_H
              v_{min,i}(1 - y_i^+) <= v_i <= v_{max,i}(1 - y_i^+), i \in R_L
              v \in R^m
              y_i^+ , y_i^- \in [0, 1]

    References
    ----------
    .. [1] Shlomi, Tomer & N Cabili, Moran & Herrgård, Markus & Ø Palsson,
           Bernhard & Ruppin, Eytan. (2008).
           Network-based prediction of human tissue-specific metabolism.
           Nature biotechnology. 26. 1003-10.
           doi:10.1038/nbt.1487.

    """
    low_cutoff, high_cutoff = cutoff
    if low_cutoff > high_cutoff:
        raise ValueError("Low cutoff value greater than high cutoff.")

    y_vars = []
    x_vars = []
    consts = []

    with model:
        prob = model.problem
        rxn_profile = expression_profile.to_reaction_dict(condition, model)

        for rxn_id, expression in iteritems(rxn_profile):
            rxn = model.reactions.get_by_id(rxn_id)
            if expression > high_cutoff:
                y_pos = prob.Variable("y_pos_{}".format(rxn_id), type="binary")
                y_neg = prob.Variable("y_neg_{}".format(rxn_id), type="binary")
                y_vars.append((y_pos + y_neg))

                pos_const = prob.Constraint(
                    rxn.flux_expression + (y_pos * (rxn.lower_bound -
                                                    epsilon)),
                    lb=rxn.lower_bound, name="y_{}_upper".format(rxn_id))

                neg_const = prob.Constraint(
                    rxn.flux_expression + (y_neg * (rxn.upper_bound +
                                                    epsilon)),
                    ub=rxn.upper_bound, name="y_{}_lower".format(rxn_id))

                consts.extend([y_pos, y_neg, pos_const, neg_const])

            elif expression < low_cutoff:
                x_var = prob.Variable("x_{}".format(rxn_id), type="binary")
                x_vars.append(x_var)

                pos_const = prob.Constraint(
                    ((1 - x_var) * rxn.upper_bound) - rxn.flux_expression,
                    lb=0, name="x_{}_upper".format(rxn_id))

                neg_const = model.problem.Constraint(
                    ((1 - x_var) * rxn.lower_bound) - rxn.flux_expression,
                    ub=0, name="x_{}_lower".format(rxn_id))

                consts.extend([x_var, pos_const, neg_const])

        model.add_cons_vars(consts)
        model.objective = prob.Objective(add(y_vars) + add(x_vars),
                                         sloppy=True, direction="max")

        sol = model.optimize()
        return model, sol
