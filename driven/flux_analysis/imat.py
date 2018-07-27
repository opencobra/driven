# -*- coding: utf-8 -*-
"""Contains function to perform iMAT"""

from __future__ import absolute_import

import six
import cobra
from optlang.symbolics import add

from driven.data_sets import ExpressionProfile


def imat(model, expression_profile, cutoff, epsilon=1, condition=0):
    """
    Integrative Metabolic Analysis Tool[1]

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

    References
    ----------
    .. [1] Shlomi, Tomer & N Cabili, Moran & Herrgård, Markus & Ø Palsson,
           Bernhard & Ruppin, Eytan. (2008).
           Network-based prediction of human tissue-specific metabolism.
           Nature biotechnology. 26. 1003-10.
           doi:10.1038/nbt.1487.
    """

    assert isinstance(model, cobra.Model)
    assert isinstance(expression_profile, ExpressionProfile)
    assert isinstance(cutoff, tuple)
    low_cutoff, high_cutoff = cutoff
    assert isinstance(low_cutoff, float)
    assert isinstance(high_cutoff, float)
    if low_cutoff > high_cutoff:
        raise ValueError("Low cutoff value greater than high cutoff.")

    y_vars = []
    x_vars = []
    consts = []

    with model:
        prob = model.problem
        rxn_profile = expression_profile.to_reaction_dict(condition, model)

        for rxn_id, expression in six.iteritems(rxn_profile):
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
        return (model, sol)
