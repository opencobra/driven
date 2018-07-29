# -*- coding: utf-8 -*-
"""Contains function to perform GIM3E"""

from __future__ import absolute_import

import six
import cobra
from optlang.symbolics import Zero
from cobra.util import fix_objective_as_constraint

from driven.data_sets import ExpressionProfile


def gim3e(model, expression_profile, metabolite_profile, media,
          optimum_fraction=0.9, penalty_fraction=1.01, condition=0):
    """
    Gene Inactivation Moderated by Metabolism, Metabolomics and Expression[1]

    Parameters
    ----------
    model: cobra.Model
        A constraint-based model to perform GIMME on.
    expression_profile: ExpressionProfile
        An expression profile to integrate in the model.
    metabolite_profile: dict
        A dictionary of medium condition as keys and list of
        detected metabolites as values.
    media: str
        The key(media) for the metabolite_profile being used.
    optimum_fraction: float, optional (default 0.9)
        The fraction of the Required Metabolic Functionalities.
    penalty_fraction: float, optional (default 1.01)
        The fraction of the penalty bound.
    condition: str or int, optional (default 0)
        The condition from the expression profile.
        If None, the first condition is used.

    Returns
    -------
    context-specfic model: cobra.Model
    solution: cobra.Solution

    References:
    ----------
    .. [1] Brian J. Schmidt, Ali Ebrahim, Thomas O. Metz, Joshua N. Adkins,
           Bernhard Ø. Palsson, Daniel R. Hyduke;
           GIM3E: condition-specific models of cellular metabolism developed
           from metabolomics and expression data, Bioinformatics, Volume 29,
           Issue 22, 15 November 2013, Pages 2900–2908,
           https://doi.org/10.1093/bioinformatics/btt493
    """
    assert isinstance(model, cobra.Model)
    assert isinstance(expression_profile, ExpressionProfile)
    assert isinstance(metabolite_profile, dict)
    assert isinstance(media, str)
    assert isinstance(optimum_fraction, float)
    assert isinstance(penalty_fraction, float)
    assert isinstance(condition, (str, int))

    exp_max, _ = expression_profile.minmax()
    tolerance = model.solver.configuration.tolerances.feasibility

    objective = model.objective
    objective_dir = model.objective_direction

    rxn_profile = expression_profile.to_reaction_dict(condition, model)
    media_metabolites = metabolite_profile[media]

    prob = model.problem

    # Fix objective
    fix_objective_as_constraint(model, fraction=optimum_fraction)

    # Add turnover metabolites
    turnover_mets = []
    for met in model.metabolites:
        turnover_met = cobra.Metabolite("TM_{}".format(met.id))
        turnover_mets.append(turnover_met)
    model.add_metabolites(turnover_mets)

    # Add turnover metabolites to reactions
    for rxn in model.reactions:
        mets_to_add = {}
        for key, value in six.iteritems(rxn.metabolites):
            turnover_met = model.metabolites.get_by_id("TM_{}".format(key.id))
            mets_to_add[turnover_met] = abs(value)
        rxn.add_metabolites(mets_to_add)

    # Add sink reactions
    for met in model.metabolites.query("^TM_"):
        sink_rxn = cobra.Reaction("SK_{}".format(met.id))
        # TAKEN FROM THE ORIGINAL IMPLEMENTATION:
        # Since both creation and consumption of
        # the real metabolite drive the turnover,
        # we require 2 units of metabolite per
        # 1 unit of flux sink so this matches
        # the turnover through the real metabolite.
        # METABOLOMICS DATA:
        sink_rxn.add_metabolites({met: -2})
        if met.id[3:] in media_metabolites:
            sink_rxn.lower_bound = penalty_fraction * tolerance
        else:
            sink_rxn.lower_bound = 0.0
        model.add_reaction(sink_rxn)
    # TRANSCRIPTOMICS DATA:
    # Add penalty coefficients
    coefficients = {rxn_id: exp_max - expression for rxn_id, expression in
                    six.iteritems(rxn_profile)}
    # Determine penalty bounds
    obj_vars = []
    for rxn_id, coefficient in six.iteritems(coefficients):
        rxn = model.reactions.get_by_id(rxn_id)
        obj_vars.append((rxn.forward_variable, coefficient))
        obj_vars.append((rxn.reverse_variable, coefficient))

    model.objective = prob.Objective(Zero, sloppy=True, direction="min")
    model.objective.set_linear_coefficients({v: c for v, c in obj_vars})
    # Penalty bound constraint
    penalty_obj = model.slim_optimize()
    penalty_bound_const = prob.Constraint(model.objective.expression,
                                          name="penalty_objective_bound",
                                          ub=penalty_fraction * penalty_obj)
    model.add_cons_vars(penalty_bound_const)

    model.objective = objective
    model.objective_direction = objective_dir
    sol = model.optimize()

    return (model, sol)
