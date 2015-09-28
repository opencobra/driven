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


from functools import partial
import numbers
import traceback
from cameo import fba
from cameo.flux_analysis.simulation import FluxDistributionResult
from cameo.core.solver_based_model import SolverBasedModel
from cameo.flux_analysis.analysis import flux_variability_analysis as fva
from cameo.util import TimeMachine
import six

from sympy import Add, RealNumber
from driven.data_sets.expression_profile import ExpressionProfile
from driven.data_sets.normalization import or2min_and2max
from driven.flux_analysis.results import GimmeResult


def gimme(model, expression_profile=None, cutoff=None, objective=None, objective_dist=None, fraction_of_optimum=0.9,
          normalization=or2min_and2max, condition=None, *args, **kwargs):
    """
    Gene Inactivity Moderated by Metabolism and Expression (GIMME)[1]

    Parameters
    ----------
    model: SolverBased Model
        A constraint based model
    expression_profile: ExpressionProfile
        An expression profile
    cutoff: float
        inactivity threshold
    objective: str or other cameo compatible objective
        The Minimal Required Functionalities (MRF)
    objective_dist: FluxDistributionResult
        A predetermined flux distribution for the objective can be provided (optional)
    fraction_of_optimum: float
        The fraction of the MRF
    normalization: function
        expression profile normalization function
    condition: str, int or None
        The condition from the expression profile. If None (default), the first condition on the profile will be used.
    normalization: function
        The normalization function to convert the gene expression profile into reaction expression profile (default: max)

    Returns
    -------
    GimmeResult

    References
    ----------
    .. [1] Becker, S. a and Palsson, B. O. (2008). Context-specific metabolic networks are consistent with experiments.
       PLoS Computational Biology, 4(5), e1000082. doi:10.1371/journal.pcbi.1000082
    """

    assert isinstance(model, SolverBasedModel)
    assert isinstance(expression_profile, ExpressionProfile)
    assert isinstance(fraction_of_optimum, numbers.Number)
    assert isinstance(cutoff, numbers.Number)

    if objective is None:
        objective = model.objective

    if objective_dist is None:
        objective_dist = fba(model, objective)

    assert isinstance(objective_dist, FluxDistributionResult)

    if objective.direction == 'max':
        fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                               lb=fraction_of_optimum * objective_dist.objective_value,
                                                               name="required metabolic functionalities")
    else:
        fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                               ub=fraction_of_optimum * objective_dist.objective_value,
                                                               name="required metabolic functionalities")
    objective_terms = list()

    condition = expression_profile.conditions[0] if condition is None else condition

    reaction_profile = expression_profile.to_reaction_dict(condition, model, cutoff, normalization)
    coefficients = {r: cutoff - exp if cutoff > exp else 0 for r, exp in six.iteritems(reaction_profile)}

    for rid, coefficient in six.iteritems(coefficients):
        reaction = model.reactions.get_by_id(rid)
        objective_terms.append(coefficient * (reaction.forward_variable + reaction.reverse_variable))

    with TimeMachine() as tm:
        gimme_objective = model.solver.interface.Objective(Add(*objective_terms), direction="min")
        tm(do=partial(setattr, model, "objective", gimme_objective),
           undo=partial(setattr, model, "objective", model.objective.expression))
        tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
           undo=partial(model.solver._remove_constraints, [fix_obj_constraint]))
        solution = model.solve()
        return GimmeResult(solution.fluxes, solution.f, objective_dist.fluxes, reaction_profile, cutoff)


def imat(model, expression_profile=None, low_cutoff=0.25, high_cutoff=0.85, epsilon=0.1, condition=None,
         normalization=max, fraction_of_optimum=0.99, *args, **kwargs):
    """
    Integrative Metabolic Analysis Tool

    Parameters
    ----------
    model: SolverBasedModel
        A constraint-based model
    expression_profile: ExpressionProfile
        The expression profile
    low_cutoff: number
        The cut off value for low expression values
    high_cutoff: number
        The cut off value for high expression values
    epsilon: float
    """

    assert isinstance(model, SolverBasedModel)
    assert isinstance(expression_profile, ExpressionProfile)
    assert isinstance(high_cutoff, numbers.Number)
    assert isinstance(low_cutoff, numbers.Number)

    condition = expression_profile.conditions[0] if condition is None else condition

    reaction_profile = expression_profile.to_reaction_dict(condition, model, normalization)

    y_variables = list()
    x_variables = list()
    constraints = list()
    try:
        for rid, expression in six.iteritems(reaction_profile):
            fva_res = fva(model, reactions=[rid], fraction_of_optimum=fraction_of_optimum)
            if expression > high_cutoff:
                reaction = model.reactions.get_by_id(rid)
                y_pos = model.solver.interface.Variable("y_%s+" % rid, type="binary")
                y_neg = model.solver.interface.Variable("y_%s-" % rid, type="binary")

                y_variables.append([y_neg, y_pos])

                pos_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_pos * (fva_res["lower_bound"][rid] - epsilon),
                    lb=fva_res["lower_bound"][rid],
                    name="pos_highly_%s" % rid)

                neg_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_neg * (fva_res["upper_bound"][rid] + epsilon),
                    ub=fva_res["upper_bound"][rid],
                    name="neg_highly_%s" % rid)

                constraints.extend([pos_constraint, neg_constraint])

            if expression < low_cutoff:
                reaction = model.reactions.get_by_id(rid)
                x = model.solver.interface.Variable("x_%s" % rid, type="binary")
                x_variables.append(x)

                pos_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression - (1 - x) * fva_res["upper_bound"][rid],
                    ub=0,
                    name="x_%s_upper" % rid)

                neg_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression - (1 - x) * fva_res["lower_bound"][rid],
                    lb=0,
                    name="x_%s_lower" % rid)

                constraints.extend([pos_constraint, neg_constraint])

        for variable in x_variables:
            model.solver._add_variable(variable)

        for variables in y_variables:
            model.solver._add_variable(variables[0])
            model.solver._add_variable(variables[1])

        for constraint in constraints:
            model.solver._add_constraint(constraint)

        objective = model.solver.interface.Objective(Add(*[(v[0] + v[1]) for v in y_variables]) + Add(*x_variables),
                                                     direction="max")

        with TimeMachine() as tm:
            tm(do=partial(setattr, model, "objective", objective),
               undo=partial(setattr, model, "objective", model.objective))

            solution = model.solve()
            return FluxDistributionResult(solution.fluxes, solution.f)
    except:
        traceback.print_exc()

    finally:
        model.solver._remove_variables(x_variables)
        model.solver._remove_variables([var for pair in y_variables for var in pair])
        model.solver._remove_constraints(constraints)


def made(model, objective=None, *args, **kwargs):
    raise NotImplementedError


def eflux(model, objective=None, *args, **kwargs):
    raise NotImplementedError


def relatch(model, objective=None, *args, **kwargs):
    raise NotImplementedError


def gx_fba(model, objective=None, *args, **kwargs):
    raise NotImplementedError


def prom(model, objective=None, *args, **kwargs):
    raise NotImplementedError