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
from functools import partial
import traceback
from cameo import fba
from cameo.core.result import FluxDistributionResult
from cameo.util import TimeMachine
import six

from sympy import Add, RealNumber
from driven.flux_analysis.results import GimmeResult


def gimme(model, expression_profile=None, cutoff=None, objective=None, fraction_of_optimum=0.9, *args, **kwargs):
    """

    :param model: SolverBased Model
    :param expression_profile: ExpressionProfile
    :param cutoff: float
    :param reference_objective: str or other cameo compatible objective
    :param fraction_of_optimum
    :param args:
    :param kwargs:
    :return:
    """

    if objective is None:
        objective = model.objective

    fba_result = fba(model, objective)

    if objective.direction == 'max':
        fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                               lb=fraction_of_optimum * fba_result.objective_value,
                                                               name="required metabolic functionalities")
    else:
        fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                               ub=fraction_of_optimum * fba_result.objective_value,
                                                               name="required metabolic functionalities")
    objective_terms = list()

    coefficients = {rid: 0 if expression <= cutoff else cutoff - expression
                    for rid, expression in six.iteritems(expression_profile)}

    for rid, coefficient in six.iteritems(coefficients):
        reaction = model.reactions.get_by_id(rid)
        objective_terms.append(coefficient * (reaction.forward_variable + reaction.reverse_variable))

    with TimeMachine() as tm:
        tm(do=partial(setattr, model, "objective", Add(*objective_terms)),
           undo=partial(setattr, model, "objective", model.objective.expression))
        tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
           undo=partial(model.solver._remove_constraints, [fix_obj_constraint]))

        return GimmeResult(model.solve(), fba_result.fluxes, expression_profile, cutoff)


def imat(model, expression_profile=None, low_cutoff=0.25, high_cutoff=0.85, epsilon=0.001, *args, **kwargs):
    y_variables = list()
    x_variables = list()
    constraints = list()
    try:
        for rid, expression in six.iteritems(expression_profile):
            if expression > high_cutoff:
                print "Highly %s" % rid
                reaction = model.reactions.get_by_id(rid)
                y_pos = model.solver.interface.Variable("y_%s_pos" % rid, type="binary")
                y_neg = model.solver.interface.Variable("y_%s_neg" % rid, type="binary")

                y_variables.append([y_neg, y_pos])

                pos_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_pos * (reaction.lower_bound - epsilon),
                    lb=reaction.lower_bound,
                    name="pos_highly_expressed_constraint_%s" % rid)

                neg_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_neg * (reaction.upper_bound + epsilon),
                    ub=reaction.upper_bound,
                    name="neg_highly_expressed_constraint_%s" % rid)

                constraints.extend([pos_constraint, neg_constraint])

            if expression < low_cutoff:
                reaction = model.reactions.get_by_id(rid)
                x = model.solver.interface.Variable("x_%s" % rid, type="binary")
                x_variables.append(x)

                pos_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression - (1 - x) * reaction.upper_bound,
                    ub=0,
                    name="x_%s_upper_limit" % rid)

                neg_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression - (1 - x) * RealNumber(reaction.lower_bound),
                    lb=0,
                    name="x_%s_lower_limit" % rid)

                constraints.extend([pos_constraint, neg_constraint])

        for variable in x_variables:
            model.solver._add_variable(variable)

        for variables in y_variables:
            model.solver._add_variable(variables[0])
            model.solver._add_variable(variables[1])

        for constraint in constraints:
            model.solver._add_constraint(constraint)

        objective = model.solver.interface.Objective(
            Add(Add(*[v[0] + v[1] for v in y_variables]), Add(*x_variables)),
            direction="max"
        )

        with TimeMachine() as tm:
            tm(do=partial(setattr, model, "objective", objective),
               undo=partial(setattr, model, "objective", model.objective))
            return FluxDistributionResult(model.solve())
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