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
from cameo.core.result import FluxDistributionResult
from cameo.util import TimeMachine

from sympy import Add

add = Add._from_args


def gimme(model, objective=None, *args, **kwargs):
    raise NotImplementedError


def imat(model, r_highly=None, r_lowly=None, epsilon=0.001, *args, **kwargs):
    y_variables = list()
    x_variables = list()
    constraints = list()
    try:
        for rid in r_highly:
            reaction = model.reactions.get_by_id(rid)
            y_pos = model.solver.interface.Variable("y_%s_pos" % rid, type="binary")
            y_neg = model.solver.interface.Variable("y_%s_neg" % rid, type="binary")

            y_variables.append([y_neg, y_pos])

            pos_constraint = model.solver.interface.Constraint(
                reaction.flux_expression + y_pos * (reaction.lower_bound - epsilon),
                lb=reaction.lower_bound,
                name="pos_highly_expressed_constraint_%s" % rid,
                sloppy=True)

            neg_constraint = model.solver.interface.Constraint(
                reaction.flux_expression + y_neg * (reaction.upper_bound + epsilon),
                ub=reaction.upper_bound,
                name="neg_highly_expressed_constraint_%s" % rid,
                sloppy=True)

            constraints.extend([pos_constraint, neg_constraint])

        for rid in r_lowly:
            reaction = model.reactions.get_by_id(rid)
            x = model.solver.interface.Variable("x_%s" % rid, type="binary")
            x_variables.append(x)

            pos_constraint = model.solver.interface.Constraint(
                reaction.flux_expression - (1 - x) * reaction.upper_bound,
                ub=0,
                name="x_%s_upper_limit",
                sloppy=True
            )

            neg_constraint = model.solver.interface.Constraint(
                reaction.flux_expression - (1 - x) * reaction.lower_bound,
                lb=0,
                name="x_%s_lower_limit",
                sloppy=True
            )

            constraints.extend([pos_constraint, neg_constraint])

            for variable in x_variables:
                model.solver._add_variable(variable)

            for variables in y_variables:
                model.solver._add_variable(variables[0])
                model.solver._add_variable(variables[1])

            for constraint in constraints:
                model.solver._add_constraint(constraint)

            objective = model.solver.interface.Objective(
                add([add([add([v[0], v[1]]) for v in y_variables]), add(x_variables)]),
                direction="max"
            )

            with TimeMachine() as tm:
                tm(do=partial(setattr, model, "objective", objective),
                   undo=partial(setattr, model, "objective", model.objective))

                solution = model.solve()
                return FluxDistributionResult(solution)

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