# -*- coding: utf-8 -*-
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
from cameo.util import TimeMachine


def tmfa(model, objective=None, delta_gs=None, K=None, *args, **kwargs):
    """
    Thermodynamics-based Metabolic Flux Analysis. [1]

    Arguments
    ---------

    objective: a compatible objective
        The objective function of the simulation
    delta_gs: dict
        The Gibbs energy of the reactions
    K: number
        The feasible Gibbs energy

    References
    ----------
    .. [1] Henry, C. S., Broadbelt, L. J. and Hatzimanikatis, V. (2007). Thermodynamics-based metabolic flux analysis.
    Biophysical Journal, 92(5), 1792–1805. doi:10.1529/biophysj.106.093138
    """
    variables = []
    constraints = []
    original_objective = model.objective.expression
    try:
        for reaction in model.reactions:
            if reaction.id in delta_gs:
                z = model.solver.interface.Variable("z_%s" % reaction.id, type='binary')
                variables.append(z)

                k_constraint = model.solver.interface.Constraint(delta_gs[reaction.id] - K + K*z,
                                                                 ub=0,
                                                                 name="second_law_constraint_%s" % reaction.id,
                                                                 sloppy=True)

                flux_constraint = model.solver.interface.Constraint(z * reaction.upper_bound - reaction.flux_expression,
                                                                    lb=0,
                                                                    name="thermo_flux_constraint_%s" % reaction.id,
                                                                    sloppy=True)

                constraints.append([k_constraint, flux_constraint])
        with TimeMachine() as tm:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective.expression))
            tm(do=partial(model.solver.add, variables),
               undo=partial(model.solver.remove, variables))
            tm(do=partial(model.solver.add, constraints),
               undo=partial(model.solver.remove, constraints))

    finally:
        model.solver._remove_constraints(constraints)
        model.solver._remove_variables(variables)
        model.objective = original_objective