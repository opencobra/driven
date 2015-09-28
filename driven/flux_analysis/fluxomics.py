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
from cameo.solver_based_model import SolverBasedModel
from cameo.util import TimeMachine
from driven.data_sets.fluxes import FluxConstraints
from driven.flux_analysis.results import FluxBasedFluxDistribution


def fba(model, objective=None, distribution=None, relax=0.01, *args, **kwargs):
    """
    Runs a Flux Balance Analysis using a set of measured fluxes.

    Arguments
    ---------

    model: SolverBasedModel
        A constraint-based model
    distribution: FluxConstraints
        A set of experimental or computational determined flux constraints
    relax: float
        A relax value to make the computation feasible
    """
    if not isinstance(model, SolverBasedModel):
        raise ValueError("Argument model must be instance of SolverBasedModel, not %s" % model.__class__)
    if not isinstance(distribution, FluxConstraints):
        raise ValueError("Argument distribution must be instance of FluxConstraints, not %s" % distribution.__class__)

    with TimeMachine() as tm:
        for reaction_id in distribution:
            reaction = model.reactions.get_by_id(reaction_id)
            bounds = distribution[reaction_id]
            tm(do=partial(setattr, reaction, "upper_bound", bounds[1] + bounds[1] * relax),
               undo=partial(setattr, reaction, "upper_bound", reaction.upper_bound))
            tm(do=partial(setattr, reaction, "lower_bound", bounds[0] - bounds[0] * relax),
               undo=partial(setattr, reaction, "lower_bound", reaction.lower_bound))

        if objective is not None:
            tm(do=partial(setattr, model, "objective", objective),
               undo=partial(setattr, model, "objective", model.objective))

        solution = model.solve()

    return FluxBasedFluxDistribution(solution.fluxes, distribution)