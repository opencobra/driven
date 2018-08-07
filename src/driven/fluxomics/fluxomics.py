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

"""Constrain models using metabolic flux measurements."""

from __future__ import absolute_import

from driven.flux_analysis.results import FluxBasedFluxDistribution


def fba(model, objective=None, distribution=None, relax=0.01, *args, **kwargs):
    """
    Run a Flux Balance Analysis using a set of measured fluxes.

    Parameters
    ----------
    model : cobra.Model
        A constraint-based model
    distribution : FluxConstraints
        A set of experimental or computational determined flux constraints
    objective : objective
        Optional objective for the model
    relax : float
        A relax value to make the computation feasible

    """
    with model:
        for reaction_id in distribution:
            reaction = model.reactions.get_by_id(reaction_id)
            bounds = distribution[reaction_id]
            reaction.bounds = (bounds[0] - bounds[0] * relax,
                               bounds[1] + bounds[1] * relax)
        if objective is not None:
            model.objective = objective

        solution = model.optimize()
    # TODO: bug here, missing arg.
    return FluxBasedFluxDistribution(solution.fluxes, distribution)
