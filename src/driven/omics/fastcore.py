# -*- coding: utf-8 -*-

# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Implement the FASTCORE method."""

from __future__ import absolute_import


def fastcore(model, core_reactions, flux_threshold=1.0):
    """
    Apply FASTCORE [1]_ to obtain context-specific metabolic network.

    #TODO: description
    Parameters
    ----------
    model: cobra.Model
        The consistent constraint-based model.
    core_reactions: list of cobra.Reaction.id
        The core reactions that are essential in the context.
    flux_threshold: float, optional (default 1e-9)
        The flux threshold to consider.

    Returns
    -------
    context_specific_model: cobra.Model
        The context-specific constraint-based model.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """

    with model:
        j = [rxn for rxn in core_reactions
             if model.reactions.get_by_id(rxn).reversibility is False]
        p = [rxn for rxn in model.reactions if rxn.id not in core_reactions]
        flipped = False
        a = find_sparse_mode(model, j, p)
        j = [rxn for rxn in core_reactions if rxn not in a]
        while j is not None:
            p = [rxn for rxn in j if rxn not in a]
            a = a.extend(find_sparse_mode(j, p, model))
            if not [rxn for rxn in j if rxn in a]:
                j = [rxn for rxn in j if rxn not in a]
                flipped = False
            else:
                if flipped:
                    flipped = False
                else:
                    flipped = True
                    j_prime = j
                    for i in [rxn for rxn in j if model.reactions.get_by_id(rxn).reversibility is False]:
                        ...
    return a


def find_sparse_mode(model, core_reactions, penalty_reactions,
                     zero_cutoff=1e-9, flux_threshold=1.0):
    """
    Obtain dense reaction set as flux mode.

    #TODO: update docstring
    Parameters
    ----------
    model: cobra.Model
        The constraint-based model.
    core_reactions: list of cobra.Reaction.id
        The reactions to be included in the final network.
    penalty_reactions: list of cobra.Reaction.id
        The reactions not considered core but may be added.

    Returns
    -------
    list of cobra.Reaction.id

    """
    prob = model.problem

    if core_reactions:
        with model:
            obj_vars = []
            vars_and_cons = []
            for reaction in core_reactions:
                rxn = model.reactions.get_by_id(reaction)
                var = prob.Variable("auxillary_{}".format(rxn.id),
                                    lb=0.0, ub=flux_threshold)
                const = prob.Constraint(rxn.forward_variable + rxn.reverse_variable
                                        - var, name="constraint_{}".format(rxn.id),
                                        lb=0.0)
                vars_and_cons.extend([var, const])
                obj_vars.append(var)

            model.add_cons_vars(vars_and_cons)
            model.objective = prob.Objective(Zero, sloppy=True, direction="max")
            model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
            sol = model.optimize()

            rxns_to_consider = [rxn_id for rxn_id, flux in sol.fluxes.iteritems() if
                                abs(flux) > zero_cutoff]

            if rxns_to_consider:
                obj_vars = []
                vars_and_cons = []
                for reaction in penalty_reactions:
                    rxn = model.reactions.get_by_id(reaction)
                    var = prob.Variable("auxillary_penalty_{}".format(rxn.id),
                                        lb...)
