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

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero


def _find_sparse_mode(model, rxns, penalty_rxns, flux_threshold, zero_cutoff):
    """Perform the LPs required for FASTCORE.

    Parameters
    ----------
    model: cobra.core.Model
        The cobra model to perform FASTCC on.
    rxns: list of cobra.Reaction
        The reactions to use for LP.
    penalty_rxns: list of cobra.Reaction
        The reactions to keep as penalty reactions.
    flux_threshold: float
        The upper threshold an auxiliary variable can have.
    zero_cutoff: float
        The cutoff below which flux is considered zero.

    Returns
    -------
    result: list
        The list of reactions to consider as consistent.

    """

    if rxns:
        # perform LP-7 according to the publication
        obj_vars = []
        vars_and_cons = []
        prob = model.problem

        for rxn in rxns:
            var = prob.Variable("auxiliary_{}".format(rxn.id),
                                lb=0.0, ub=flux_threshold)
            const = prob.Constraint(rxn.forward_variable +
                                    rxn.reverse_variable -
                                    var, name="constraint_{}".format(rxn.id),
                                    lb=0.0)
            vars_and_cons.extend([var, const])
            obj_vars.append(var)

        model.add_cons_vars(vars_and_cons)
        model.objective = prob.Objective(Zero, sloppy=True)
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})

        model.optimize(objective_sense="max")

        interim_result = [rxn for rxn in rxns if rxn.flux > zero_cutoff]

        if interim_result:
            # perform LP-10 according to the publication
            obj_vars = []
            vars_and_cons = []
            prob = model.problem

            for rxn in penalty_rxns:
                var = prob.Variable("auxiliary_{}".format(rxn.id),
                                    lb=0.0, ub=flux_threshold)
                const_1 = prob.Constraint(
                    rxn.forward_variable + rxn.reverse_variable - var,
                    name="constraint_1_{}".format(rxn.id), ub=0.0
                )
                const_2 = prob.Constraint(
                    rxn.forward_variable + rxn.reverse_variable + var,
                    name="constraint_2_{}".format(rxn.id), lb=0.0
                )
                vars_and_cons.extend([var, const_1, const_2])
                obj_vars.append(var)

            for rxn in interim_result:
                const = prob.Constraint(
                    rxn.forward_variable - rxn.reverse_variable - zero_cutoff,
                    name="constraint_interim_var_{}".format(rxn.id), lb=0.0
                )
                vars_and_cons.extend([const])

            model.add_cons_vars(vars_and_cons)
            model.objective = prob.Objective(Zero, sloppy=True)
            model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})

            model.optimize(objective_sense="min")

            result = [rxn for rxn in model.reactions
                      if abs(rxn.flux) > zero_cutoff]
        else:
            result = []
    else:
        result = []

    return result


def _flip_coefficients(model, rxns):
    """Flip the coefficients for optimizing in reverse direction."""

    # flip reactions
    for rxn in rxns:
        const = model.constraints.get("constraint_{}".format(rxn.id))
        var = model.variables.get("auxiliary_{}".format(rxn.id))
        coefs = const.get_linear_coefficients(const.variables)
        const.set_linear_coefficients({k: -v for k, v in coefs.items()
                                       if k is not var})

    # flip objective
    objective = model.objective
    objective_coefs = objective.get_linear_coefficients(objective.variables)
    objective.set_linear_coefficients({k: -v for k, v in
                                       objective_coefs.items()})


def fastcore(model, active_reactions, flux_threshold=1.0, zero_cutoff=None):
    """
    Generate a context-specific model using FASTCORE [1]_.

    FASTCORE is a rapid and efficient algorithm to generate reconstructions
    of consistent models. This implementation requires that you provide a
    consistent model generated using FASTCC [1]_ or similar algorithms.
    FASTCORE is a generic algorithm as it relies only on a preselected set of
    reactions and a simple mathematical objective (flux consistency). Thus, it
    allows easy integration of multi-omics data and can be used in different
    contexts. For more details on FASTCORE, please check [1]_.

    Parameters
    ----------
    model: cobra.Model
        The consistent constraint-based model to operate on.
    active_reactions: list of str
        The reaction (IDs) that must be active in reconstruction.
    flux_threshold: float, optional (default 1.0)
        The flux threshold to consider.
    zero_cutoff: float, optional
        The cutoff to consider for zero flux (default model.tolerance).

    Returns
    -------
    context-specific model: cobra.Model
        The reconstructed context-specific COBRA model.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """

    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    active_reactions = [model.reactions.get_by_id(rxn)
                        for rxn in active_reactions]

    irreversible_rxns = [rxn for rxn in model.reactions
                         if not rxn.reversibility]
    rxns_to_check = list(set(active_reactions).intersection(irreversible_rxns))
    penalty_rxns = list(set(model.reactions).difference(active_reactions))

    with model:
        rxns_to_keep = _find_sparse_mode(
            model,
            rxns_to_check,
            penalty_rxns,
            flux_threshold,
            zero_cutoff
        )

    rxns_to_check = list(set(active_reactions).difference(rxns_to_keep))

    while rxns_to_check:
        penalty_rxns = list(set(penalty_rxns).difference(rxns_to_keep))

        with model:
            new_rxns = _find_sparse_mode(
                model,
                rxns_to_check,
                penalty_rxns,
                flux_threshold,
                zero_cutoff
            )
            rxns_to_keep.extend(new_rxns)

            # this condition will be valid for all but the last iteration
            if list(set(rxns_to_check).intersection(rxns_to_keep)):
                rxns_to_check = list(
                    set(rxns_to_check).difference(rxns_to_keep)
                )

            else:
                rxns_to_flip = list(
                    set(rxns_to_check).difference(irreversible_rxns)
                )
                _flip_coefficients(model, rxns_to_flip)
                sol = model.optimize(min)
                to_add_rxns = sol.fluxes.index[sol.fluxes.abs() > zero_cutoff]\
                                        .tolist()
                rxns_to_keep.extend([model.reactions.get_by_id(rxn)
                                     for rxn in to_add_rxns])
                # since this is the last iteration, it needs to break or else
                # it will run forever since rxns_to_check won't be empty
                break

    consistent_rxns = set(rxns_to_keep)
    # need the ids since Reaction objects are created fresh with model.copy()
    rxns_to_remove = [rxn.id for rxn in
                      set(model.reactions).difference(consistent_rxns)]

    consistent_model = model.copy()
    consistent_model.remove_reactions(rxns_to_remove, remove_orphans=True)

    return consistent_model
