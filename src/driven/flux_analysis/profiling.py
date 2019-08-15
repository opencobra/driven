# _-*- coding: utf-8 -*-

"""Define classes for knockout profilers."""

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

from __future__ import absolute_import

# from cobra.manipulation.delete import find_gene_knockout_reactions


# __all__ = ["ReactionKnockoutProfiler", "GeneKnockoutProfiler"]


# class KnockoutProfiler(object):
#     def __init__(self, model=None, simulation_method=None):
#         self._model = model
#         self._simulation_method = simulation_method

#     def profile(self, elements, *args, **kwargs):
#         return {e.id: self.simulate_knockout(e, *args, **kwargs)
#                 for e in elements}

#     def simulate_knockout(self, to_knockout, *args, **kwargs):
#         raise NotImplementedError


# class ReactionKnockoutProfiler(KnockoutProfiler):
#     def __init__(self, *args, **kwargs):
#         super(ReactionKnockoutProfiler, self).__init__(*args, **kwargs)

#     def simulate_knockout(self, to_knockout, *args, **kwargs):
#         cache = {
#             "variables": {},
#             "constraints": {},
#             "first_run": True
#         }

#         with self._model:
#             to_knockout.knock_out()
#             return self._simulation_method(
#                 self._model,
#                 volatile=False,
#                 cache=cache,
#                 *args,
#                 **kwargs)[cache['original_objective']]


# class GeneKnockoutProfiler(KnockoutProfiler):
#     def __init__(self, *args, **kwargs):
#         super(GeneKnockoutProfiler, self).__init__(*args, **kwargs)

#     def simulate_knockout(self, to_knockout, *args, **kwargs):
#         reactions = find_gene_knockout_reactions(self._model, [to_knockout])
#         cache = {
#             "variables": {},
#             "constraints": {},
#             "first_run": True
#         }

#         with self._model:
#             for reaction in reactions:
#                 reaction.knock_out()

#             return self._simulation_method(
#                 self._model,
#                 volatile=False,
#                 cache=cache,
#                 *args,
#                 **kwargs)[cache['original_objective']]
