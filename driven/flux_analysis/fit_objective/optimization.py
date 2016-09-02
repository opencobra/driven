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
from random import Random

import sympy

from inspyred.ec.emo import Pareto, NSGA2


from cameo import config

from driven.flux_analysis.fit_objective.generators import zero_one_binary_generator, zero_one_linear_generator
from driven.flux_analysis.fit_objective.variators import zero_one_binary_variator, zero_one_linear_variator

__all__ = ["FitProfileStrategy"]

One = sympy.singleton.S.One
Add = sympy.Add._from_args
Mul = sympy.Mul._from_args
Real = sympy.RealNumber


def evaluate(model, candidates, coefficients, profiles, evaluators, **kwargs):
    if all(map(lambda v : v == 0, coefficients)):
        return Pareto([1000000] + [100000 for _ in profiles])
    else:
        fitness_list = [f(model, coefficients, candidates, profile, **kwargs) for f, profile in zip(evaluators, profiles)]
        coefficient_sum = sum(coefficients)
        return Pareto([coefficient_sum] + fitness_list)


def zero_one_bounder(candidate, args):
    for i, v in enumerate(candidate):
        candidate[i] = min(v, 1) if v >= 0 else 0
    return candidate


class FitProfileStrategy(object):
    def __init__(self, profiles=[], evaluators=[], binary=True, use_reactions=True, model=None,
                 heuristic_method=NSGA2, **kwargs):
        self._heuristic_method = heuristic_method(Random())
        self.model = model

        if use_reactions:
            self.candidates = [r for r in self.model.exchanges if r.lower_bound >= 0]
        else:
            self.candidates = self.model.metabolites

        self.kwargs = dict(use_reactions=use_reactions)
        self.kwargs.update(kwargs)
        self.profiles = profiles
        self.evaluators = evaluators
        self.binary = binary
        if binary:
            self._heuristic_method.variator = zero_one_binary_variator
            self._generator = zero_one_binary_generator
        else:
            self._heuristic_method.variator = zero_one_linear_variator
            self._generator = zero_one_linear_generator

    def _evaluator(self, candidates, args):
        return [evaluate(self.model, self.candidates, coeffs, self.profiles, self.evaluators, **self.kwargs)
                for coeffs in candidates]

    def run(self, view=config.default_view, maximize=False, **kwargs):
        res = self._heuristic_method.evolve(generator=self._generator,
                                            maximize=maximize,
                                            representation=self.candidates,
                                            bounder=zero_one_bounder,
                                            evaluator=self._evaluator,
                                            binary=self.binary,
                                            view=view,
                                            **kwargs)
        return res
