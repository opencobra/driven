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


from cameo import config
from cameo.strain_design.heuristic.optimization import HeuristicOptimization
from inspyred.ec.emo import Pareto, NSGA2
import sympy
from driven.flux_analysis.fit_objective.generators import zero_one_binary_generator, zero_one_linear_generator
from driven.flux_analysis.fit_objective.variators import zero_one_binary_variator, zero_one_linear_variator


One = sympy.singleton.S.One
Add = sympy.Add._from_args
Mul = sympy.Mul._from_args
Real = sympy.RealNumber


def evaluate(model, reactions, coefficients, profiles, evaluators):
    if all(map(lambda v : v == 0, coefficients)):
        return Pareto([1000000] + [100000 for _ in profiles])
    else:
        fitness_list = [f(model, coefficients, reactions, profile) for f, profile in zip(evaluators, profiles)]
        coefficient_sum = sum(coefficients)
        return Pareto([coefficient_sum] + fitness_list)


def zero_one_bounder(candidate, args):
    for i, v in enumerate(candidate):
        candidate[i] = min(v, 1) if v >= 0 else 0
    return candidate


class FitProfileStrategy(HeuristicOptimization):
    def __init__(self, profiles=[], evaluators=[], binary=True, heuristic_method=NSGA2, *args, **kwargs):
        HeuristicOptimization.__init__(self, heuristic_method=heuristic_method, *args, **kwargs)
        self.reactions = [r for r in self.model.exchanges if r.lower_bound >= 0]

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
        return [evaluate(self.model, self.reactions, coeffs, self.profiles, self.evaluators) for coeffs in candidates]

    def run(self, view=config.default_view, maximize=False, **kwargs):
        return super(FitProfileStrategy, self).run(view=view,
                                                   maximize=maximize,
                                                   _representation=self.reactions,
                                                   bounder=zero_one_bounder,
                                                   binary=self.binary,
                                                   **kwargs)
