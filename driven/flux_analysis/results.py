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

import six

from math import sqrt
from numpy import log2, nan, zeros

from cameo.core.result import Result
from cameo.flux_analysis.simulation import FluxDistributionResult

from pandas import DataFrame


def _compare_flux_distributions(flux_dist1, flux_dist2, self_key="A", other_key="B"):
    assert isinstance(flux_dist1, FluxDistributionResult)
    assert isinstance(flux_dist2, FluxDistributionResult)
    return FluxDistributionDiff(flux_dist1, flux_dist2, self_key, other_key)

FluxDistributionResult.__sub__ = lambda self, other: _compare_flux_distributions(self, other)


class FluxBasedFluxDistribution(FluxDistributionResult):
    def __init__(self, fluxes, objective_value, c13_flux_distribution, *args, **kwargs):
        super(FluxBasedFluxDistribution, self).__init__(fluxes, objective_value, *args, **kwargs)
        self._c13_fluxes = c13_flux_distribution

    @property
    def data_frame(self):
        index = list(self.fluxes.keys())
        data = zeros((len(self._fluxes.keys()), 3))
        data[:, 0] = [self._fluxes[r] for r in index]
        data[:, 1] = [self._c13_fluxes.get(r, [nan, nan])[0] for r in index]
        data[:, 2] = [self._c13_fluxes.get(r, [nan, nan])[1] for r in index]
        return DataFrame(data, index=index, columns=["fluxes", "c13_lower_limit", "c13_upper_limit"])


class ExpressionBasedResult(FluxDistributionResult):
    def __init__(self, fluxes, objective_value, expression, *args, **kwargs):
        super(ExpressionBasedResult, self).__init__(fluxes, objective_value, *args, **kwargs)
        self.expression = expression


class GimmeResult(ExpressionBasedResult):
    def __init__(self, fluxes, objective_value, fba_fluxes, reaction_expression, cutoff, *args, **kwargs):
        super(GimmeResult, self).__init__(fluxes, objective_value, reaction_expression, *args, **kwargs)
        self._fba_fluxes = fba_fluxes
        self.cutoff = cutoff

    @property
    def data_frame(self):
        index = list(self.fluxes.keys())
        data = zeros((len(self._fluxes.keys()), 4))
        data[:, 0] = [self._fluxes[r] for r in index]
        data[:, 1] = [self._fba_fluxes[r] for r in index]
        data[:, 2] = [self.expression.get(r, nan) for r in index]
        data[:, 3] = [self.reaction_inconsistency_score(r) for r in index]
        return DataFrame(data, index=index, columns=["gimme_fluxes", "fba_fluxes", "expression", "inconsistency_scores"])

    def reaction_inconsistency_score(self, reaction):
        if reaction in self.expression:
            return abs(self._fluxes[reaction]) * max(self.cutoff - self.expression[reaction], 0)
        else:
            return 0

    @property
    def distance(self):
        return sum([abs(self.fluxes[r]-self._fba_fluxes[r]) for r in self.fluxes.keys()])

    @property
    def inconsistency_score(self):
        return sum([self.reaction_inconsistency_score(reaction) for reaction in self.expression])

    def trim_model(self, model, tm=None):
        for r_id, flux in six.iteritems(self.fluxes):
            expression = self.expression.get(r_id, self.cutoff+1)
            if abs(flux) == 0 and expression < self.cutoff and self.reaction_inconsistency_score(r_id) <= 0:
                model.reactions.get_by_id(r_id).knock_out(tm)


class IMATResult(ExpressionBasedResult):
    def __init__(self, fluxes, objective_value, expression, lower_cutoff, higher_cutoff, epsilon, *args, **kwargs):
        super(IMATResult, self).__init__(fluxes, objective_value, expression, *args, **kwargs)
        self.lower_cutoff = lower_cutoff
        self.higher_cutoff = higher_cutoff
        self.epsilon = epsilon

    @property
    def highly_express(self):
        return {r: self._highly_expressed(r) for r in self.fluxes.keys()}

    def _highly_expressed(self, value):
        if value in self.expression:
            return self.expression[value] >= self.higher_cutoff
        else:
            return nan

    @property
    def lowly_express(self):
        return {r: self._lowly_expressed(r) for r in self.fluxes.keys()}

    def _lowly_expressed(self, value):
        if value in self.expression:
            return self.expression[value] < self.lower_cutoff
        else:
            return nan

    @property
    def data_frame(self):
        index = list(self.fluxes.keys())
        columns = ["fluxes", "expression", "highly_express", "lowly_expressed"]
        data = zeros((len(self._fluxes.keys()), 4))
        data[:, 0] = [self.fluxes[r] for r in index]
        data[:, 1] = [self.expression.get(r, nan) for r in index]
        data[:, 2] = [self._highly_expressed(r) for r in index]
        data[:, 3] = [self._lowly_expressed(r) for r in index]
        return DataFrame(data, columns=columns, index=index)


class FluxDistributionDiff(Result):
    def __init__(self, flux_dist_a, flux_dist_b, a_key="A", b_key="B", *args, **kwargs):
        super(FluxDistributionDiff, self).__init__(*args, **kwargs)
        assert isinstance(flux_dist_a, FluxDistributionResult)
        assert isinstance(flux_dist_b, FluxDistributionResult)
        assert all([rid in flux_dist_a.fluxes for rid in flux_dist_b.fluxes]) and \
               all([rid in flux_dist_b.fluxes for rid in flux_dist_a.fluxes])

        self._a_key = a_key
        self._fluxes_a = flux_dist_a

        self._b_key = b_key
        self._fluxes_b = flux_dist_b

    def _manhattan_distance(self, value):
        return abs(self._fluxes_a[value] - self._fluxes_b[value])

    @property
    def manhattan_distance(self):
        return sum(self._manhattan_distance(rid) for rid in self._fluxes_a.keys())

    def _euclidean_distance(self, value):
        return (self._fluxes_a[value] - self._fluxes_b[value])**2

    @property
    def euclidean_distance(self):
        return sqrt(sum(self._euclidean_distance(rid) for rid in self._fluxes_a.keys()))

    def _activity(self, value, threshold=1e-6):
        value_a = abs(self._fluxes_a[value])
        value_b = abs(self._fluxes_b[value])

        if value_a < threshold and value_b < threshold:
            return None

        else:
            value_a = 1 if value_a > threshold else 0
            value_b = 1 if value_b > threshold else 0

        return float(value_a - value_b)

    @property
    def activity_profile(self):
        return {rid: self._activity(rid) for rid in self._fluxes_a.keys()}

    def _fold_change(self, value):
        value_a = self._fluxes_a[value]
        value_b = self._fluxes_b[value]
        try:
            return (value_a - value_b)/value_a
        except ZeroDivisionError:
            return None

    @property
    def data_frame(self):
        index = list(self._fluxes_a.keys())
        columns = ["fluxes_%s" % self._a_key, "fluxes_%s" % self._b_key,
                   "manhattan_distance", "euclidean_distance",
                   "activity_profile", "fold_change"]
        data = zeros((len(self._fluxes_a.keys()), 6))
        data[:, 0] = [self._fluxes_a.fluxes[r] for r in index]
        data[:, 1] = [self._fluxes_b.fluxes[r] for r in index]
        data[:, 2] = [self._manhattan_distance(r) for r in index]
        data[:, 3] = [self._euclidean_distance(r) for r in index]
        data[:, 4] = [self._activity(r) for r in index]
        data[:, 5] = [self._fold_change(r) for r in index]
        return DataFrame(data, index=index, columns=columns)
