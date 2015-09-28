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
from collections import OrderedDict
from math import sqrt

from IPython.core.display import display

import warnings
from cameo.core.result import Result
from driven.generic.normalization import log_plus_one
from driven.vizualization.escher_viewer import EscherViewer

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from IPython.html.widgets import Dropdown
    except ImportError:
        from ipywidgets import Dropdown


from cameo.flux_analysis.simulation import FluxDistributionResult
from pandas import DataFrame
import numpy as np
import six


class FluxBasedFluxDistribution(FluxDistributionResult):
    def __init__(self, fluxes, objective_value, c13_flux_distribution, *args, **kwargs):
        super(FluxBasedFluxDistribution, self).__init__(fluxes, objective_value, *args, **kwargs)
        self._c13_fluxes = c13_flux_distribution

    @property
    def data_frame(self):
        index = list(self.fluxes.keys())
        data = np.zeros((len(self._fluxes.keys()), 3))
        data[:, 0] = [self._fluxes[r] for r in index]
        data[:, 1] = [self._c13_fluxes.get(r, [np.nan, np.nan])[0] for r in index]
        data[:, 2] = [self._c13_fluxes.get(r, [np.nan, np.nan])[1] for r in index]
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
        data = np.zeros((len(self._fluxes.keys()), 4))
        data[:, 0] = [self._fluxes[r] for r in index]
        data[:, 1] = [self._fba_fluxes[r] for r in index]
        data[:, 2] = [self.expression.get(r, float("nan")) for r in index]
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
            if abs(flux) == 0 and self.expression.get(r_id, self.cutoff+1) < self.cutoff:
                model.reactions.get_by_id(r_id).knock_out(tm)

    def display_on_map(self, map_name=None):
        color_scales = {
            "gimme_fluxes":         [dict(type='min', color="yellow", size=20),
                                     dict(type='value', value=0, color="green", size=7),
                                     dict(type='max', color='blue', size=20)],
            "expression":           [dict(type='value', value=log_plus_one(self.cutoff), color="green", size=10),
                                     dict(type='max', color='blue', size=20),
                                     dict(type='min', color='yellow', size=5)],
            "inconsistency_scores": [dict(type='value', value=0, color="yellow", size=10),
                                     dict(type='max', color='green', size=20)]
        }

        normalization_functions = {
            "gimme_fluxes": log_plus_one,
            "expression": log_plus_one,
            "inconsistency_scores": log_plus_one
        }

        viewer = EscherViewer(self.data_frame, map_name, color_scales, normalization_functions)
        drop_down = Dropdown()
        drop_down.options = {
            "Fluxes": "gimme_fluxes",
            "Expression": "expression",
            "Inconsistency Score": "inconsistency_scores"
        }
        drop_down.on_trait_change(lambda x: viewer(drop_down.get_state("value")["value"]))
        display(drop_down)
        viewer("gimme_fluxes")

    def compare(self, other_result, self_key="A", other_key="B"):
        assert isinstance(other_result, FluxDistributionResult)
        return FluxDistributionComparison(self, other_result, self_key, other_key)


class FluxDistributionComparison(Result):
    def __init__(self, flux_dist_a, flux_dist_b, a_key="A", b_key="B", *args, **kwargs):
        super(FluxDistributionComparison, self).__init__(*args, **kwargs)
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
        return {rid: self._manhattan_distance(rid) for rid in self._fluxes_a.keys()}

    def _euclidean_distance(self, value):
        return sqrt((self._fluxes_a[value] - self._fluxes_b[value])**2)

    @property
    def euclidean_distance(self):
        return {rid: self._euclidean_distance(rid) for rid in self._fluxes_a.keys()}

    def _activity(self, value, threshold=1e-6):
        value_a = abs(self._fluxes_a[value])
        value_b = abs(self._fluxes_b[value])

        if value_a < threshold and value_b < threshold:
            return None

        else:
            value_a = 1 if value_a > threshold else 0
            value_b = 1 if value_b > threshold else 0

        return value_a - value_b

    @property
    def activity_profile(self):
        return {rid: self._activity(rid) for rid in self._fluxes_a.keys()}

    @property
    def data_frame(self):
        index = list(self._fluxes_a.keys())
        columns = ["fluxes_%s" % self._a_key, "fluxes_%s" % self._b_key,
                   "manhattan_distance", "euclidean_distance", "activity_profile"]
        data = np.zeros((len(self._fluxes_a.keys()), 5))
        data[:, 0] = [self._fluxes_a.fluxes[r] for r in index]
        data[:, 1] = [self._fluxes_b.fluxes[r] for r in index]
        data[:, 2] = [self._manhattan_distance(r) for r in index]
        data[:, 3] = [self._euclidean_distance(r) for r in index]
        data[:, 4] = [self._activity(r) for r in index]
        return DataFrame(data, index=index, columns=columns)

    def plot(self, grid=None, width=None, height=None, title=None):
        pass

    def display_on_map(self, map_name):
        color_scales = {
            "fluxes_%s" % self._a_key: [dict(type='min', color="yellow", size=20),
                                        dict(type='value', value=0, color="blue", size=7),
                                        dict(type='max', color='green', size=20)],
            "fluxes_%s" % self._b_key: [dict(type='min', color="yellow", size=20),
                                        dict(type='value', value=0, color="blue", size=7),
                                        dict(type='max', color='green', size=20)],
            "manhattan_distance":      [dict(type='min', color="yellow", size=20),
                                        dict(type='value', value=0, color="blue", size=7),
                                        dict(type='max', color='green', size=20)],
            "euclidean_distance":      [dict(type='min', color="yellow", size=20),
                                        dict(type='value', value=0, color="blue", size=7),
                                        dict(type='max', color='green', size=20)],
            "activity_profile":        [dict(type='value', value=-1, color="yellow", size=10),
                                        dict(type='value', value=0, color="green", size=10),
                                        dict(type='value', value=1, color='blue', size=10)],
        }

        normalization_functions = {
            "fluxes_%s" % self._a_key: log_plus_one,
            "fluxes_%s" % self._b_key: log_plus_one,
            "manhattan_distance": float,
            "euclidean_distance": float,
            "activity_profile": int
        }

        viewer = EscherViewer(self.data_frame, map_name, color_scales, normalization_functions)
        drop_down = Dropdown()
        drop_down.options = OrderedDict({
            "Flux Distribution %s" % self._a_key: "fluxes_%s" % self._a_key,
            "Flux Distribution %s" % self._b_key: "fluxes_%s" % self._b_key,
            "Manhattan Distance": "manhattan_distance",
            "Euclidean Distance": "euclidean_distance",
            "Activity Profile": "activity_profile"
        })
        drop_down.value = "fluxes_%s" % self._a_key
        drop_down.on_trait_change(lambda x: viewer(drop_down.get_state("value")["value"]))
        display(drop_down)
        viewer("fluxes_%s" % self._a_key)
