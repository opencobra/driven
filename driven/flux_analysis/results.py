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
from cameo.core.result import FluxDistributionResult
from pandas import DataFrame
import numpy as np


class GimmeResult(FluxDistributionResult):
    def __init__(self, solution, fba_fluxes, expression_profile, cutoff, *args, **kwargs):
        super(GimmeResult, self).__init__(solution, *args, **kwargs)
        self._fba_fluxes = fba_fluxes
        self.expression_profile = expression_profile
        self.cutoff = cutoff

    @property
    def data_frame(self):
        index = list(self.fluxes.keys())
        data = np.zeros((8, 3))
        data[:, 0] = list(self._fluxes.values())
        data[:, 1] = list(self._fba_fluxes.values())
        data[:, 2] = [self.reaction_inconsistency_score(r) for r in index]
        return DataFrame(data, index=index, columns=["gimme_fluxes", "fba_fluxes", "inconsistency_scores"])

    def reaction_inconsistency_score(self, reaction):
        if reaction in self.expression_profile:
            return abs(self._fluxes[reaction]) * (self.cutoff - self.expression_profile[reaction])
        else:
            return 0

    @property
    def inconsistency_score(self):
        return sum([self.reaction_inconsistency_score(reaction) for reaction in self.expression_profile])