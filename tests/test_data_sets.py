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
import unittest

import numpy as np

from driven.data_sets.expression_profile import ExpressionProfile
from driven.data_sets.fluxes import FluxConstraints


class ExpressionProfileTestCase(unittest.TestCase):
    def test_difference(self):
        genes = ["G1"]
        conditions = ["T1", "T2", "T3", "T4"]

        expression = np.zeros((1, 4))
        expression[0] = [10, 11, 65, 109]
        pvalues = np.zeros((1, 3))
        pvalues[0] = [0.02, 0.048, 0.0012]

        profile = ExpressionProfile(genes, conditions, expression, pvalues)

        self.assertEqual(profile.differences(), {"G1": [0, 0, 1]})

    def test_export_import(self):
        genes = ["G1"]

        conditions = ["T1", "T2", "T3", "T4"]

        expression = np.zeros((1, 4))
        expression[0] = [10, 11, 65, 109]

        profile = ExpressionProfile(genes, conditions, expression)
        data_frame = profile.data_frame

        new_profile = ExpressionProfile.from_data_frame(data_frame)

        self.assertEqual(profile, new_profile)


class FluxConstraintsTestCase(unittest.TestCase):
    def test_export_import(self):
        reaction_ids = ["R1", "R2", "R3"]
        limits = np.zeros((3, 2))
        limits[0] = [0, 10]
        limits[1] = [0.5, 0.7]
        limits[2] = [5.1, 5.2]
        flux_constraints = FluxConstraints(reaction_ids, limits)

        new_flux_constraints = FluxConstraints.from_data_frame(flux_constraints.data_frame, type="constraints")

        print(new_flux_constraints.data_frame)

        self.assertEqual(flux_constraints, new_flux_constraints)

