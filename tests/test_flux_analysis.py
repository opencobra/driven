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


class TranscriptomicsTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_gimme(model, objective=None, *args, **kwargs):
        raise NotImplementedError

    def test_imat(model, objective=None, *args, **kwargs):
        raise NotImplementedError

    def test_made(model, objective=None, *args, **kwargs):
        raise NotImplementedError

    def test_eflux(model, objective=None, *args, **kwargs):
        raise NotImplementedError

    def test_relatch(model, objective=None, *args, **kwargs):
        raise NotImplementedError

    def test_gx_fba(model, objective=None, *args, **kwargs):
        raise NotImplementedError


class ThermodynamicsTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_(self):
        pass


class C13TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_fba(self):
        pass