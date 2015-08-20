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


def zero_one_binary_generator(random, args):
    max_objectives = args.get('max_objectives', 5)
    representation = args.get('_representation')
    individual = [0 for _ in representation]
    for i in range(len(representation) - 1):
        individual[i] = 0 if random.rand() < 0.5 else 1 if individual.count(1) <= max_objectives else 0

    return individual


def zero_one_linear_generator(random, args):
    max_objectives = args.get('max_objectives', 5)
    representation = args.get('_representation')
    individual = [0 for _ in representation]
    for i in range(len(representation) - 1):
        individual[i] = 0 if random.rand() < 0.5 else random.rand() if individual.count(1) <= max_objectives else 0

    return individual