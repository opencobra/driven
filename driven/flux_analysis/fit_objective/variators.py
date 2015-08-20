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


def zero_one_binary_variator(random, candidates, args):
    mutation_rate = args.get("mutation_rate", 0.15)
    new_candidates = [None for _ in candidates]
    for i, c in enumerate(candidates):
        new_candidate = [0 for _ in c]
        for j, v in enumerate(c):
            new_candidate[j] = v if random.rand() < mutation_rate else 1 if random.rand() < 0.5 else 0
        new_candidates[i] = new_candidate

    return new_candidates


def zero_one_linear_variator(random, candidates, args):
    mutation_rate = args.get("mutation_rate", 0.15)
    new_candidates = [None for _ in candidates]
    for i, c in enumerate(candidates):
        new_candidates[i] = [v if random.rand() < mutation_rate else random.rand() for j, v in enumerate(c)]

    return new_candidates