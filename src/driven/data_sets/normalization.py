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

from sympy import Add, Mul
from sympy.functions.elementary.miscellaneous import Max, Min

from cobra import Reaction


def or2min_and2max(reaction, gene_values):
    assert isinstance(reaction, Reaction)
    assert isinstance(gene_values, dict)
    expression = reaction.gene_expression()
    expression = expression.replace(Mul, Max).replace(Add, Min)
    return expression.subs(gene_values).evalf()
