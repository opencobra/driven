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
from cobra import Reaction
from sympy import Symbol, Add, Mul
from sympy.parsing.ast_parser import parse_expr
from sympy.functions.elementary.miscellaneous import Max, Min


def or2min_and2max(reaction, gene_expression):
    assert isinstance(reaction, Reaction)
    assert isinstance(gene_expression, dict)

    local_dict = {g.id: Symbol(g.id) for g in reaction.genes}
    gene_expression = {gid: gene_expression.get(gid, 0) for gid in local_dict.keys()}
    rule = reaction.gene_reaction_rule.replace("and", "+").replace("or", "*")
    expression = parse_expr(rule, local_dict)
    expression = expression.replace(Mul, Max).replace(Add, Min)
    return expression.evalf(subs=gene_expression)

