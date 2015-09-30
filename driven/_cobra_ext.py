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


from cobra import Metabolite, Reaction
from sympy import Symbol
from sympy.parsing.ast_parser import parse_expr


def __gene_to_expression__(reaction):
    local_dict = {g.id: Symbol(g.id) for g in reaction.genes}
    rule = reaction._gene_reaction_rule.replace("and", "+").replace("or", "*")
    return parse_expr(rule, local_dict)


def __str__metabolite__(self):
    if self.formula is not None and len(self.formula.elements) > 0:
        return "%s (%s)" % (self.name, self.formula)
    else:
        return self.name


def _repr_html_metabolite_(self):
    return """
    <table>
        <tbody>
            <tr><td><b>Id</b></td><td>%s</td></tr>
            <tr><td><b>Name</b></td><td>%s</td></tr>
            <tr><td><b>Formula</b></td><td>%s/td></tr>
        </tbody>
    </table>
    """ % (self.id, self.name, self.formula)


Metabolite.__str__ = __str__metabolite__
Metabolite.__repr__ = __str__metabolite__
Metabolite._repr_html_ = _repr_html_metabolite_
Reaction.gene_expression = __gene_to_expression__

