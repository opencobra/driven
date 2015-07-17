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

from numpy import ndarray
import six


class ExpressionProfile(object):
    def __init__(self, genes, conditions, expression):
        assert isinstance(genes, list)
        assert isinstance(conditions, list)
        assert isinstance(expression, ndarray)

        self.conditions = conditions
        self._condition_index = dict((c, i) for i, c in enumerate(conditions))
        self.genes = genes
        self._gene_index = dict((g, i) for i, g in enumerate(genes))
        self.expression = expression

    def split(self, thresholds=[]):
        pass

    def to_dict(self, condition):
        return dict(zip(self.genes, self.expression[:, self._condition_index[condition]]))

    def to_reaction_dict(self, condition, model, normalization=max):
        reaction_expression = {}
        for g in self.genes:
            gene = model.genes.get_by_id(g)
            expression = self.expression[self._gene_index[g], self._condition_index[condition]]
            for reaction in gene.reactions:
                if reaction.id not in reaction_expression:
                    reaction_expression[reaction.id] = []
                reaction_expression[reaction.id].append(expression)

        return {rid: normalization(exp) for rid, exp in six.iteritems(reaction_expression)}

    def __getitem__(self, item):
        if not isinstance(item, tuple):
            raise AttributeError(
                "Non supported slicing method. E.g. profile[1,2] or profile[\"gene1\", \"condition_a\"]")

        if isinstance(item[0], str):
            i = self._gene_index[item[0]]
        elif isinstance(item[0], (slice, int)):
            i = item[0]
        else:
            raise AttributeError(
                "Non supported slicing value. E.g. profile[1,2] or profile[\"gene1\", \"condition_a\"]")

        if isinstance(item[1], str):
            j = self._condition_index[item[1]]
        elif isinstance(item[1], (slice, int)):
            j = item[1]
        else:
            raise AttributeError(
                "Non supported slicing method. E.g. profile[1,2] or profile[\"gene1\", \"condition_a\"]")

        return self.expression[i, j]

    @property
    def differences(self):
        diff = {}
        for index, gene in enumerate(self.genes):
            diff[gene] = []
            for i in range(1, len(self.conditions)):
                start, end = self.expression[index, i-1: i+1]
                if start < end:
                    diff[gene].append(+1)
                elif start > end:
                    diff[gene].append(-1)
                else:
                    diff[gene].append(0)

        return diff