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

from itertools import combinations
from bokeh.charts import Histogram
from numpy import ndarray, vectorize

from driven.data_sets.normalization import or2min_and2max
from pandas import DataFrame


class ExpressionProfile(object):
    """
    Representation of an Expression profile. It can be RNA-Seq, Proteomics, TNSeq or any other
    profile that links genes/proteins to a value (continuous or discrete).

    It the storage of single or multiple conditions as well as p-values.

    Attributes
    ----------

    identifiers: list
        The gene or protein ids
    conditions: list
        The conditions in the expression profile (time points, media conditions, etc...)
    expression: numpy.ndarray
        An 2 dimensional array (nxm) where n is the number of genes and m the number of conditions.
    p_values: numpy.ndarray
        The p-values between conditions.
    """
    @classmethod
    def from_csv(cls, file_path):
        """
        Reads and expression profile from a Comma Separated Values (CSV) file.

        Parameters
        ----------
        file_path: str
            The path to load.

        Returns
        -------
        ExpressionProfile
            An expression profile built from the file.

        """
        data_frame = DataFrame.from_csv(file_path)
        return cls.from_data_frame(data_frame)

    @classmethod
    def from_data_frame(cls, data_frame):
        """
        Reads and expression profile from a pandas.DataFrame.

        Parameters
        ----------
        data_frame: pandas.DataFrame
            A DataFrame containing the genes as index and the columns as conditions.
            For more information about p-values see @ExpressionProfile.p_value_columns

        Returns
        -------
        ExpressionProfile
            An expression profile built from the DataFrame.
        """
        columns = list(data_frame.columns)
        conditions = [c for c in columns if "p-value" not in c]

        p_value_keys = [c for c in columns if "p-value" not in c]
        if len(p_value_keys) > 0:
            p_values = data_frame[p_value_keys].values
        else:
            p_values = None

        expression = data_frame[conditions].values
        identifiers = list(data_frame.index)
        return ExpressionProfile(identifiers, conditions, expression, p_values)

    def __init__(self, identifiers, conditions, expression, p_values=None):
        assert isinstance(identifiers, list)
        assert isinstance(conditions, list)
        assert isinstance(expression, ndarray)
        assert expression.shape == (len(identifiers), len(conditions))

        self.conditions = conditions
        self._condition_index = dict((c, i) for i, c in enumerate(conditions))
        self.identifiers = identifiers
        self._gene_index = dict((g, i) for i, g in enumerate(identifiers))
        self.expression = expression
        self._p_values = p_values

    def __getitem__(self, item):
        if not isinstance(item, tuple):
            raise AttributeError(
                "Non supported slicing method. E.g. profile[1,2] or profile[\"id1\", \"condition_a\"]")

        if isinstance(item[0], str):
            i = self._gene_index[item[0]]
        elif isinstance(item[0], (slice, int)):
            i = item[0]
        else:
            raise AttributeError(
                "Non supported slicing value. E.g. profile[1,2] or profile[\"id1\", \"condition_a\"]")

        if isinstance(item[1], str):
            j = self._condition_index[item[1]]
        elif isinstance(item[1], (slice, int)):
            j = item[1]
        else:
            raise AttributeError(
                "Non supported slicing method. E.g. profile[1,2] or profile[\"id1\", \"condition_a\"]")

        return self.expression[i, j]

    def __cmp__(self, other):
        if not isinstance(other, ExpressionProfile):
            return False
        else:
            return self.identifiers == other.identifiers and \
                   self.conditions == other.conditions and \
                   self._p_values == other._p_values and \
                   self.expression == other.expression

    def _repr_html_(self):
        return self.data_frame._repr_html_()

    @property
    def data_frame(self):
        """
        Builds a pandas.DataFrame from the ExpressionProfile.

        Returns
        -------
        pandas.DataFrame
            A DataFrame
        """
        if self._p_values is None:
            return DataFrame(self.expression,
                             index=self.identifiers,
                             columns=self.conditions)

        else:
            return DataFrame(self.expression+self.p_values, index=self.identifiers,
                             columns=self.conditions+self.p_value_columns)



    @property
    def p_value_columns(self):
        """
        Generates the p-value column names. The p-values are between conditions.

        Returns
        -------
        list
            A list with p-value column headers.
        """
        return ["%s %s p-value" % c for c in combinations(self.conditions, 2)]

    @property
    def p_values(self):
        return self._p_values

    @p_values.setter
    def p_values(self, p_values):
        assert isinstance(p_values, (ndarray, type(None)))
        if p_values is not None:
            if p_values.shape[1] != len(self.p_value_columns):
                raise ValueError("Argument p_values does not cover all conditions (expected %s)" % self.p_value_columns)

        self._p_values = p_values

    @p_values.deleter
    def p_values(self):
        self._p_values = None

    def histogram(self, transform=float, bins=100):
        df = self.data_frame
        transform = vectorize(transform)

        return Histogram(DataFrame(transform(df.values), index=df.index, columns=df.columns), bins=bins, legend=True)

    def to_dict(self, condition):
        """
        Builds a dict with genes as keys and the expression value for the selected condition.

        Parameters
        ----------
        condition: str or int
            The condition or the index.

        Returns
        -------
        dict
        """
        if isinstance(condition, int):
            index = condition
        else:
            index = self._condition_index[condition]
        return dict(zip(self.identifiers, self.expression[:, index]))

    def to_reaction_dict(self, condition, model, cutoff, normalization=or2min_and2max):
        gene_exp = self.to_dict(condition)
        reaction_exp = {}
        for r in model.reactions:
            reaction_genes = r.genes
            if len(reaction_genes) > 0 and any([identifier.id in self.identifiers for identifier in reaction_genes]):
                reaction_exp[r.id] = normalization(r, {g.id: gene_exp.get(g.id, cutoff) for g in reaction_genes})

        return reaction_exp

    def differences(self, p_value=0.005):
        diff = {}
        for index, gene in enumerate(self.identifiers):
            diff[gene] = []
            for i in range(1, len(self.conditions)):
                start, end = self.expression[index, i-1: i+1]
                p = self.p_values[index, i-1]
                if p <= p_value:
                    if start < end:
                        diff[gene].append(+1)
                    elif start > end:
                        diff[gene].append(-1)
                    else:
                        diff[gene].append(0)
                else:
                    diff[gene].append(0)

        return diff