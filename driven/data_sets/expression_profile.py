# -*- coding: utf8 -*-
""" Contains class definition for gene/protein expression data storage."""
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
from __future__ import absolute_import

from itertools import combinations

from numpy import ndarray, log2, amin, amax
from pandas import DataFrame
import seaborn as sns
from sympy.parsing.ast_parser import parse_expr
from sympy import Add, Mul, Max, Min, Symbol

# from driven.stats import freedman_diaconis
from driven.utils import get_common_start


class ExpressionProfile(object):
    """
    Representation of an Expression profile. It can be RNA-Seq, Proteomics,
    TNSeq or any other profile that links genes/proteins to a value
    (continuous or discrete).

    It is the storage of single or multiple conditions as well as p-values.

    Attributes
    ----------
    identifiers: list
        The gene/protein IDs.
    conditions: list
        The conditions in the expression profile (e.g., time points,
        media conditions, etc.)
    expression: numpy.ndarray
        A 2-dimensional array having no. of genes/proteins as rows and
        no. of conditions as the columns.
    p_values: numpy.ndarray
        The p-values between conditions.
    """
    def __init__(self, identifiers, conditions, expression, p_values=None):
        array_error = "Not an array-like structure."
        dimension_error = "Expression data and label dimensions don't match."
        assert isinstance(identifiers, (list, ndarray)), array_error
        assert isinstance(conditions, (list, ndarray)), array_error
        assert isinstance(expression, ndarray), array_error
        dimension = (len(identifiers), len(conditions))
        if len(identifiers) == 1:
            assert expression.shape[0] == len(conditions)
        else:
            assert expression.shape == dimension, dimension_error

        self.identifiers = identifiers
        self.identifier_index = {iden: idx for idx, iden in
                                 enumerate(identifiers)}
        self.conditions = conditions
        self._condition_index = {cond: idx for idx, cond in
                                 enumerate(conditions)}
        self.expression = expression
        self._p_values = p_values

    def __getitem__(self, item):
        attribute_error = ("Non-supported indexing method. "
                           "Use profile['gene', 'condition'] "
                           "or profile[1, 2].")
        if not isinstance(item, tuple):
            raise AttributeError(attribute_error)

        if isinstance(item[0], str) and isinstance(item[1], str):
            i = self.identifier_index[item[0]]
            j = self._condition_index[item[1]]
        elif isinstance(item[0], (slice, int)) and isinstance(item[1], (slice, int)):
            i = item[0]
            j = item[1]
        else:
            raise AttributeError(attribute_error)

        return self.expression[i, j]

    def __eq__(self, other):
        if not isinstance(other, ExpressionProfile):
            return False
        else:
            if self._p_values is None and other.p_values is None:
                return self.identifiers == other.identifiers and \
                       self.conditions == other.conditions and \
                       self._p_values == other._p_values and \
                       (self.expression == other.expression).all()
            else:
                return self.identifiers == other.identifiers and \
                       self.conditions == other.conditions and \
                       (self._p_values == other._p_values).all() and \
                       (self.expression == other.expression).all()

    @classmethod
    def from_data_frame(cls, data_frame):
        """
        Reads expression data from a pandas DataFrame.

        Parameters
        ----------
        data_frame: pandas.DataFrame
            A DataFrame containing the genes as index and
            the columns as conditions.
            For more information about p-values see
            @ExpressionProfile.p_value_columns

        Returns
        -------
        ExpressionProfile
            The expression profile built from the DataFrame.
        """
        columns = data_frame.columns.tolist()
        conditions = [c for c in columns if "p-value" not in c]
        p_value_keys = [c for c in columns if "p-value" in c]
        if p_value_keys:
            p_values = data_frame[p_value_keys].values
        else:
            p_values = None

        expression = data_frame[conditions].values
        identifiers = data_frame.index.tolist()
        return ExpressionProfile(identifiers, conditions, expression, p_values)

    @classmethod
    def from_csv(cls, file_path, sep=",", replicates=None):
        """
        Reads expression data from a comma separated values (csv) file.

        Parameters
        ----------
        file_path: str
            The file path.
        sep: str, optional (default ",")
            The delimiter to use for file parsing.
        replicates: int, optional (default None)
            Number of replicates. (uses median of replicates if not None).

        Returns
        -------
        ExpressionProfile
            The expression profile built from the file.

        """
        data = DataFrame.from_csv(file_path, sep=sep)
        if replicates:
            columns = data.columns
            data = DataFrame([data[columns[i:i+replicates]].median(axis=1)
                              for i in range(0, len(columns), replicates)]).transpose()
            data.columns = [get_common_start(*columns[i:i+replicates].tolist())
                            for i in range(0, len(columns), replicates)]
        return cls.from_data_frame(data)

    @property
    def data_frame(self):
        """
        Builds a pandas DataFrame from the ExpressionProfile.

        Returns
        -------
        pandas.DataFrame
        """
        if self._p_values is None:
            return DataFrame(self.expression,
                             index=self.identifiers,
                             columns=self.conditions)

        else:
            return DataFrame(self.expression + self.p_values,
                             index=self.identifiers,
                             columns=self.conditions + self.p_value_columns)

    @property
    def p_value_columns(self):
        """
        Generates the p-value column names.
        The p-values are between conditions.

        Returns
        -------
        list
        """
        return ["%s %s p-value" % c for c in combinations(self.conditions, 2)]

    @property
    def p_values(self):
        """
        Returns the p-values between the conditions
        defined for the identifiers.

        Returns
        -------
        ndarray
        """
        if not self._p_values:
            raise ValueError("No p-values defined.")
        else:
            return self._p_values

    @p_values.setter
    def p_values(self, p_values):
        assert isinstance(p_values, (ndarray, type(None)))
        if p_values is not None:
            if p_values.shape[1] != len(self.p_value_columns):
                raise ValueError("Argument p-values do not cover all \
                                  conditions")

        self._p_values = p_values

    @p_values.deleter
    def p_values(self):
        self._p_values = None

    def differences(self, p_value=0.005):
        """
        Calculates the differences based on the MADE method.

        Parameters
        ----------
        p_value: float, optional (default 0.005)
            A p-value to set as cutoff for calculation.

        Returns
        -------
        dict
        """
        diff = {}
        for idx, iden in enumerate(self.identifiers):
            diff[iden] = []
            for i in range(1, len(self.conditions)):
                start, end = self.expression[idx, i-1: i+1]
                p = self.p_values[idx, i-1]
                if p <= p_value:
                    if start < end:
                        diff[iden].append(+1)
                    elif start > end:
                        diff[iden].append(-1)
                    else:
                        diff[iden].append(0)
                else:
                    diff[iden].append(0)
        return diff

    def map_gene_to_rxn(self, reaction, gene_values, by):
        """
        Maps gene data to reactions (using GPR associations).

        Parameters
        ----------
        reaction: cobra.Reaction
            The reaction to map gene data to.
        gene_values: dict
            The dict of gene IDs as keys and expression values as values.
        by: str
            The method to use for mapping the gene data.
            "or2max_and2min": if more than one gene contribute as a complex,
                              the min of the gene values is taken and if the
                              genes act as isozymes, the max is taken.
            "or2sum_and2min": if more than one gene contribute as a complex,
                              the min of the gene values is taken and if the
                              genes act as isozymes, the sum is taken.

        Returns
        -------
        float
        """
        local_dict = {gene.id: Symbol(gene.id) for gene in reaction.genes}
        rule = reaction.gene_reaction_rule.replace("and", "*").replace("or", "+")
        expression = parse_expr(rule, local_dict)
        if by == "or2max_and2min":
            expression = expression.replace(Mul, Min).replace(Add, Max)
        elif by == "or2sum_and2min":
            expression = expression.replace(Mul, Min)
        return expression.subs(gene_values).evalf()

    def to_dict(self, condition):
        """
        Builds a dict with identifiers as keys and the expression values for
        the selected condition as values.

        Parameters
        ----------
        condition: str or int
            The condition or the column index.

        Returns
        -------
        dict
        """
        if isinstance(condition, int):
            index = condition
        else:
            index = self._condition_index[condition]
        return {identifier: expression for identifier, expression in
                zip(self.identifiers, self.expression[:, index])}

    def to_reaction_dict(self, condition, model, map_by="or2max_and2min",
                         cutoff=0.0):
        """
        Builds a dict with reaction as keys and the expression values for the
        selected condition as values. If identifiers in the profile are genes,
        GPR associations are respected. If the gene expression data is missing
        in the profile, the reaction(s) related to that gene, is(are) supplied
        with the cutoff value.

        Parameters
        ----------
        condition: str or int
            The condition or the column index.
        model: cobra.Model
            The constraint-based model to obtain the GPR associations from.
        map_by: str, optional (default "or2max_and2min")
            The method to use for mapping gene values to reactions.
        cutoff: float, optional (default 0.0)
            The value to set for reaction(s) whose gene lacks a value in
            the profile.

        Returns
        -------
        dict
        """
        gene_exp = self.to_dict(condition)
        rxn_exp = {}
        for rxn in model.reactions:
            if rxn.genes is not None and any([gene.id in self.identifiers
                                              for gene in rxn.genes]):
                gene_values = {gene.id: gene_exp.get(gene.id, cutoff)
                               for gene in rxn.genes}
                rxn_exp[rxn.id] = self.map_gene_to_rxn(rxn, gene_values,
                                                       by=map_by)
        return rxn_exp

# TODO: Refine interface for a better graph; fix the bin to auto
    def hist(self, bins=10, figsize=(6, 6)):
        """
        Plots an intensity histogram of the expression data.

        Parameters
        ----------
        bins: int, optional (default "auto")
            The number of bins to use.
        figsize: tuple (width, height) in inches
            The figure size of the plot.

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
        """
        hist = self.data_frame.plot.hist(bins=bins, figsize=figsize,
                                         legend=True, alpha=0.5,
                                         title="Histogram of expression \
                                                values")
        return hist

    def scatter(self, x, y, figsize=(6, 6), color="#AFDCEC"):
        """
        Generates a scatter plot for comparison between two conditions.

        Parameters
        ----------
        x: str
            The condition to plot on X-axis.
        y: str
            The condition to plot on Y-axis.
        figsize: tuple (width, height) in inches
            The figure size of the plot.
        color: str
            Either name of color like 'red' or hex code like '#AFDCEC'

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot

        Raises
        ------
        AssertionError
        """
        if len(self.conditions) == 1:
            raise AssertionError("Cannot build a scatter plot with only one \
                                  condition.")

        if isinstance(x, str) and isinstance(y, str):
            scatter = self.data_frame.plot.scatter(x, y, c=color,
                                                   figsize=figsize,
                                                   title="Expression values \
                                                          {0} vs. \
                                                          {1}".format(x, y))
            return scatter
        else:
            raise AssertionError("Column names should be of type str.")

    def heatmap(self):
        """
        Generates a heatmap from the expression profile.

        Parameters
        ----------
        figsize: tuple (width, height) in inches
            The figure size of the plot.

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
        """
        heatmap = sns.heatmap(self.data_frame, square=True, cmap="RdYlGn_r")
        heatmap.set_title("Expression profile heatmap")
        return heatmap

    def box(self, by=None, figsize=(6, 6)):
        """
        Generates a box plot for comparison between conditions.

        Parameters
        ----------
        by: str or None (default None)
            The condition to group by.
        figsize: tuple (width, height) in inches
            The figure size of the plot.

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
        """
        box = self.data_frame.plot.box(by=by, figsize=figsize,
                                       title="Box plot of expression values")
        return box

    def minmax(self, condition=None):
        """
        Returns the minimum and maximum values for the
        specified condition.

        Parameters
        ----------
        condition: str or int or None, optional (default None)
            The condition to obtain the min and max values for.

        Returns
        -------
        tuple of (min, max)
        """
        if condition is None:
            values = self[:, :]
        elif isinstance(condition, int):
            values = self[:, condition]
        elif isinstance(condition, str):
            values = self[:, self._condition_index[condition]]
        return amin(values), amax(values)

    def normalize(self, using="log2"):
        """
        Normalizes the expression data using the chosen method.
        By default uses log2 for the purpose.

        Parameters
        ---------
        using: str, optional (default "log2")
            The method to use for normalizing the expression data.
        """
        if using == "log2":
            self.expression = log2(self.expression)

    def _repr_html_(self):
        return self.data_frame._repr_html_()

#     def bin_width(self, condition=None, min_val=None, max_val=None):
#         if condition is None:
#             values = self[:, :]
#         else:
#             values = self[:, condition]

#         if min_val:
#             values = values[values >= min_val]
#         if max_val:
#             values = values[values <= max_val]

#         values = values[:, condition]
#         return freedman_diaconis(values)
