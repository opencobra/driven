# -*- coding: utf-8 -*-

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

"""Define a general gene or protein expression data store."""

from __future__ import absolute_import

from itertools import combinations
from ast import And, Or, NodeTransformer

import altair as alt
import numpy as np
import pandas as pd
from sympy import Min, Max, Add, Mul, Symbol
from cobra.core.gene import parse_gpr
from driven.utils import get_common_start


class Transform(NodeTransformer):

    def visit_Name(self, node):
        return Symbol(node.id)
    def visit_BoolOp(self, node):
        values = [self.visit(v) for v in node.values]
        if isinstance(node.op, And):
            return Mul(*values)
        if isinstance(node.op, Or):
            return Add(*values)
        return node

def parse_expr(s):
    a, genes = parse_gpr(s.strip())
    a = Transform().visit(a)
    return a.body


class ExpressionProfile(object):
    """
    Representation of an expression profile.

    It can be RNA-Seq, Proteomics, TNSeq or any other profile that
    links genes/proteins to a value (continuous or discrete). It is
    the storage of single or multiple conditions as well as p-values.

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
        """
        Instantiate a new ExpressionProfile.

        Parameters
        ----------
        identifiers: list or array-like
            Genes/Proteins for the dataset.
        conditions: list or array-like
            Conditions for the dataset.
        expression: numpy.ndarray
            Expression data being used.
        p_values: numpy.ndarray, optional (default None)
            p-values obtained from the expresion data.

        """
        array_error = "Not an array-like structure."
        dimension_error = "Expression data and label dimensions don't match."
        assert isinstance(identifiers, (list, np.ndarray)), array_error
        assert isinstance(conditions, (list, np.ndarray)), array_error
        assert isinstance(expression, np.ndarray), array_error
        assert isinstance(p_values, (type(None), np.ndarray)), array_error
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
        """
        Index the ExpressionProfile.

        Parameters
        ----------
        item: [str, str] or [int, int] or [slice, slice]

        """
        attribute_error = ("Non-supported indexing method. "
                           "Use profile['gene', 'condition'] "
                           "or profile[1, 2].")
        if not isinstance(item, tuple):
            raise AttributeError(attribute_error)

        if isinstance(item[0], str) and isinstance(item[1], str):
            i = self.identifier_index[item[0]]
            j = self._condition_index[item[1]]
        elif isinstance(item[0], (slice, int)) and isinstance(item[1],
                                                              (slice, int)):
            i = item[0]
            j = item[1]
        else:
            raise AttributeError(attribute_error)

        return self.expression[i, j]

    def __eq__(self, other):
        """
        Check equality of two ExpressionProfiles.

        Parameters
        ----------
        other: ExpressionProfile or other object
            The object to check equality against.

        """
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
        Read expression data from a pandas.DataFrame.

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
    def from_csv(cls, file_path, replicates=None, **kwargs):
        """
        Read expression data from a comma separated values (csv) file.

        Parameters
        ----------
        file_path: str
            The file path.
        replicates: int, optional (default None)
            Number of replicates. (uses median of replicates if not None).

        Returns
        -------
        ExpressionProfile

        """
        data = pd.read_csv(file_path, **kwargs)
        if replicates:
            columns = data.columns
            data = pd.DataFrame([data[columns[i:i+replicates]].median(axis=1)
                                 for i in range(0, len(columns), replicates)]
                                ).transpose()
            data.columns = [get_common_start(*columns[i:i+replicates].tolist())
                            for i in range(0, len(columns), replicates)]
        return cls.from_data_frame(data)

    @property
    def data_frame(self):
        """
        Build a pandas.DataFrame from the ExpressionProfile.

        Returns
        -------
        pandas.DataFrame

        """
        if self._p_values is None:
            expression = self.expression
            conditions = self.conditions
        else:
            expression = np.concatenate((self.expression, self.p_values),
                                        axis=1)
            conditions = self.conditions + self.p_value_columns

        return pd.DataFrame(expression,
                            index=self.identifiers,
                            columns=conditions)

    @property
    def p_value_columns(self):
        """
        Generate the p-value column names.

        Returns
        -------
        list

        """
        return ["{} {} p-value".format(c[0], c[1])
                for c in combinations(self.conditions, 2)]

    @property
    def p_values(self):
        """
        Return the p-values between the conditions.

        Returns
        -------
        numpy.ndarray

        Raises
        ------
        ValueError

        """
        if not self._p_values.all():
            raise ValueError("No p-values defined.")
        else:
            return self._p_values

    @p_values.setter
    def p_values(self, p_values):
        """
        Set p_values.

        Parameters
        ----------
        p_values: numpy.ndarray

        """
        assert isinstance(p_values, (np.ndarray, type(None)))
        if p_values is not None:
            if p_values.shape[1] != len(self.p_value_columns):
                raise ValueError("Argument p-values do not cover all \
                                  conditions")

        self._p_values = p_values

    @p_values.deleter
    def p_values(self):
        """Delete p_values."""
        self._p_values = None

    def differences(self, p_value=0.005):
        """
        Calculate the differences based on the MADE method.

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
                p_val = self.p_values[idx, i-1]
                if p_val <= p_value:
                    if start < end:
                        diff[iden].append(+1)
                    elif start > end:
                        diff[iden].append(-1)
                    else:
                        diff[iden].append(0)
                else:
                    diff[iden].append(0)
        return diff

    def _map_gene_to_rxn(self, reaction, gene_values, by):
        """
        Map gene data to reactions (using GPR associations).

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
        rule = reaction.gene_reaction_rule
        expression = parse_expr(rule)
        if by == "or2max_and2min":
            expression = expression.replace(Mul, Min).replace(Add, Max)
        elif by == "or2sum_and2min":
            expression = expression.replace(Mul, Min)
        return expression.subs(gene_values).evalf()

    def to_dict(self, condition):
        """
        Build a dict with identifiers and the expression for the condition.

        Parameters
        ----------
        condition: str or int
            The condition label or the column index.

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

    def to_reaction_dict(self, condition, model, **kwargs):
        """
        Build a dict with reactions and the expression for the condition.

        If identifiers in the profile are genes, GPR associations are
        respected. If the gene expression data is missing in the profile,
        the reaction(s) related to that gene, is(are) supplied with the cutoff
        value.

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
        cutoff = 0.0
        map_by = "or2max_and2min"
        gene_exp = self.to_dict(condition)
        rxn_exp = {}
        for rxn in model.reactions:
            if rxn.genes is not None and any([gene.id in self.identifiers
                                              for gene in rxn.genes]):
                gene_values = {gene.id: gene_exp.get(gene.id, cutoff)
                               for gene in rxn.genes}
                rxn_exp[rxn.id] = self._map_gene_to_rxn(rxn, gene_values,
                                                        by=map_by)
        return rxn_exp

    def hist(self, **kwargs):
        """
        Plot an intensity histogram of the expression data.

        Returns
        -------
        altair.Chart

        """
        df = self.data_frame.melt()
        hist = alt.Chart(df, title="Histogram of expression values", width=400,
                         height=300).mark_bar().encode(
                             x=alt.X("value:Q", bin=True,
                                     title="Expression values"),
                             y="count()")
        rule = alt.Chart(df).mark_rule(color='red').encode(
            x="mean(value):Q",
            size=alt.value(3),
            tooltip="mean(value)")
        return hist + rule

    def scatter(self, x, y, **kwargs):
        """
        Generate a scatter plot for comparison between two conditions.

        Parameters
        ----------
        x: str
            The condition to plot on X-axis.
        y: str
            The condition to plot on Y-axis.

        Returns
        -------
        altair.Chart

        Raises
        ------
        AssertionError

        """
        if len(self.conditions) == 1:
            raise AssertionError("Cannot build a scatter plot with only one \
                                  condition.")

        if isinstance(x, str) and isinstance(y, str):
            df = self.data_frame.reset_index()
            scatter = alt.Chart(df, title="Expression values {} vs. \
                                {}".format(x, y), width=400, height=300
                                ).mark_point().encode(
                                    x="{}:{}".format(x, "Q"),
                                    y="{}:{}".format(y, "Q"),
                                    tooltip=[x, y],
                                    color="index:N",
                                    shape=alt.Color("index:N",
                                                    title="Identifiers")
                                ).interactive()
            return scatter
        else:
            raise AssertionError("Column names should be of type str.")

    def heatmap(self, **kwargs):
        """
        Generate a heatmap from the expression profile.

        Returns
        -------
        altair.Chart

        """
        df = self.data_frame.reset_index().melt(id_vars='index',
                                                var_name='x')
        heatmap = alt.Chart(df, title="Expression profile heatmap", width=400,
                            height=300).mark_rect().encode(
                                x=alt.X('x:O', title="Conditions"),
                                y=alt.Y('index:O', title="Identifiers"),
                                color='value:Q')
        return heatmap

    def box(self, **kwargs):
        """
        Generate a box plot for comparison between conditions.

        Returns
        -------
        altair.Chart

        """
        df = self.data_frame.melt()
        lower_box = 'q1(value):Q'
        lower_whisker = 'min(value):Q'
        upper_box = 'q3(value):Q'
        upper_whisker = 'max(value):Q'

        lower_plot = alt.Chart(df, title="Box plot of expression values",
                               width=400, height=300).mark_rule().encode(
                                   y=alt.Y(lower_whisker,
                                           axis=alt.Axis(title="Expression \
                                                         value")),
                                   y2=lower_box,
                                   x=alt.X('variable:O',
                                           axis=alt.Axis(title="Conditions")))

        middle_plot = alt.Chart(df).mark_bar(size=10.0).encode(
            y=lower_box,
            y2=upper_box,
            x='variable:O')

        upper_plot = alt.Chart(df).mark_rule().encode(
            y=upper_whisker,
            y2=upper_box,
            x='variable:O')

        middle_tick = alt.Chart(df).mark_tick(
            color='white',
            size=5.0).encode(
            y='median(value):Q',
            x='variable:O')

        return lower_plot + middle_plot + upper_plot + middle_tick

    def minmax(self, condition=None):
        """
        Return the min and max values for the specified condition.

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
        return np.amin(values), np.amax(values)

    def _repr_html_(self):
        return self.data_frame._repr_html_()
