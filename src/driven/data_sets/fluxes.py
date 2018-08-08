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

"""Define a general flux constraint data store."""

from __future__ import absolute_import

import numpy as np
import pandas as pd


class FluxConstraints(object):
    """
    Representation of flux constraint data store.

    Attributes
    ----------
    reaction_ids: list
        The reaction IDs to store data for.
    limits: numpy.ndarray
        The limits of the reactions.

    """

    def __init__(self, reaction_ids, limits):
        """
        Instantiate FluxConstraints.

        Parameters
        ----------
        reaction_ids: list
            The list of reaction IDs.
        limits: numpy.ndarray
            The bounds for the reactions.

        """
        array_error = "Not an array-like structure."
        dimension_error = "Dimensions of reactions and limits don't match."
        assert isinstance(reaction_ids, list), array_error
        assert isinstance(limits, np.ndarray), array_error
        assert limits.shape == (len(reaction_ids), 2), dimension_error

        self.reaction_ids = reaction_ids
        self._reaction_id_index = {rxn_id: idx
                                   for idx, rxn_id in enumerate(reaction_ids)}
        self.limits = limits

    def __getitem__(self, item):
        """
        Index FluxConstraints.

        Parameters
        ----------
        item: str or int
            The item to search for.

        """
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self._reaction_id_index[item]
        else:
            raise ValueError("Invalid value {} for retrieving reaction \
                             limits. Use either index or \
                             key".format(item))

        return self.limits[index, :]

    def __iter__(self):
        """Iterate FluxConstraints."""
        return iter(self.reaction_ids)

    def __eq__(self, other):
        """
        Check equality with other object.

        Parameters
        ----------
        other: FluxConstraints or object

        """
        if not isinstance(other, FluxConstraints):
            equality = False
        else:
            equality = self.reaction_ids == other.reaction_ids and \
                (self.limits == other.limits).all()
        return equality

    @classmethod
    def from_data_frame(cls, data_frame, data_type="measurement"):
        """
        Instantiate FluxConstraints from pandas.DataFrame.

        Parameters
        ----------
        data_frame: pandas.DataFrame
            The DataFrame to obtain data from.

        Returns
        -------
        FluxConstraints

        """
        reaction_ids = data_frame.index.tolist()
        if data_type == "measurement":
            limits = data_frame.apply(lambda x: [x["value"] - x["deviation"],
                                                 x["value"] + x["deviation"]],
                                      axis=1).values
        elif data_type == "constraints":
            limits = data_frame.apply(lambda x: [x["lower_limit"],
                                                 x["upper_limit"]],
                                      axis=1).values
        else:
            raise ValueError("Invalid input type {}".format(data_type))

        return FluxConstraints(reaction_ids, limits)

    @classmethod
    def from_csv(cls, file_path, data_type="measurement", **kwargs):
        """
        Instantiate FluxConstraints from CSV file.

        Parameters
        ----------
        file_path: str
            The file path of the .csv file.
        data_type: str, optional (default "measurement")
            The type of data being imported.

        Returns
        -------
        FluxConstraints

        """
        return cls.from_data_frame(pd.read_csv(file_path, **kwargs),
                                   data_type=data_type)

    @property
    def data_frame(self):
        """Return a pandas.DataFrame."""
        return pd.DataFrame(self.limits, index=self.reaction_ids,
                            columns=["lower_limit", "upper_limit"])

    def _repr_html_(self):
        return self.data_frame._repr_html_()
