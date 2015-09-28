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

from numpy import ndarray
from pandas import DataFrame


class FluxConstraints(object):

    @classmethod
    def from_csv(cls, file_path, type="measurement"):
        return cls.from_data_frame(DataFrame.from_csv(file_path), type=type)

    @classmethod
    def from_data_frame(cls, data_frame, type="measurement"):
        reaction_ids = list(data_frame.index.values)
        if type == "measurement":
            limits = data_frame.apply(lambda x: [x["value"]-x["deviation"], x["value"]+x["deviation"]], axis=1).values
        elif type == "constraints":
            limits = data_frame.apply(lambda x: [x["lower_limit"], x["upper_limit"]], axis=1).values
        else:
            raise ValueError("Invalid input type %s" % type)

        return FluxConstraints(reaction_ids, limits)

    def __init__(self, reaction_ids, limits):
        assert isinstance(reaction_ids, list), "Invalid class for reactions_ids %s" % type(reaction_ids)
        assert isinstance(limits, ndarray), "Invalid class for limits %s" % type(limits)
        assert limits.shape == (len(reaction_ids), 2)
        self.reaction_ids = reaction_ids
        self._reaction_id_index = {rid: i for i, rid in enumerate(reaction_ids)}
        self.limits = limits

    def __getitem__(self, item):
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self._reaction_id_index[item]
        else:
            raise ValueError("Invalid value %s for retrieving reaction limits. Use either index or key" % item)

        return self.limits[index, :]

    def __iter__(self):
        return iter(self.reaction_ids)

    def __eq__(self, other):
        if not isinstance(other, FluxConstraints):
            return False
        else:
            return self.reaction_ids == other.reaction_ids and (self.limits == other.limits).all()

    def _repr_html_(self):
        return self.data_frame._repr_html_()

    @property
    def data_frame(self):
        return DataFrame(self.limits, index=self.reaction_ids, columns=["lower_limit", "upper_limit"])