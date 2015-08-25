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

from cameo.visualization.escher_ext import NotebookBuilder
from IPython.display import display
import os
import six


class EscherViewer(object):
    def __init__(self, data_frame, map_name, color_scales, normalization_functions):
        self.data_frame = data_frame
        self.map_name = map_name
        self.builder = None
        self.color_scales = color_scales
        self.normalization_functions = normalization_functions

    def __call__(self, column):
        reaction_data = dict(self.data_frame[column].apply(self.normalization_functions[column]))
        reaction_data = {r: v for r, v in six.iteritems(reaction_data) if v is not None}
        reaction_scale = self.color_scales[column]
        if self.builder is None:
            self._init_builder(reaction_data, reaction_scale)
        else:
            self.builder.update(reaction_data=reaction_data, reaction_scale=reaction_scale)

    def _init_builder(self, reaction_data, reaction_scale):
        if os.path.isfile(self.map_name):
            self.builder = NotebookBuilder(map_json=self.map_name,
                                           reaction_data=reaction_data,
                                           reaction_scale=reaction_scale,
                                           reaction_no_data_color="lightgray",
                                           reaction_no_data_size=5)
        else:
            self.builder = NotebookBuilder(map_name=self.map_name,
                                           reaction_data=reaction_data,
                                           reaction_scale=reaction_scale,
                                           reaction_no_data_color="lightgray",
                                           reaction_no_data_size=5)
        display(self.builder.display_in_notebook())
