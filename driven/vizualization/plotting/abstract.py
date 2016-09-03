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
from driven.vizualization.utils import golden_ratio


class Plotter(object):
    __default_options__ = {
        'palette': 'Spectral',
        'width': 800,
        'color': "#AFDCEC",
    }

    def __init__(self, **defaults):
        self.__default_options__.update(defaults)

    def _palette(self, palette, *args, **kwargs):
        raise NotImplementedError

    def histogram(self, dataframe, bins=100, width=None, height=None, palette=None, title='Histogram', values=None,
                  groups=None, legend=True):
        raise NotImplementedError

    def scatter(self, dataframe, x=None, y=None, width=None, height=None, color=None, title='Scatter',
                xaxis_label=None, yaxis_label=None, label=None):
        raise NotImplementedError

    def heatmap(self, dataframe, y=None, x=None, values=None, width=None, height=None,
                max_color=None, min_color=None, mid_color=None, title='Heatmap'):
        raise NotImplementedError

    def line(self, dataframe, x=None, y=None, width=None, height=None, groups=None, title="Line"):
        raise NotImplementedError

    def boxplot(self, dataframe, values='value', groups=None, width=None, height=None, palette=None,
                title="BoxPlot", legend=True):
        raise NotImplementedError

    @staticmethod
    def _width_height(width, height):
        if width is None or height is None:
            return golden_ratio(width, height)
        else:
            return width, height

    @classmethod
    def display(cls, plot):
        raise NotImplementedError