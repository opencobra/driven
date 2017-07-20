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
from warnings import warn

from plotly.graph_objs import Figure, Layout, Scatter

from driven.vizualization.plotting import Plotter


class PlotlyPlotter(Plotter):

    def __init__(self, **defaults):
        warn("Plotly requires configuration before start (https://plot.ly/python/getting-started/)")
        super(PlotlyPlotter, self).__init__(**defaults)

    def scatter(self, dataframe, x=None, y=None, width=None, height=None, color=None, title='Scatter', xaxis_label=None,
                yaxis_label=None, label=None):

        color = self.__default_options__.get('color', None) if color is None else color
        width = self.__default_options__.get('width', None) if width is None else width

        scatter = Scatter(x=dataframe[x],
                          y=dataframe[y],
                          mode='markers',
                          marker=dict(color=color))

        if label:
            scatter['text'] = dataframe[label]

        width, height = self._width_height(width, height)

        layout = Layout(title=title,
                        width=width,
                        height=height)

        if xaxis_label:
            layout['xaxis'] = dict(title=xaxis_label)

        if yaxis_label:
            layout['yaxis'] = dict(title=yaxis_label)

        return Figure(data=[scatter], layout=layout)

    def histogram(self, dataframe, bins=100, width=None, height=None, palette=None, title='Histogram', values=None,
                  groups=None, legend=True):
        pass

    @classmethod
    def display(cls, plot):
        pass

    def heatmap(self, dataframe, y=None, x=None, values=None, width=None, height=None, max_color=None, min_color=None,
                mid_color=None, title='Heatmap'):
        pass

    def _palette(self, palette, *args, **kwargs):
        pass
