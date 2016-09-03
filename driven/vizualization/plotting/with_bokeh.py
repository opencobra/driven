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

import collections
import six

from bokeh.charts import Histogram, Scatter, HeatMap, Line, BoxPlot
from bokeh.models import HoverTool, GlyphRenderer
from bokeh.io import show
from bokeh.palettes import brewer

from driven.vizualization.plotting.abstract import Plotter


TOOLS = 'pan,box_zoom,wheel_zoom,resize,reset,save'


class BokehPlotter(Plotter):
    def __init__(self, **defaults):
        super(BokehPlotter, self).__init__(**defaults)

    def _palette(self, palette, number, **kwargs):
        if isinstance(palette, six.string_types):
            n = 3 if number < 3 else number
            return brewer[palette][n]
        elif isinstance(palette, collections.Iterable):
            return palette
        else:
            raise ValueError("Invalid palette %s" % palette)

    def histogram(self, dataframe, bins=None, width=None, height=None, palette=None, title='Histogram', values=None,
                  groups=None, legend=True):
        palette = self.__default_options__.get('palette', None) if palette is None else palette
        width = self.__default_options__.get('width', None) if width is None else width

        if values:
            unique_values = dataframe[groups].unique()
            palette = self._palette(palette, len(unique_values))
        else:
            palette = None

        width, height = self._width_height(width, height)

        histogram = Histogram(dataframe, values=values, color=groups, bins=bins, legend=True, width=width,
                              height=height, palette=palette, title=title)

        return histogram

    def scatter(self, dataframe, x=None, y=None, width=None, height=None, color=None, title=None,
                xaxis_label=None, yaxis_label=None, label=None):
        color = self.__default_options__.get('palette', None) if color is None else color
        width = self.__default_options__.get('width', None) if width is None else width

        width, height = self._width_height(width, height)

        scatter = Scatter(dataframe, x=x, y=y, width=width, height=height, color=color, title=title,
                          tools=TOOLS + ',hover' if label else '')

        if label:
            hover = scatter.select_one(dict(type=HoverTool))
            hover.tooltips = [("Id", "@%s" % label)]
            renderer = scatter.select_one(dict(type=GlyphRenderer))
            renderer.data_source.data[label] = dataframe[label].tolist()

        if xaxis_label:
            scatter._xaxis.axis_label = xaxis_label
        if yaxis_label:
            scatter._yaxis.axis_label = yaxis_label

        return scatter

    def heatmap(self, dataframe, y=None, x=None, values=None, width=None, height=None, palette=None,
                max_color=None, min_color=None, mid_color=None, title='Heatmap', xaxis_label=None, yaxis_label=None):

        palette = self.__default_options__.get('palette', None) if palette is None else palette
        width = self.__default_options__.get('width', None) if width is None else width

        width, height = self._width_height(width, height)

        heatmap = HeatMap(dataframe, x=x, y=y, values=values, width=width, height=height, palette=palette, title=title)
        if xaxis_label:
            heatmap._xaxis.axis_label = xaxis_label
        if yaxis_label:
            heatmap._yaxis.axis_label = yaxis_label

        return heatmap

    def line(self, dataframe, x=None, y=None, width=None, height=None, groups=None, palette=None, title="Line",
             xaxis_label=None, yaxis_label=None):
        palette = self.__default_options__.get('palette', None) if palette is None else palette
        width = self.__default_options__.get('width', None) if width is None else width

        width, height = self._width_height(width, height)

        line = Line(dataframe, x=x, y=y, color=groups, width=width, height=height, palette=palette, title=title,
                    legend=True)

        if xaxis_label:
            line._xaxis.axis_label = xaxis_label
        if yaxis_label:
            line._yaxis.axis_label = yaxis_label

        return line

    def boxplot(self, dataframe, values='value', groups=None, width=None, height=None, palette=None,
                title="BoxPlot", legend=True):

        palette = self.__default_options__.get('palette', None) if palette is None else palette
        width = self.__default_options__.get('width', None) if width is None else width

        if values:
            unique_values = dataframe[groups].unique()
            palette = self._palette(palette, len(unique_values))
        else:
            palette = None

        width, height = self._width_height(width, height)

        boxplot = BoxPlot(dataframe, values=values, label=groups, color=groups, legend=True, width=width,
                          height=height, palette=palette, title=title)

        return boxplot

    @classmethod
    def display(cls, plot):
        show(plot)
