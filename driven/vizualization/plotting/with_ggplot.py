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
from matplotlib.ticker import Locator

from driven.vizualization.plotting.abstract import Plotter
from ggplot import ggplot, aes
from ggplot import geom_histogram, geom_tile, geom_point
from ggplot import scale_color_brewer, scale_colour_manual
from ggplot import ggtitle
from ggplot import scale_x_continuous, scale_y_continuous, scale_colour_gradient2

Locator.MAXTICKS = 15000

class gradient:
    def __init__(self, low, mid, high):
        self.low = low
        self.mid = mid
        self.high = high


class GGPlotPlotter(Plotter):
    def __init__(self, **defaults):
        super(GGPlotPlotter, self).__init__(**defaults)

    def histogram(self, dataframe, bins=100, width=None, height=None, palette=None, title='Histogram', values=None,
                  groups=None, legend=True):
        palette = self.__default_options__.get('palette', None) if palette is None else palette

        return ggplot(dataframe, aes(x=values, fill=groups, color=groups)) + \
               geom_histogram(alpha=0.6, breaks=bins, position="fill") + \
               self._palette(palette) + \
               ggtitle(title) + \
               scale_y_continuous(name="Count (%s)" % values)

    def scatter(self, dataframe, x=None, y=None, width=None, height=None, color=None, title='Scatter', xaxis_label=None,
                yaxis_label=None, label=None):
        color = self.__default_options__.get('palette', None) if color is None else color
        width = self.__default_options__.get('width', None) if width is None else width

        gg = ggplot(dataframe, aes(x, y)) + geom_point(color=color, alpha=0.6) + ggtitle(title)
        if xaxis_label:
            gg += scale_x_continuous(name=xaxis_label)
        if yaxis_label:
            gg += scale_y_continuous(name=xaxis_label)

        return gg

    def heatmap(self, dataframe, y=None, x=None, values=None, width=None, height=None,
                max_color=None, min_color=None, mid_color=None, title='Heatmap'):
        max_color = self.__default_options__.get('max_color', None) if max_color is None else max_color
        min_color = self.__default_options__.get('min_color', None) if min_color is None else min_color
        mid_color = self.__default_options__.get('mid_color', None) if mid_color is None else mid_color
        width = self.__default_options__.get('width', None) if width is None else width

        palette = gradient(min_color, mid_color, max_color)
        return ggplot(dataframe, aes(x=x, y=y, fill=values)) + \
               geom_tile() + \
               self._palette(palette, "div")

    def _palette(self, palette, type="seq", **kwargs):
        if isinstance(palette, six.string_types):
            return scale_color_brewer(type=type, palette=palette)
        elif isinstance(palette, gradient):
            return scale_colour_gradient2(low=palette.low, mid=palette.mid, high=palette.high)
        elif isinstance(palette, collections.Iterable):
            return scale_colour_manual(values=palette)

    @classmethod
    def display(cls, plot):
        pass