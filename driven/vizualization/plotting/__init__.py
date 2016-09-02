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

import six

from driven.vizualization.plotting.abstract import Plotter

__all__ = ["plotting"]

_engine = None
_engines = {}

try:
    from driven.vizualization.plotting.with_bokeh import BokehPlotter
    _engine = BokehPlotter()
    _engines["bokeh"] = BokehPlotter
except ImportError:
    pass

try:
    from driven.vizualization.plotting.with_ggplot import GGPlotPlotter
    _engine = GGPlotPlotter() if _engine is None else _engine
    _engines["ggplot"] = GGPlotPlotter
except (ImportError, RuntimeError):
    pass


class _plotting:
    def __init__(self, engine):
        self.__dict__['_engine'] = engine

    def __getattr__(self, item):
        if item not in ["_engine", "engine"]:
            return getattr(self.__dict__['_engine'], item)
        else:
            return self.__dict__['_engine']

    def __setattr__(self, key, item):
        if key not in ["_engine", "engine"]:
            raise KeyError(key)
        else:
            if isinstance(item, six.string_types):
                item = _engines[item]()

            if not isinstance(item, Plotter):
                raise AssertionError("Invalid engine %s" % item)

            self.__dict__['_engine'] = item


plotting = _plotting(_engine)
