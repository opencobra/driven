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


from functools import partial
from cameo import pfba
from cameo.util import TimeMachine


def fba(model, objective=None, distribution=None, *args, **kwargs):
    with TimeMachine() as tm:
        for reaction_id in distribution:
            reaction = model.reactions.get_by_id(reaction_id)
            tm(do=partial(setattr, reaction, "upper_bound", distribution[reaction_id][1]),
               undo=partial(setattr, reaction, "upper_bound", reaction.upper_bound))
            tm(do=partial(setattr, reaction, "lower_bound", distribution[reaction_id][0]),
               undo=partial(setattr, reaction, "lower_bound", reaction.lower_bound))
        return pfba(model, objective, *args, **kwargs)