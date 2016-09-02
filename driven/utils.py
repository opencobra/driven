# Copyright 2016 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


def all_same(seq):
    """
    Determines whether all the elements in a sequence are the same.

    seq: list
    """

    # Compare all the elements to the first in the sequence,
    # then do a logical and (min) to see if they all matched.
    return min([elem == seq[0] for elem in seq]+[True])


def get_common_start(*seq_list):
    """
    Returns the common prefix of a list of sequences.

    Raises
        an exception if the list is empty.
    """
    m = [all_same(seq) for seq in zip(*seq_list)]  # Map the matching elements
    m.append(False)                                # Guard in case all the sequences match
    return seq_list[0][0:m.index(False)]           # Truncate before first mismatch
