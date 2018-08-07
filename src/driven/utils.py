# -*- coding: utf-8 -*-

# Copyright (c) 2016 Novo Nordisk Foundation Center for Biosustainability,
# Technical University Denmark
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University Denmark
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Provide general utility functions."""

from __future__ import absolute_import

from depinfo import print_dependencies


__all__ = ("show_versions",)


def show_versions():
    """Print system and Python environment dependency information."""
    print_dependencies("driven")


def all_same(seq):
    """Determine whether all the elements in a sequence are the same."""
    # Compare all the elements to the first in the sequence.
    return all(elem == seq[0] for elem in seq)


def get_common_start(*seq_list):
    """Return the common prefix of a list of sequences."""
    # Map the matching elements.
    m = [all_same(seq) for seq in zip(*seq_list)]
    # Append a guard in case all the sequences match.
    m.append(False)
    # Truncate the sequence before first mismatch.
    return seq_list[0][0:m.index(False)]
