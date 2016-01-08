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
import math
import numpy as np


def iqr(sample):
    q75, q25 = np.percentile(sample, [75, 25])
    return q75 - q25


def freedman_diaconis(sample):
    n = len(sample)
    return 2 * iqr(sample) * math.pow(n, -1/3)
