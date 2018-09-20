#!/usr/bin/env python3

__author__ = "Pramit Barua"
__copyright__ = "Copyright 2018, INT, KIT"
__credits__ = ["Pramit Barua"]
__license__ = "INT, KIT"
__version__ = "1"
__maintainer__ = "Pramit Barua"
__email__ = ["pramit.barua@student.kit.edu", "pramit.barua@gmail.com"]

'''

'''

import numpy as np


def distance(center_ss, coordinate):
    unit_cell = coordinate[center_ss]
    result = coordinate - unit_cell
    return result


def alpha_cal(hamiltonian, kx, kz, distance):
    result = []
    for itemKx in kx:
        buffer = []
        for itemKz in kz:
            sum_value = 0
            for idxH in range(len(hamiltonian)):
                kr = itemKx*distance[idxH][0] + itemKz*distance[idxH][2]
                sum_value = sum_value + hamiltonian[idxH]*np.exp(1j*kr)
            buffer.extend(sum_value)
        buffer = np.array(buffer)
        result.extend(np.transpose(buffer))
    return np.array(result)


def greens_fun(energy, hamiltonian, overlap):
    result = []
    for energy_item in energy:
        result.append(np.linalg.inv((energy_item*overlap)-hamiltonian))
    return np.array(result)


def trace_cal(array):
    result = []
    for idx in range(len(array)):
        result.append(np.trace(array[idx]))
    return np.array(result)
