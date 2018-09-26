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
    len_kx = len(kx)
    len_kz = len(kz)
    H_shape = hamiltonian[0].shape

    result = np.zeros((len_kx, len_kz, H_shape[0], H_shape[1]), dtype=complex)

    for idxKx, itemKx in enumerate(kx):
        for idxKz, itemKz in enumerate(kz):
            for idxH, itemH in enumerate(hamiltonian):
                kr = itemKx*distance[idxH][0] + itemKz*distance[idxH][2]
                result[idxKx][idxKz] += itemH*np.exp(1j*kr)

    return result


def greens_fun(energy, hamiltonian, overlap):
    H_shape = hamiltonian.shape
    result_shape = energy.shape + (H_shape[2], H_shape[3])

    result = np.zeros(result_shape, dtype=complex)

    for idE, energy_item in enumerate(energy):
        value_inv = np.zeros(H_shape, dtype=complex)
        for idxKx in range(H_shape[0]):
            for idxKz in range(H_shape[1]):
                value = (energy_item*overlap[idxKx][idxKz]) - hamiltonian[idxKx][idxKz]
                value_inv[idxKx][idxKz] = np.linalg.inv(value)
        sum_value = sum_cal(value_inv)
        result[idE] = sum_value/(H_shape[0]*H_shape[1])
    return result


def sum_cal(array):
    shape_array = array.shape
    value = 0

    for idx in range(shape_array[0]):
        for idy in range(shape_array[1]):
            value += array[idx][idy]

    return value


def trace_cal(array):
    shape_array = array.shape
    result = np.zeros(shape_array[0], dtype=complex)
    for idx in range(len(array)):
        result[idx] = np.trace(array[idx])
    return result
