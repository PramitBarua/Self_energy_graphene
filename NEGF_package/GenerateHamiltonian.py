#!/usr/bin/env python3

__author__ = "Pramit Barua"
__copyright__ = "Copyright 2018, INT, KIT"
__credits__ = ["Pramit Barua"]
__license__ = "INT, KIT"
__version__ = "1"
__maintainer__ = "Pramit Barua"
__email__ = ["pramit.barua@student.kit.edu", "pramit.barua@gmail.com"]

'''
methods:
    matrix_slice: returns Hn1, H00, Hp1 matrixes from KS matrix
            or    returns sn1, s00, sp1 matrixes from overlap matrix

    matrix_block: 

'''

import numpy as np


def matrix_slice(matrix, matrix_name, **kargs):
    '''
    matrix_slice(matrix, matrix_name, num_unit=3)
        returns Hn1, H00, Hp1 matrixes from KS matrix

        Parameter:
            matrix: KS Matrix or Overlap Matrix
            matrix_name: 'KS_matrix' or 'Overlap_matrix'
                this defines which matrix the method is dealing
            num_unit: int, optional
                this will divide the matrix into NxN matrix
                default is 3
    '''

    def matrix_split(matrix, num_unit):
        matrix_shape = matrix.shape
        if (matrix.shape[0] % num_unit) == 0 and (matrix.shape[1] % num_unit) == 0:
            A = matrix[:int(matrix_shape[0]/num_unit),
                       :int(matrix_shape[0]/num_unit)]
            K = matrix[int(matrix_shape[0]*2/num_unit):,
                       int(matrix_shape[0]*2/num_unit):]
            up_middle_matrix = matrix[:int(matrix_shape[0]/num_unit),
                                      int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit)]
            left_matrix = matrix[int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit),
                                 :int(matrix_shape[0]/num_unit)]
            middle_matrix = matrix[int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit),
                                   int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit)]
            right_matrix = matrix[int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit),
                                  int(matrix_shape[0]*2/num_unit):]
            low_middle_matrix = matrix[int(matrix_shape[0]*2/num_unit):,
                                       int(matrix_shape[0]/num_unit):int(matrix_shape[0]*2/num_unit)] 

            if np.all(up_middle_matrix == np.matrix.getH(left_matrix)):
                print('B and D+ is the same')
            else:
                print('WARNNING: B and D+ is not the same')

            if np.all(right_matrix == np.matrix.getH(low_middle_matrix)):
                print('F and H+ is the same')
            else:
                print('WARNNING: F and H+ is not the same')

            if np.all(middle_matrix == A) and np.all(middle_matrix == K):
                print('Diagonal elements are the same')
            else:
                print('WARNNING: Diagonal elements are not the same')

            status = True
            return left_matrix, middle_matrix, right_matrix, status
        else:
            left_matrix = np.array([])
            middle_matrix = np.array([])
            right_matrix = np.array([])
            status = False
            return left_matrix, middle_matrix, right_matrix, status

    if matrix_name == 'KS_matrix':
        variable_name = ['hn1', 'h00', 'hp1']
    elif matrix_name == 'Overlap_matrix':
        variable_name = ['sn1', 's00', 'sp1']

    if 'num_unit' in kargs:
        num_unit = kargs['num_unit']
    else:
        num_unit = 3

    if matrix.size != 0:
        hn1, h00, hp1, status = matrix_split(matrix, num_unit)
        if status:
            print('=== The matrix is sliced successfully ===')
        else:
            print('=== Error: Matrix dimension is not divisible by ',
                  str(num_unit), ' ===')

        return hn1, h00, hp1


def matrix_block(num_unit_cell, map_file, dat_file, **kargs):
    '''
    num_unit_cell: numpy array
        which super-cells have to be extracted
        This argument should not contain more than 2 'int'
        If both of the 'int' are the same then the code will
        produce a matrix of the following form
            [[a, 0], [0, a]]
    num_orbital: int, optioal
        number of orbital of the system
    map_file: number of atom of SS (which atom belongs to which SS, the map of the system)
    dat_file: KS Matrix or Overlap matrix
    '''

    def atom_index(mapping, unit_cells):
        # this method accumulates the atom index 
        # from unit cells
        index = []
        for atom in mapping[unit_cells]:
            index.append(int(atom))
        return index

    def extract(value, num_orbital, extract_form):
        # 'extract_form' is a numpy ndarray and this contains
        # all the values (from where to extract)
        # value is an array that contains the number of
        # atom we are interested in (what to extract)
        # This method extracts value from 'extract_form' based
        # on the value and num_orbital
        extract_value = []
        for target_atom in value[0]:
            row_num = ((target_atom-1)*num_orbital)
            for idx1 in range(num_orbital):
                buffer = []
#                 print('row_num = ', row_num)
                for atom in value[1]:
                    column_num = ((atom-1)*num_orbital)
#                     print('column_num = ', column_num)
                    buffer.extend(extract_form[row_num+idx1, column_num:column_num+num_orbital])
                extract_value.append(np.array(buffer))
        return np.array(extract_value)

#     def zero_append(array):
#         array_size = 2
#         return np.kron(np.eye(array_size), array)

    if 'num_orbital' in kargs:
        num_orbital = kargs['num_orbital']
    else:
        num_orbital = 4

#     if num_unit_cell[0] == num_unit_cell[1]:
#         num_unit_cell = np.array([num_unit_cell[0]])

    all_atom = []
    for index in num_unit_cell:
        all_atom.append(atom_index(map_file, index))
    extract_value = extract(all_atom, num_orbital, dat_file)

#     if len(num_unit_cell) == 1:
#         extract_value = zero_append(extract_value)
    return extract_value


def pz_block(array, **kargs):
    if 'num_orbital' in kargs:
        num_orbital = kargs['num_orbital']
    else:
        num_orbital = 4

    array_shape = array.shape
    result_shape = (int(array_shape[0]/num_orbital), int(array_shape[1]/num_orbital))
    result = np.zeros(result_shape)

    count = 0
    for idx in range(3, array_shape[0], num_orbital):
        result[count] = array[idx][num_orbital-1::num_orbital]
        count += 1
    return result