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

import os
import numpy as np
import yaml


def global_write(location, file_name, **kargs):
    r''' read me:
    This is the global function to write data into file
    ideally this function should handle any sort of data
    this function can write any type of file (.dat, .out, .xyz)

    Data type can be:
        1. purely numeric
        2. purely string

    arguments:
        location: where data file will be saved (directory address)
        file_name: name of the file where data will be stored
        *kargs: it contains 2 types of key
            num_data: if data is purely numerical (numpy array)
            message: if data is string (wrapped in a list)
                message is in append mode which means if the file exists
                then program will append data instead of over writing.

    note: in any situation if data contains number as well as
    string then convert data into string and wrapped each line
    in a list. However, this wrapping should be done outside of
    this function
    '''
    if not os.path.isdir(location):
        os.makedirs(location)

    if 'num_data' in kargs:
        np.savetxt(os.path.join(location, file_name),
                   kargs['num_data'].view(float), delimiter=' ')

    if 'message' in kargs:
        with open(os.path.join(location, file_name), "a+") as f:
            for line in kargs['message']:
                f.write(str(line) + "\n")


def ao_file_loader(input_location, file_name):
    r''' read me
    this method takes ao file and return overlap matrix and kohn-sham 
    matrix into numpy array

    Note: directory should not contain more than one ao file

    input argument:
        input_location: ao file directory
        file_name: name of the ao file
    output:
        kohn-sham matrix
        overlap matrix
    '''

    overlap_matrix = []
    ks_matrix = []
    file = os.path.join(input_location, file_name)

    if not os.path.isfile(file):
        message = (
            '=== Directory path is invalid or '
            + file_name
            + ' does not exist in the directory ===')
        global_write(input_location, 'output.out', message=[message])

    else:
        message = ('=== ' + file_name + ' is found ===')
        global_write(input_location, 'output.out', message=[message])

        overlap_matrix_visited = False
        ks_matrix_visited = False
        overlap_matrix_write = False
        ks_matrix_write = False

        with open(file, "r") as raw_data:
            for line in raw_data:
                if line == ' OVERLAP MATRIX\n':
                    if not overlap_matrix_visited:
                        overlap_matrix_visited = True
                        overlap_matrix_write = True
                        ks_matrix_write = False
                    elif ks_matrix_visited:
                        break
                elif line == ' KOHN-SHAM MATRIX\n':
                    if not ks_matrix_visited:
                        ks_matrix_visited = True
                        ks_matrix_write = True
                        overlap_matrix_write = False
                    elif overlap_matrix_visited:
                        break

                if overlap_matrix_visited and overlap_matrix_write:
                    if line == '\n':
                        pass
                    else:
                        line = line.split()
                        if len(line) > 4:
                            data_value = []
                            for item in line[4:]:
                                data_value.append(float(item))
                            try:
                                overlap_matrix[int(line[0])-1].extend(data_value)
                            except IndexError:
                                overlap_matrix.append(data_value)

                elif ks_matrix_visited and ks_matrix_write:
                    if line == '\n':
                        pass
                    else:
                        line = line.split()
                        if len(line) > 4:
                            data_value = []
                            for item in line[4:]:
                                data_value.append(float(item)*27.2114)
                            try:
                                ks_matrix[int(line[0])-1].extend(data_value)
                            except IndexError:
                                ks_matrix.append(data_value)

    return np.array(ks_matrix), np.array(overlap_matrix)


def csv_file_loader(location, file_name):
    r''' read me
    this method returns atom name and coordinate of atoms from xyz file
    atom name is in list and coordinate is in numpy array
    input argument:
        location:
            directory of where the csv file is
            example: C:\Users\Desktop
        file_name:
            the file name
            example: nt-4-0-3.csv

    output/return:
        contant
    '''

    contant = []
    file = os.path.join(location, file_name)
    if not os.path.isfile(file):
        message = (
            '=== Directory path is invalid or' 
            + file_name
            + ' does not exist in the directory ===')
        global_write(location, 'output.out', message=[message])
    else:
        message = ('=== ' + file_name + ' is found ===')
        global_write(location, 'output.out', message=[message])

        with open(file, "r+") as file_contant:
            # skip first 2 line of the file
            for i in range(2):
                file_contant.readline()
            for line in file_contant:
                line_contant = line.split()
#                 contant.append(line_contant[1:])
                contant.append(np.array([float(item) for item in line_contant[1:]]))

    return contant


def xyz_file_loader(location, file_name):
    r''' read me
    this method returns atom name and coordinate of atoms from xyz file
    atom name is in list and coordinate is in numpy array
    input argument:
        location: 
            directory of where the xyz file is
            example: C:\Users\Desktop
        file_name:
            the file name
            example: nt-4-0-3.xyz

    output/return:
        atom name
        x, y and z coordinates
    '''
    atoms = []
    coordinates = []

    file = os.path.join(location, file_name)

    if not os.path.isfile(file):
        message = (
            '=== Directory path is invalid or '
            + file_name
            + ' does not exist in the directory ===')
        global_write(location, 'output.out', message=[message])

    else:
        message = ('=== ' + file_name + ' is found ===')
        global_write(location, 'output.out', message=[message])

        with open(file, "r+") as file_contant:
            n_atoms = int(file_contant.readline())
            title = file_contant.readline()
            for line in file_contant:
                try:
                    atom, x, y, z = line.split()
                    atoms.append(atom)
                    coordinates.append([float(x), float(y), float(z)])
                except ValueError:
                    pass

        if n_atoms != len(coordinates):
            message = (
                '=== Number of atom in the file does not match'
                + ' with the number of atom mentioned in the file ===')
            global_write(location, 'output.out', message=[message])
            atoms = []
            coordinates = []

    return atoms, np.array(coordinates)


def yaml_file_loader(location, file_name):
    ''' read me
    '''
    file = os.path.join(location, file_name)
    if not os.path.isfile(file):
        message = (
            '=== Directory path is invalid or '
            + file_name
            + ' does not exist in the directory ===')
        global_write(location, 'output.out', message=[message])
        return {}
    else:
        message = ('=== ' + file_name + ' is found ===')
        global_write(location, 'output.out', message=[message])
        if os.path.isfile(file):
            with open(file, "r") as raw_data:
                data = yaml.load(raw_data)
                return data['Input']



def display(location, file_name, varX, varY, **kargs):
    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    fig_location = os.path.join(location, 'fig')
    if not os.path.isdir(fig_location):
        os.makedirs(fig_location)

    fig1 = plt.figure(1)

    if 'label' in kargs:
        plt.plot(varX, varY, label=kargs['label'])
        plt.legend()
    else:
        plt.plot(varX, varY)

    if 'xlabel' in kargs:
        plt.xlabel(kargs['xlabel'])

    if 'ylabel' in kargs:
        plt.ylabel(kargs['ylabel'])

    if 'title' in kargs:
        plt.title(kargs['title'])

    plt.grid()
    fig1.savefig(os.path.join(fig_location, file_name))
    message = ('=== Figure saved ===')
    global_write(location, 'output.out', message=[message])
