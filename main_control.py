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
import time
import argparse
import numpy as np


from NEGF_package import NEGFGlobal
from NEGF_package import GenerateHamiltonian
from NEGF_package import Alpha

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("Folder_Name",
                        help="address of the folder that contains 'main_control_parameter.yml' file")
    args = parser.parse_args()

    input_parameter = NEGFGlobal.yaml_file_loader(args.Folder_Name,
                                                  'main_control_parameter.yml')

    location = input_parameter['File_location']

    message = ['=== Program Start ===']
    NEGFGlobal.global_write(location, 'output.out', message=message)

    a = 1.42
    ax = np.sqrt(3)*a
    bz = 3*a
    NK = 1000
    NE = 100
    start_time = time.time()

    file_name = 'coordinate_' + input_parameter['system_name'] + '.ao'

    ks_matrix, overlap_matrix = NEGFGlobal.ao_file_loader(location, file_name)

#     if input_parameter['Save_file']:
#         NEGFGlobal.global_write(location, 'ks_matrix.dat', num_data=ks_matrix)
#         NEGFGlobal.global_write(location, 'overlap_matrix.dat', num_data=overlap_matrix)

    num_unit_cell = [np.array(item.split(' '), dtype='int')
                     for item in input_parameter['Num_unit_cell'].split(', ')]

    file_name = 'map_' + input_parameter['system_name'] + '.csv'
    map_file = NEGFGlobal.csv_file_loader(location, file_name)

    hamiltonian_block = []
    overlap_block = []
    for idx in range(len(num_unit_cell)):
        hamiltonian = GenerateHamiltonian.matrix_block(num_unit_cell[idx], map_file, ks_matrix)
#         hamiltonian = GenerateHamiltonian.pz_block(hamiltonian)
        hamiltonian_block.append(hamiltonian)

        overlap = GenerateHamiltonian.matrix_block(num_unit_cell[idx], map_file, overlap_matrix)
#         overlap = GenerateHamiltonian.pz_block(overlap)
        overlap_block.append(overlap)

#         file_name = 'H' + str(idx) + '.dat'
#         NEGFGlobal.global_write(location, file_name, num_data=hamiltonian_block[idx])
#
#         file_name = 'S' + str(idx) + '.dat'
#         NEGFGlobal.global_write(location, file_name, num_data=overlap_block[idx])

    file_name = 'map_coordinate_' + input_parameter['system_name'] + '.csv'
    map_coordinate_file = NEGFGlobal.csv_file_loader(location, file_name)

    center_ss = input_parameter['Center_ss']
    distance = Alpha.distance(center_ss, map_coordinate_file)
#     NEGFGlobal.global_write(location, 'distance.dat', num_data=distance)

    kx = np.linspace(-np.pi/ax, np.pi/ax, NK)
    kz = np.linspace(-np.pi/bz, np.pi/bz, NK)
    alpha_H = Alpha.alpha_cal(np.array(hamiltonian_block), kx, kz, distance)
    alpha_O = Alpha.alpha_cal(np.array(overlap_block), kx, kz, distance)

#     NEGFGlobal.global_write(location, 'alpha_H.dat', num_data=alpha_H)
#     NEGFGlobal.global_write(location, 'alpha_O.dat', num_data=alpha_O)

    energy = np.linspace(-10, 10, NE)
    energy = energy + 1j*0.001

    GR = Alpha.greens_fun(energy, alpha_H, alpha_O)
    GR_trace = Alpha.trace_cal(GR)

    NEGFGlobal.global_write(location, 'GR_trace.dat', num_data=GR_trace)
    NEGFGlobal.global_write(location, 'energy.dat', num_data=energy)

    energy = energy - (-0.06633352904223*27.2114)
    NEGFGlobal.display(location, 'DOS.png', energy, -(1/np.pi)*np.imag(GR_trace), xlabel='Energy', ylabel='DOS')
