#!/usr/bin/python3
import os
import subprocess

import argparse as argp


def read_args(cfg_file):
    parser = argp.ArgumentParser(description="")

    parser.add_argument('-c', '--config', default=cfg_file, help='File containing definition for relaxation groups and additional information (like coordinates, formulas, etc.)')
    parser.add_argument('-t', '--topol', default=None, type=str, required=True, dest='topology', help='Gromacs generated topology file')
    parser.add_argument('-m', '--mdtr', default=[], type=str, dest='mdtrace', nargs='*', help='MD trace in .pdb or .trr format')
    parser.add_argument('-o', default='./savings/checkpoint.npy', type=str, dest='save_file', help='Name of out file')
    parser.add_argument('-l', '--load', default=None, type=str, dest='load_file')
    parser.add_argument('-name', default='md2nmr', dest='name', help='all files generataed with this program will have this name')
    # parser.add_argument('-g', '--group-types', default='all', type=str, dest='group_types', nargs='*', choices= ['all', 'NH', 'NH2', 'CHAl', 'CHAr', 'CH2', 'CH3'], help='Relaxation group types for calculation (for ex. NH CH)')
    return parser.parse_args()


def check_args(ns, logger):
    for name, files in vars(ns).items():
        if type(files) is list:
            for file in files:
                if not os.path.exists(file):
                    print("ERROR! File {} not found!".format(file))
                    return False
        elif files != None and name != 'save_file' and name != 'name':
            if not os.path.exists(files):
                print("ERROR! File {} not found!".format(files))
                return False

    if ns.mdtrace != None:
        for file in ns.mdtrace:
            filename, file_extension = os.path.splitext(file)
            if file_extension != '.pdb':
                pdb_file = filename + '.pdb'
                print('Convert {} to .pdb format'.format(file))
                subprocess.run( '/usr/local/gromacs/bin/gmx trjconv -f {} -o {} -s ./topol.tpr'.format(file, pdb_file))
                file = pdb_file

    return True
