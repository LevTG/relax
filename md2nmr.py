#!/usr/bin/python3 -u

import os
import sys
import time
import copy
import subprocess

import numpy as np

import _pickle as pickle
import multiprocessing.dummy as mp

from arg_parse import read_args, check_args
# from utils import join_lists

from classes.config import Config
from classes.topology import Topology
from classes.mdtrace import Mdtrace
from classes.logger import create_logger
from classes.relaxation_groups import Relaxation_Group_Types_Names

# test_dir = os.path.join('.', 'src/')
# sys.path.insert(0, test_dir)

script_path = os.path.abspath(__file__)
script_dir = os.path.split(script_path)[0]
cfg_file = os.path.join(script_dir, './md2nmr.cfg')


def main():
    args = read_args(cfg_file)
    logger = create_logger(args)
    if not check_args(args, logger):
        exit()

    top = args.topology
    traces = args.mdtrace
    md_count = len(traces)

    # groups = input('Please input groups which you want to acf or ccf should be calculated: ').split(sep=' ')

    # if 'all' in args.group_types:
    #     groups = Relaxation_Group_Types_Names
    # else:
    #     groups = args.group_types

    groups = Relaxation_Group_Types_Names


    try:
        if args.load_file:
            print('Loading...')
            ss = time.time()
            mol = Topology()
            mol = mol.load(args.load_file)
            print(time.time() - ss)
        else:
            config = Config()
            config.read_config(cfg_file)
            mol = Topology()
            mol.read_topology(config, top)
            mol.read_coords()

    except Exception as e:
        print("Undefined exception while loading: {} {}".format(type(e), e))
        raise

    cmol = mol.molecules[0]

    calc_groups = cmol.get_calc_groups(groups)
    print(len(calc_groups))
    # calc_groups = calc_groups

    cmol.calculate_accf(calc_groups)
    cmol.save_accfs(calc_groups)
    # cmol.load_accfs(calc_groups)

    # for group in calc_groups:
        # if group.simple_type in groups:
            # group.trace.plot_accf()

#    try:
#        for i in range(2, 6):
#           cmol.fit(calc_groups, i)
#
#            fd = open('exppng/exp{}.txt'.format(i), 'w')
#            cmol.output(calc_groups, fd)
#            fd.close()
#
#            for group in calc_groups:
#                if group.simple_type in groups:
#                    if group.trace.fit and group.trace.fmodel.model_fit.covar is not None:
#                        group.trace.plot_fit()
#    except:
#        print('Smth get wrong')
#        raise

#    fd = open('exp2.txt', 'w')

#    cmol.output(calc_groups, fd)
#    fd.close()

#    for group in calc_groups:
#        if group.simple_type in groups:
#            if group.trace.fit and group.trace.fit.model_fit.success:
#                group.trace.plot_fit()

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print("Undefined exception while main: {} {}".format(type(e), e))
        raise
