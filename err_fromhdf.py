#!/usr/bin/python3 -u

import h5py
import numpy as np

from argparse import ArgumentParser

from fitter.fit_res import REFERENCE_PARAMS_SEQ


def main():
    parser = ArgumentParser()
    parser.add_argument("filename", type=str)
    parser.add_argument('-i', '--idata', type=int, help='index in data array for fit, e.g. 0-1')
    parser.add_argument('-s', '--exp-start', default=4, type=int, help='Number of exponents to start from')
    parser.add_argument('-f', '--exp-finish', default=6, type=int, help='Number of exponents when finish')
    parser.add_argument('-g', '--group', default='NH', help='Which group you want to fit. Need to fit data from hdf')
    parser.add_argument('--tcf', default='acf', help='Need to fit data from hdf')
    parser.add_argument('-o', '--output', default='out.hdf', help='filename for saving results')
    args = parser.parse_args()

    print("Parameters and errors in HDF file {}".format(args.filename))
    fid = h5py.File(args.filename, 'r')
    tcf = fid[args.group][args.tcf]
    for exp in range(args.exp_start, args.exp_finish + 1):
        print("data {}, exp{}:".format(args.idata,exp))
        values = tcf['exp{}'.format(exp)]['params'][args.idata]
        covar = tcf['exp{}'.format(exp)]['covar'][args.idata]
        errors = np.sqrt(np.diag(covar))

        for name, val, err in zip(REFERENCE_PARAMS_SEQ, values, errors):
            print("{:11} : {:13.6f} +/- {:13.6f}".format(name, val, err))
        print('\n')

if __name__ == '__main__':
    main()

