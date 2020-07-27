import h5py
import numpy as np

from argparse import ArgumentParser

PARAMS_SEQ = ['c', 'aamplitude', 'adecay', 'bamplitude', 'bdecay', 'camplitude', 'cdecay',
            'damplitude', 'ddecay', 'eamplitude', 'edecay', 'famplitude', 'fdecay',
            'gamplitude', 'gdecay','hamplitude', 'hdecay','iamplitude', 'idecay']

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

    fid = h5py.File(args.filename, 'r')
    tcf = fid[args.group][args.tcf]
    for exp in range(args.exp_start, args.exp_finish + 1):
        cov = tcf['exp{}'.format(exp)]['covar'][args.idata]
        err = np.sqrt(np.diag(cov))
        for err, par in zip(err, PARAMS_SEQ):
            print("{:11} : {:.4f}".format(par, err))
        print('\n')

if __name__ == '__main__':
    main()

