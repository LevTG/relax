#!/usr/bin/python3 -u
### Extraction
#
import h5py

import numpy as np

from hdf5_API import hdfAPI
from argparse import ArgumentParser


def reform(file, output, tcf = '', gname = ''):
    ownfile = False
    out = h5py.File(output, 'w')

    if not hasattr(file, 'read'):
        ownfile = True
        file = hdfAPI(file, 'r')

    tcfs = tcf if tcf else file.get_tcfList()
    groups = gname if gname else file.get_groupList()

    out.create_dataset('time', data=file.get_time())
    for name in groups:
        group = out.create_group(name)

        for tcf in tcfs:
            if not file._get_zeroGroup(tcf, name):
                continue
            gtcf = group.create_group(tcf)
            gtcf.attrs['group_size'] = file._get_zeroGroup(tcf, name)['group_size'][()]
            gtcf.attrs['names'] = file.get_names(tcf, name)

            mean, std = file.mean_tcf(tcf, name)
            gtcf.create_dataset("mean", data=mean)
            gtcf.create_dataset("errs", data=std)

    if ownfile:
        file.close()
    out.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("filename", type=str)
    parser.add_argument('-o', '--output', required=False, type=str, default='out')
    parser.add_argument('--tcf', required=False, nargs='*', default='')
    parser.add_argument('-g', '--gname', required=False, nargs='*', default='')
    args = parser.parse_args()
    reform(args.filename, args.output, args.tcf, args.gname)

if __name__ == '__main__':
    main()