#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hdfplot.py
#
#  Copyright 2020 Svetlana Nolde <sveta@elk.localdomain>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import h5py
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser

class HdfPlotError(Exception):
    """docstring for HdfPlotError"""
    def __init__(self, msg):
        super(HdfPlotError, self).__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


def param_autodetect(params):
    coeff = params[0]
    ampl = params[1::2]
    tau  = params[2::2]
    return tau, ampl, coeff, None


def expall(param, arrtime):
    tau, ampl, coeff, tcf_type = param_autodetect(param)
    out = np.copy(coeff)
    for i in range(tau.shape[0]):
        out = out + ampl[i]*np.exp(arrtime/(-tau[i]))
    return out


def main(args):
    parser = ArgumentParser()
    parser.add_argument("-d", '--data-file', required=True, type=str, help='file which content original data')
    parser.add_argument("-f", '--fit-file', required=True, type=str, help='file which content fit result')
    parser.add_argument('-g', '--group', default='NH', help='Which group you want to plot')
    parser.add_argument('--tcf', default='acf', help='Type of correlation function. Can be *acf* or *ccf*')
    parser.add_argument('--timelines', action='store_true', help='Show grid line for each time value')
    parser.add_argument('--chisqr', action='store_true', help='Shoe grid line for each time value')
    parser.add_argument('--exp-start', default=4, type=int, help='Number of exponents to start from')
    parser.add_argument('--exp-finish', default=6, type=int, help='Number of exponents when finish')
    parser.add_argument('--istart', default=0, type=int, help='index to start plot from')
    parser.add_argument('--ifinish', default=float('inf'), type=int, help='index to finish plot at')
    args = parser.parse_args()

    with h5py.File(args.data_file, 'r') as fd:
        time = fd['time'][:]
        print("Time line: {} ({}..{} ps)".format(time.shape, time[0], time[-1]))
        if not args.group in fd.keys():
            raise HdfPlotError("Group {} not found in file {}".format(args.group, args.fit_file))
        elif not args.tcf in fd[args.group].keys():
            raise HdfPlotError("Tcf {} not found in file {}".format(args.tcf, args.fit_file))
        mean = fd[args.group][args.tcf]['mean'][:]
        errs = fd[args.group][args.tcf]['errs'][:]
        #print("Data shape: {}".format(mean.shape))
        names = fd[args.group][args.tcf].attrs['names']
        names = names.splitlines()

    params = {}
    stats = {}
    taus = {}
    ampls = {}
    with h5py.File(args.fit_file, 'r') as fd:
        group_list = [g for g in fd]
        if not args.group in group_list:
            raise HdfPlotError("Group {} not found in file {}".format(args.group, args.fit_file))
        item_list = [item for item in fd[args.group]]
        if args.tcf in item_list:
            fd_exps = fd[args.group][args.tcf]
        else:
            # Backward compatibility - no 'acf/ccf' was in previous versions
            fd_exps = fd[args.group]
        exps_range = range(args.exp_start, args.exp_finish + 1)
        exps_args = ['exp{}'.format(e) for e in exps_range]
        exps_infile = [exp for exp in fd_exps if exp in exps_args]
        if not exps_infile:
            raise HdfPlotError("Exponents '{}' not found in {}/{}".format(" ".join(exps_args),
                                                             args.fit_file, args.group))
        print("exps in file: {}".format(exps_infile))
        for exp in exps_infile:
            params[exp] = fd_exps[exp]['params'][:]
            stats[exp] = fd_exps[exp]['stats'][:]
            taus[exp], ampls[exp], coeff, tcf_type = param_autodetect(params[exp][0, :])

            # print("params[{}].shape: {}".format(exp, params[exp].shape))
            # print("stats[{}].shape: {}".format(exp, stats[exp].shape))
        # exp4=np.array(fd.get('{}/exp4/params'.format(args.group)))
        # exp5=np.array(fd.get('{}/exp5/params'.format(args.group)))
        # exp6=np.array(fd.get('{}/exp6/params'.format(args.group)))
        # stat4 = np.array(fd.get('{}/exp4/stats'.format(args.group)))
        # stat5 = np.array(fd.get('{}/exp5/stats'.format(args.group)))
        # stat6 = np.array(fd.get('{}/exp6/stats'.format(args.group)))

    # tau4 = exp4[:,1::2]
    # tau5 = exp5[:,1::2]
    # tau6 = exp6[:,1::2]
    #chisqr4 = stat4[:][2]
    #chisqr5 = stat5[:][2]
    #chisqr6 = stat6[:][2]

    # print("exp4.shape={}".format(exp4.shape))
    # print("tau4.shape={}".format(tau4.shape))

    start = max(0, args.istart)
    finish = min(mean.shape[0], args.ifinish)
    for res in range(start, finish):
        labels = {}
        for exp in exps_infile:
            if args.chisqr:
                #stats_exp = stats[exp]
                # print("stats_exp: {}".format(stats_exp))
                #stats_exp_res = stats_exp[res]
                # print("stats_exp_res: {}".format(stats_exp_res))
                labels[exp] = '{} {:>8.2f}'.format(exp, stats[exp][res][2])
            else:
                labels[exp] = exp
            # label4 = '4 exps {:>7.2f}'.format(stat4[res][2])
            # label5 = '5 exps {:>7.2f}'.format(stat5[res][2])
            # label6 = '6 exps {:>7.2f}'.format(stat6[res][2])
        # else:
        #     (label4, label5, label6) = ('4 exps', '5 exps', '6 exps')
        print("labels: {}".format(labels))

        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(names[res])
        ax[0].errorbar((time+1)[:7], mean[res, :7],
                       yerr=errs[res, :7], fmt=' ', zorder=0, color='darkgrey')
        for xsel in (slice(7, 100, None), slice(100, 1000, 5), slice(1000, 2000, 10),
                     slice(2000, 4000, 50), slice(4000, 8000, 200),
                     slice(8000, None, 1000)):
            ax[0].errorbar((time+1)[xsel], mean[res, xsel],
                           yerr=errs[res, xsel], fmt=' ', zorder=0, color='C0')
        ax[0].set_xscale('log')

        line_styles_list = ['dashdot', 'dashed', 'dotted']
        colors, styles, c = {}, {}, 0
        for exp in exps_infile:
            colors[exp] = 'C{}'.format(c+1)
            styles[exp] = line_styles_list[c % len(line_styles_list)]
            c += 1
            if max(params[exp][res, :]) > 0:
                ax[0].plot(time + 1, expall(params[exp][res], time),
                           label=labels[exp], zorder=10, color=colors[exp])
        # for expdata, explabel, col in zip((exp4, exp5, exp6),
        #                                   (label4, label5, label6), ('C1', 'C2', 'C3')):
        #     if(expdata[res,-1]>0):
        #         ax[0].plot(time+1, expall(expdata[res],time),
        #                    label=explabel, zorder=10, color=col)


        ntime = 200
        if args.timelines:
            for exp in exps_infile:
            # for tau, col, style in zip(
            #         (tau4[res], tau5[res], tau6[res]),
            #         ('C1', 'C2', 'C3'),
            #         ('dashdot', 'dashed', 'dotted')):
                for t in taus[exp][res]:
                    if t < time[-1]:
                        ax[0].axvline(t, linestyle=styles[exp], color=colors[exp])
                    if t < time[ntime+1]:
                        ax[1].axvline(t, linestyle=styles[exp], color=colors[exp])

        ax[0].legend()


        ax[1].errorbar(time[1:7]+1, mean[res, 1:7],
                       yerr=errs[res, 1:7], fmt=' ', zorder=0, color='darkgrey')
        ax[1].errorbar(time[7:ntime]+1, mean[res, 7:ntime],
                       yerr=errs[res, 7:ntime], fmt=' ', zorder=0)
        ax[1].set_xscale('log')
        ylim = ax[1].get_ylim()
        for exp in exps_infile:
        # for expdata, explabel, col  in zip((exp4, exp5, exp6),
        #                                    (label4, label5, label6), ('C1', 'C2', 'C3')):
            if max(params[exp][res, :]) > 0:
                ax[1].plot((time[1:ntime]+1),
                           expall(params[exp][res, :ntime],
                                  time[:ntime])[1:],
                           label=labels[exp], zorder=10, color=colors[exp])
        ax[1].set_ylim(ylim)
        ax[1].legend()
        plt.show()
        plt.close()
    return 0


if __name__ == '__main__':
    import sys
    try:
        sys.exit(main(sys.argv))
    except HdfPlotError as e:
        print(e, file=sys.stderr)
        exit(-1)
