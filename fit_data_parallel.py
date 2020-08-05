#!/usr/bin/python3 -u
import sys
import h5py
import copy
from threadpoolctl import threadpool_limits

import time as t
import numpy as np
import utils as u

from fitter.fitting import Fitter
from counter import Counter

from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count

STAT_PARAMS_NAMES = ('aic', 'bic', 'chisqr', 'redchi')
NCPU = cpu_count()

#os.system("taskset -p 0xff %d" % os.getpid())

def load_data(args, group):
    (time, func, errs, names) = (None, None, None, None)
    if args.type == 'npy':
        with open(args.filename, 'rb') as fd:

            time = np.load(fd)
            func = np.load(fd)
            errs = np.load(fd)

            #errs[:, 0] = errs[:, 1]
            names = [str(i) for i in range(func.shape[0])]


    elif args.type == 'csv':
        data = np.loadtxt(args.filename, delimiter=',')
        time = data[:, 0]
        func = [ data[:, 1] ]
        errs = [ np.sqrt(data[:, 2]) ]
        names = [args.filename]
        # ОЧЕНЬ ВАЖНО!!
        errs[0] = errs[1]

    elif args.type == 'hdf':
        with h5py.File(args.filename, 'r') as fd:
            if group not in fd.keys():
                raise Exception("Wrong groupname")
            time = fd['time'][:]
            func = fd[group][args.tcf]['mean'][:]
            errs = fd[group][args.tcf]['errs'][:]
            names = fd[group][args.tcf].attrs['names']
            names = names.splitlines()
            errs[:, 0] = errs[:, 1]
    return time, func, errs, names

def prepare_data(time, data, errs, time_cut, sum_one_flag=False):
    last_point_to_del = 0
    for i, t in enumerate(time):
        if t >= time_cut:
            last_point_to_del = i
            break
    #while time[point_cut]<time_cut:
    #    point_cut += 1
    #step = time[1] - time[0]
    #space_to_del = int(time_cut // step)
    start_point_to_del = 1 if sum_one_flag else 0
    idx_del = np.arange(start_point_to_del, last_point_to_del)
    print("idx_del= {}".format(idx_del))
    time = np.delete(time, idx_del)
    data = np.delete(data, idx_del, axis=1)
    errs = np.delete(errs, idx_del, axis=1)
    return time, data, errs

def pool_fit_one(args):
    group = args['group']
    data = args['data']
    counter = args['counter']
    fitMod = args['fitter']
    names = args['names']
    errs = args['errs']
    time = args['time']
    i = args['idx']
    fitMod.logger = counter.add_fitInfo

    counter.set_curN(names)
    name_string = "{:10} {:25}".format(group, names)
    bestRes = fitMod.fit(data, errs, time, method=args['method'], name_string = name_string)
    parallelResults = (i, bestRes)
    print("{}: DONE".format(name_string))

    return counter, parallelResults, name_string

def main():
    parser = ArgumentParser()
    parser.add_argument("filename", type=str)
    parser.add_argument('-t', '--type', default='hdf', help="Type of using datafile. Can be: \'npy\', \'csv\', \'hdf\' (default)")
    parser.add_argument('-i', '--idata', type=str, help='index in data array for fit, e.g. 0-1')
    parser.add_argument('-s', '--exp-start', default=4, type=int, help='Number of exponents to start from')
    parser.add_argument('-f', '--exp-finish', default=6, type=int, help='Number of exponents when finish')
    parser.add_argument('-n', '--ntry', default=5, type=int, help='Number of tryings (for method NexpNtry)')
    parser.add_argument('-m', '--method', default='NexpNtry', type=str)
    parser.add_argument('-g', '--group', nargs='*', default=['NH'], help='Which group you want to fit. Need to fit data from hdf')
    parser.add_argument('--tcf', default='acf', help='Need to fit data from hdf')
    parser.add_argument('-o', '--output', default='out.hdf', help='filename for saving results')
    parser.add_argument('-c', '--time-cut', default=3.0, type=float,
                         help='time in ps which need to be cut from timeline')
    parser.add_argument('-1','--sum-one', action='store_true', help='Constrain the sum of all amplitudes to be exactly one (1.000) for ACFs')
    parser.add_argument('--nc',type=int, help='number of cpu (DEBUG)')
    # parser.add_argument('--lm', required=False, action='store_true', help='enable silent mode (DEBUG)')
    args = parser.parse_args()
    fid = h5py.File(args.output, 'w')
    counter = Counter()

    start=t.monotonic()
    for group in args.group:
        fitMod = Fitter(minexp=args.exp_start, maxexp=args.exp_finish,
                        ntry=args.ntry, tcf_type=args.tcf, sum_one=args.sum_one)
        counter.set_curMethod(args.method)
        counter.set_curTcf(args.tcf)
        counter.set_curGroup(group)

        try:
            time, data, errs, names = load_data(args, group)
            time, data, errs = prepare_data(time, data, errs, args.time_cut, args.sum_one)
        except Exception as e:
            print(e)
            continue

        data_size = data.shape[0]

    ## Prepare file for saving results
        grp = fid.create_group(group)
        tcf_grp = grp.create_group(args.tcf)
        exps = []
        for i in fitMod.exp_iter():
            cexp = tcf_grp.create_group('exp{}'.format(i))
            exps.append(cexp)
        for exp_grp, i in zip(exps, range(args.exp_start, args.exp_finish + 1)):
            nparams = 2 * i + 1
            exp_grp.create_dataset('params', data=np.zeros((data_size, nparams)))
            exp_grp.create_dataset('covar', data=np.zeros((data_size, nparams, nparams)))
            exp_grp.create_dataset('stats', data=u.create_nameArray(data_size, STAT_PARAMS_NAMES))

        if args.idata:
            s = args.idata.split('-')
            if len(s) == 1:
                fit_start = int(s[0])
                fit_end = fit_start + 1
            elif len(s) >= 2:
                fit_start = int(s[0])
                fit_end = int(s[1]) + 1
            else:
                fit_start = 0
                fit_end = data_size
            fit_start = 0 if fit_start < 0 else fit_start
            fit_end = data_size if fit_end > data_size else fit_end
        else:
            fit_start = 0
            fit_end = data_size

        # prepare arguments for parallel fitting
        arg_list = [{'data': data[i], 'errs': errs[i], 'time': time,
                     'idx': i, 'names': names[i], 'group': group,
                     'fitter': copy.copy(fitMod), 'method':args.method,
                     'counter': copy.deepcopy(counter)} for i in range(fit_start, fit_end)]
        # print(arg_list)
        # print("Start pool of {} CPU".format(nproc))
        with threadpool_limits(limits=1, user_api='blas'):
            pool = Pool()
            res_par = pool.map_async(pool_fit_one, arg_list)

            pool.close()
            pool.join()
            res_par.wait()
            result = res_par.get()

            try:
                for rc, rf, name_string in result:
                    counter = counter + rc
                    # print(type(rf), len(rf))
                    i, bestRes = rf
                    # print(type(bestRes), len(bestRes))
                    for group_hdf, res in zip(exps, bestRes):
                        if res and hasattr(res, 'success') and res.success:
                            # print(res)
                            group_hdf['params'][i] = res.param_vals
                            group_hdf['covar'][i]  = res.covar
                            ## !!! ОЧЕНЬ УПОРОТЫЙ МЕТОД ИЗ-ЗА НЕВОЗМОЖНОСТИ ПОЭЛЕМЕНТНОЙ ЗАМЕНЫ ЭЛЕМЕНТОВ ДАТАСЕТА.
                            stat_list = []
                            for key in STAT_PARAMS_NAMES:
                                vals = res.stats
                                stat_list += [vals[key]]

                            group_hdf['stats'][i] = tuple(stat_list)
                            ## КОНЕЦ УПОРОТОГО МОМЕНТА
                        else:
                            print("{}: fit failed".format(name_string), file=sys.stderr)

            except Exception as e:
                print("{}: ERROR in fit".format(name_string), file=sys.stderr)
                print(type(e), e, file=sys.stderr)

            print("{}: DONE".format(name_string))

    finish = t.monotonic()
    counter.set_overalltime(finish-start)
    counter.save('fitStatistic.csv')
    print(counter)
    # fitMod.save_toFile('out')

if __name__ == '__main__':
    main()
