import sys
import h5py

import numpy as np
import copy as c
import matplotlib.pyplot as plt

from fitter.exp_model import CModel

BASE_VAL = np.array([0.1, 1, 0.01, 10, 0.01, 50, 0.01, 100, 0.01, 1000, 0.01, 3000, 0.01])
BASE_KEY = ['c', 'adecay', 'aamplitude', 'bdecay', 'bamplitude', 'cdecay', 'camplitude', 'ddecay', 'damplitude', 'edecay', 'eamplitude', 'fdecay', 'famplitude']

rndmizer = np.random.RandomState()

class Fitter:
    def __init__(self, minexp=4, maxexp=6, randFactor=(.2, 5), ntry=5):
        self.model = None
        self.ntry  = ntry
        self.randFactor  = randFactor
        self.expInterval = (minexp, maxexp)

        self.data = None
        self.std  = None
        self.time = None

        self.nexp = 0
        self.params  = None
        self.lastFit = None
        self.init_values = None

        self.res = None
        self.succes = False

        self.fileToSave = h5py.File('fit_result.hdf5', 'w')

    ## Setters

    def set_nexp(self, minexp, maxexp):
        self.expInterval = (minexp, maxexp)
        self.prep_model()

    def set_time(self, timeline=None, dt=None, nframe=None):
        if timeline:
            self.time = timeline
        elif dt and nframe:
            self.time = np.array([i*dt for i in range(nframe)])
        else:
            print('ERROR!! Missing arguments. Please set timeline or dt with number of frames', file=sys.stderr)

    def set_ntry(self, ntry):
        self.ntry = ntry

    def set_randomFactor(self, randFactor):
        self.randFactor = randFactor


    ## Prepare to Fit

    def prep_model(self):
        self.nexp = self.expInterval[0]
        self.model = CModel(self.expInterval[0])
        self.init_values = dict(zip(BASE_KEY[:(2*self.nexp + 1)], c.copy(BASE_VAL[:(2*self.nexp + 1)])))




    def fit(self, mode='NexpNtry', gname, tcf, i, *args, **kwargs):
        self.clear()
        self.gname = gname
        self.tcf = tcf
        self.i = i
        if mode == 'NexpNtry':
            self.fit_NtryNexp(*args, **kwargs)
        else:
            print('ERROR!! You choose wrong mode', file=sys.stderr)

    def _fit(self, init_values=None):
        #if self.model.nexp > 4:
        #    return
        #################
        y = self.data
        length = y.shape[0]
        x = self.time
        #################
        self.lastFit = self.model.fit(y, x=x, method='least_squares', weights=1/self.std, nan_policy='omit', **init_values)

    def fit_NtryNexp(self, data, std, timeline, init_values=None):
        self.data = data
        self.std = std
        self.time = timeline
        self.prep_model()

        for self.nexp in range(self.expInterval[0], self.expInterval[1] + 1):
            self.fit_ntry()
            try:
                if not self.res:
                    self.res = self.lastFit

                elif self.res.chisqr > self.lastFit.chisqr:
                    self.res = self.lastFit

                if self.nexp == self.expInterval[1]:
                    print('res:', self.res.chisqr)
                    return self.res
            except AttributeError as e:
                print("ERROR!! res is None. It happend while fitting {} exponents. With those initial values: {}".format(self.nexp, self.init_values), file=sys.stderr)
                self.succes = False
                return

            self.model.add_exp()
            self.init_values = dict(zip(BASE_KEY[:(2*self.nexp + 1)], c.copy(BASE_VAL[:(2*self.nexp + 1)])))
        return self.res

    def fit_ntry(self):
        curFit = None

        for _ in range(self.ntry):
            self._fit(self.init_values)
            if not self.model.has_covar():
                pass
            elif not curFit:
                curFit = self.lastFit
            elif curFit and curFit.chisqr > self.lastFit.chisqr:
                curFit = self.lastFit

            if curFit:
                self.succes = True
                print(self.lastFit.chisqr, curFit.chisqr)
            self.change_init()

        self.lastFit = curFit

    def change_init(self):
        rndmizer.seed()
        new_val = c.copy(BASE_VAL)
        new_val[1::2] = new_val[1::2] * rndmizer.uniform(*self.randFactor)
        self.init_values = dict(zip(BASE_KEY[:(2*self.nexp + 1)], new_val[:(2*self.nexp + 1)]))

    def plot_fit(self, name):
        fig = plt.figure(1)
        # self.res.model_fit.plot()
        # step = self.time
        times=self.time
        length = self.data.shape[0]
        # times = np.arange(0, length*step, step)
        x = (times + 1e-16) / 1000
        # plt.xticks(times, x)
        # errticks = x[:3000] + x[3000::10000]
        plt.subplot(211)
        ##################################################
        plt.errorbar(x[::10000], self.data[::10000], yerr=self.std[::10000], fmt='none')
        ##################################################
        plt.plot(x, self.data, 'c-', label='init data')
        plt.plot(x, self.res.best_fit, 'r-', label='best fit')
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        plt.subplot(212)
        ##################################################
        plt.errorbar(x[:1000:100], self.data[:1000:100], yerr=self.std[:1000:100], fmt='none')
        plt.plot(x[:1000:100], self.data[:1000:100], 'c-', label='init data')
        plt.plot(x[:1000:100], self.res.best_fit[:1000:100], 'r-', label='best fit')
        ##################################################
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        ##################################################
        # cres = self.group.res
        plt.savefig('exp{}.png'.format(name))
        ###################################################
        plt.close()

    def get_result(self):
        return self.res.report()

    def clear(self):
        self.data = None
        self.std = None
        self.time = None

        self.nexp = 0
        self.res = None
        self.params = None
        self.init_values = None
        self.lastFit = None

        self.succes = False