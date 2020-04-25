import time as t
import numpy as np
import matplotlib
matplotlib.use('Agg')
import multiprocessing.dummy as mp
import matplotlib.pyplot as plt

from numpy.fft import fft, ifft
# from lmfit import Model
from classes.exp_model import CModel
from utils import near_2pow


def exp_sum(x, aA, aT, bA, bT, c):
    return aA * np.exp(-aT * x) + bA*np.exp(-bT * x) + c


class Mdtrace:
    def __init__(self, md_count, mol=None, res=None, group=None, atom=None):
        self.mol    = mol
        self.res    = res
        self.group  = group
        self.atom   = atom
        self.coords = [None]*md_count
        self.md_count = md_count

        if mol != None:
            self.tframe = [None]*md_count

        if group != None:
            self.accfs   = [0]*md_count
            self.accf_mean = 0
            self.std = 0

            self.vec1 = None
            self.vec2 = None
            self.vec3 = None

            self.fit = None
            self.fmodel = None
            self.params = None

    def accf(self):
        if self.group == None:
            return
        self.group.add_vecs(self.md_count)

        for trace, md_num in zip(self.coords, range(self.md_count)):
            self.group.prepare_vec(md_num)

            try:
                self.accfs[md_num] = self.acf(md_num)
            except:
                self.accfs[md_num] = self.ccf(md_num)

    def acf(self, md_num):
        vec = self.vec1[md_num]
        length = vec.shape[0]
        new_length = near_2pow(2*length)
        vec = vec / (np.linalg.norm(vec, axis=-1, keepdims=True) + 1e-16)
        ct = 0
        for x in range(3):
            for y in range(x, 3):
                ct += self.n_aa((x, y), vec, new_length) * 1.5 * ((x != y) + 1)

        ct = ct[:length]
        coef = np.linspace(length, 1, length)
        acf = ct / coef - 0.5
        return acf

    def n_aa(self, axis, vec, n):
        length = vec.shape[0]

        vec_1 = vec[:, axis[0]]
        vec_2 = vec[:, axis[1]]
        FFT = fft(vec_1 * vec_2, n=n)
        return ifft(FFT.conj()*FFT, n=n)

    def calc_mean(self):
        add = lambda *args: sum(args)
        min_length = min(map(len, self.accfs))
        min_length = 2*10**5
        for i in range(self.md_count):
            self.accfs[i] = self.accfs[i][:min_length].real
            if np.isnan(np.sum(self.accfs[i])):
                print(self.group.group_id, self.group.res.res, self.group.res.resnr)
                print(self.accfs[i])
        try:
            self.accf_mean = np.array(list(map(add, *self.accfs))) / self.md_count
            self.std = np.sqrt(np.sum((self.accfs - self.accf_mean)**2, axis=0) / self.md_count) + 1e-16
        except Exception as e:
            print("Undefined exception while calc_mean: {} {}".format(type(e), e))
            raise

    def fitting(self, init_values=None):
        try:
            #if self.fmodel.nexp > 4:
            #    return
            y = self.accf_mean
            length = y.shape[0]
            tframe = self.tframe[0]
            x = np.arange(0, length*tframe, tframe)

            self.params = self.fmodel.prep_params()
            self.fit = self.fmodel.fit(y, self.params, x=x, method='least_squares', weights=1/self.std, nan_policy='omit', **init_values)
        except Exception as e:
            print("Undefined exception while fitting: {} {}".format(type(e), e))
            raise
        # try:
        #     if self.fmodel.check_errs():
        #         init_values = self.fit.model_fit.best_values
        #         self.fmodel.add_exp()
        #         self.fmodel.prep_params()
        #         self.fitting(init_values)
        # except Exception as e:
        #     print("Undefined exception while checking: {} {}".format(type(e), e))
        #     raise

    def plot_accf(self):
        tframe = self.tframe[0]
        length = self.accf_mean.shape[0]
        x = np.arange(0, length*tframe, tframe) / 1000
        cres = self.group.res
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.plot(x, self.accf_mean, label='ACF')
        plt.legend()
        plt.savefig('./acf/accf_{}_{}_{}.png'.format(cres.resnr, cres.res, self.group.group_id))
        plt.close()

    def plot_fit(self):
        fig = plt.figure(1)
        # self.fit.model_fit.plot()
        tframe = self.tframe[0]
        length = self.accf_mean.shape[0]
        times = np.arange(0, length*tframe, tframe)
        x = (times + 1e-16) / 1000
        # plt.xticks(times, x)
        # errticks = x[:3000] + x[3000::10000]
        plt.subplot(211)
        plt.errorbar(x[::10000], self.accf_mean[::10000], yerr=self.std[::10000], fmt='none')
        plt.plot(x, self.accf_mean, 'c-', label='init data')
        plt.plot(x, self.fit.model_fit.best_fit, 'r-', label='best fit')
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        plt.subplot(212)
        plt.errorbar(x[:6000], self.accf_mean[:6000], yerr=self.std[:6000], fmt='none')
        plt.plot(x[:6000], self.accf_mean[:6000], 'c-', label='init data')
        plt.plot(x[:6000], self.fit.model_fit.best_fit[:6000], 'r-', label='best fit')
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        cres = self.group.res
        plt.savefig('./exppng/accf_{}_{}_{}_{}exp.png'.format(cres.resnr, cres.res, self.group.group_id, self.fmodel.nexp))
        plt.close()
 
    def get_result(self):
        return self.fit.report()


if __name__ == '__main__':

    accf = Mdtrace(1)
    # print('Result: ', accf.acf(vec))
