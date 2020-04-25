import matplotlib.pyplot as plt
import numpy as np

from classes.exp_model import CModel

class Fitting():
    """docstring for Accf"""
    def __init__(self, data, std, timeline):
        self.data = data
        self.std = std
        self.time = timeline
        self.fit = None
        self.fmodel = None
        self.params = None

    def fitting(self, init_values=None):
        try:
            #################
            y = self.data
            length = y.shape[0]
            x = self.time
            #################
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

    def plot_fit(self):
        fig = plt.figure(1)
        # self.fit.model_fit.plot()
        times = self.time
        length = self.data.shape[0]
        x = (times + 1e-16) / 1000

        plt.errorbar(x[1:500:10], self.data[1:500:10], yerr=self.std[1:500:10], fmt='none')
        plt.plot(x[1:500], self.data[1:500], 'c-', label='init data')
        plt.plot(x[1:500], self.fit.model_fit.best_fit[1:500], 'r-', label='best fit')
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        ##################################################
        # cres = self.group.res
        # plt.savefig('./exppng/accf_{}_{}_{}_{}exp.png'.format(cres.resnr, cres.res, self.group.group_id, self.fmodel.nexp))
        ###################################################
        # plt.close()
        plt.show()

    def plot_logfit(self):
        fig = plt.figure(1)
        times = self.time
        length = self.data.shape[0]
        xb = (times + 1e-16) / 1000
        x = np.log(xb)

        # plt.xticks(times, x)
        # errticks = x[:3000] + x[3000::10000]
        plt.errorbar(x[1:500:10], self.data[1:500:10], yerr=self.std[1:500:10], fmt='none')
        plt.plot(x[1:500], self.data[1:500], 'c-', label='init data')
        plt.plot(x[1:500], self.fit.model_fit.best_fit[1:500], 'r-', label='best fit')
        plt.xlabel('t, ns')
        plt.ylabel('C(t)')
        plt.legend()
        print(self.data - self.fit.model_fit.best_fit)
        ##################################################
        # cres = self.group.res
        # plt.savefig('./exppng/accf_{}_{}_{}_{}exp.png'.format(cres.resnr, cres.res, self.group.group_id, self.fmodel.nexp))
        ###################################################
        # plt.close()
        plt.show()

    def get_result(self):
        return self.fit.report()
