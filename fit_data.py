#!/usr/bin/python3 -u
import sys

import numpy as np
import lmfit as lm

from classes.fitting import Fitting
from classes.exp_model import CModel

fname = sys.argv[1]

data = np.loadtxt(fname, delimiter=',')

time = data[:, 0]
func = data[:, 1]
errs = data[:, 2]
# ОЧЕНЬ ВАЖНО!!
errs[0] = errs[1]

exp_mod = Fitting(func[:], np.sqrt(errs[:]), time[:])
exp_mod.fmodel = CModel(5)

initial_values = {'aamplitude' : 0.01,
                  'bamplitude' : 0.01,
                  'camplitude' : 0.01,
                  'damplitude' : 0.01,
                  'eamplitude' : 0.01,
                  # 'famplitude' : 0.01,
                  'c' : 0.1,
                  'adecay' : 1,
                  'bdecay' : 10,
                  'cdecay' : 50,
                  'ddecay' : 100,
                  'edecay' : 1000}
                  # 'fdecay' : 3000}

# initial_values = {'aamplitude' : 0.01,
#                   'c' : 0.1,
#                   'adecay' : 1}

exp_mod.fitting(init_values=initial_values)
print(exp_mod.fit.model_fit.fit_report())
# d = exp_mod.fit.model_fit.best_values

print(exp_mod.fit.model_fit.covar)

def f(x, aamplitude, bamplitude, camplitude, damplitude, c, adecay, bdecay, cdecay, ddecay, eamplitude, edecay):
    return c + aamplitude * np.exp(-x/adecay) + bamplitude * np.exp(-x/bdecay) + camplitude * np.exp(-x/cdecay) +damplitude * np.exp(-x/ddecay) + eamplitude * np.exp(-x/edecay)
# print(f(0, **d))
# exp_mod.plot_logfit()