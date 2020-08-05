import os
import sys
import h5py

import numpy as np
import utils as u

REFERENCE_PARAMS_SEQ = ['c', 'aamplitude', 'adecay', 'bamplitude', 'bdecay', 'camplitude', 'cdecay',
            'damplitude', 'ddecay', 'eamplitude', 'edecay', 'famplitude', 'fdecay',
            'gamplitude', 'gdecay','hamplitude', 'hdecay','iamplitude', 'idecay']

class FitResult():
    def __init__(self, init_values=None, params=None, covar_vars=None,
                 covar=None, stats=None, nexp=None, success=None):
        self._nexp = nexp
        self._init_values = init_values

        self._param_names = list(params.keys()) if params else list(init_values.keys())
        self._param_vals = self.sort_vals(params) if params else None
        # self._param_errs = None

        self._covar = covar
        self._covar_names = covar_vars
        self._stats = stats
        self.success = success

# В версии python 3.8 нужно через post_init для гарантии выполнения после всех присвоений
        self._param_errs = None if not success else np.sqrt(np.diag(self._covar))
        if success:
            self.sort_covar()

    # def __post_init__(self):
    #     self._param_errs = np.sqrt(np.diag(self._covar))
    #     self.sort_covar()


    def calc_paramErrs(self):
        self._param_errs = np.sqrt(np.diag(self._covar))

    def sort_covar(self):
        #print("sort_covar started")
        idx = []
        for ref_par, _ in zip(REFERENCE_PARAMS_SEQ, self._covar_names):
            idx.append(self._covar_names.index(ref_par))
        #print("sort_covar: _covar_names = {}".format(self._covar_names))
        #print("sort_covar: idx = {}".format(idx))
        new_covar = self.covar[idx]
        new_covar = new_covar[:, idx]
        self.covar = new_covar

    def sort_vals(self, params):
        par_vals = np.zeros(len(params))
        for i, par_key in enumerate(REFERENCE_PARAMS_SEQ):
            if par_key in params:
                par_vals[i] = params[par_key]
        return par_vals


    # Getters
    @property
    def nexp(self):
        return self._nexp

    @nexp.setter
    def nexp(self, nexp):
        self._nexp = nexp

    @property
    def initialValues(self):
        return self._init_values

    @initialValues.setter
    def initialValues(self, init_values):
        self._init_values = init_values

    @property
    def param_vals(self):
        return self._param_vals

    @param_vals.setter
    def param_vals(self, param_vals):
        self._param_vals = param_vals

    @property
    def param_names(self):
        return self._param_names


    @param_names.setter
    def param_names(self, param_names):
        self._param_names = param_names

    @property
    def covar(self):
        return self._covar

    @covar.setter
    def covar(self, covar_matrix):
        self._covar = covar_matrix

    @property
    def paramErrors(self):
        if not self._param_errs.all():
            self.calc_paramErrs()
        return self._param_errs

    @property
    def stats(self, stat=None):
        if not stat:
            return self._stats
        elif stat in self._stats.keys():
            return self._stats[stat]
        else:
            print('ERROR!! there no {} in fit statistic'.format(stat), file=sys.stderr)

    def __repr__(self):
        output = "covar: {}\nstats: {}\nparam vals: {}\nparam names{}".format(self._covar, self._stats, self._param_vals, self._param_names)
        return "Fit Result: " + output

    def __str__(self):
        output = "Fit results\n"
        output += 'Number of exponents: {}\n'.format(self.nexp)
        # output += 'Используемый метод: {}\n'.format(self.model)
        output += 'Parameters: {}\n'.format(self._param_names)
        output += 'Initial values of the parameters: {}\n'.format(self.initialValues)

        if self.success:
            output += 'Successfull fitting\n'
            output += 'Parameters are found: {}\n'.format(self._param_vals)
            output += 'Covariance matrix is found:\n{}\n'.format(self.covar)
            output += 'Order of parameters in covariance matrix:\n{}\n'.format(self._covar_names)
            output += 'Statistics: {}\n'.format(self.stats)
        else:
            output += 'Fitting failed\n'
        return output
