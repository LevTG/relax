import sys

import numpy as np

from lmfit import Model
from lmfit.models import ConstantModel, ExponentialModel


def exp(x, A, T):
    return A*np.exp(-T*x)


class CModel(Model):
    """docstring for CModel"""
    def __init__(self, nexp=0, sum_one=False):
        self.sum_one_flag = sum_one
        self.model = ConstantModel()
        self.nexp = 0
        self.errors = None
        self.params = None
        self.res = None
        for _ in range(nexp):
            self.add_exp()

    def add_exp(self):
        self.nexp += 1
        self.model += ExponentialModel(prefix=Prefixes[self.nexp])


    def prep_params(self, tcf_type):
        if tcf_type == 'acf':
            self.prep_acfParams()
        elif tcf_type == 'ccf':
            self.prep_ccfParams()

    def prep_ccfParams(self):
        amp = '{}amplitude'
        dec = '{}decay'
        self.model.set_param_hint('c', value=1.0)
        for i in range(1, self.nexp + 1):
            currA = amp.format(Prefixes[i])
            currT = dec.format(Prefixes[i])
            self.model.set_param_hint(currA, value=1, min=-1.5, max=1.5)
            self.model.set_param_hint(currT, value=1, min=0)

        self.params = self.model.make_params()

    def prep_acfParams(self):
        amp = '{}amplitude'
        dec = '{}decay'
        cntrl_expr = ' '
        if self.sum_one_flag:
            cntrl_expr = 'c'
        self.model.set_param_hint('c', value=1.0, min=0)
        for i in range(1, self.nexp + 1):
            currA = amp.format(Prefixes[i])
            currT = dec.format(Prefixes[i])
            self.model.set_param_hint(currA, value=1.0, min=0.0)
            self.model.set_param_hint(currT, value=1.0, min=0.0)
            if self.sum_one_flag:
                cntrl_expr += ' + ' + currA

        self.params = self.model.make_params()
        if self.sum_one_flag:
            self.params.add('cntrl', value= 1, max= 1+1e-5, min= 1-1e-5)
            self.params['cntrl'].vary = True
            self.params['cntrl'].expr = cntrl_expr


    def fit(self, tcf_type='acf', *args, **kwargs):
        try:
            self.prep_params(tcf_type)
            self.res = self.model.fit(params=self.params, *args, **kwargs)
        except Exception as e:
            print("Undefined exception while fit: {} {}".format(type(e), e), file=sys.stderr)
            raise
        return self.res

    def has_covar(self):
        covar = self.res.covar
        if covar is None or np.isnan(np.sum(covar)):
            return False
        else:
            return True

    def check_errors(self):
        for ivar, name in enumerate(self.res.var_names):
            par = self.res.params[name]
            if par.stderr / par.value > 1:
                return False
        return True

    def report(self):
        return self.res.fit_report(), self.res.covar


Prefixes = {
    1: 'a',
    2: 'b',
    3: 'c',
    4: 'd',
    5: 'e',
    6: 'f',
    7: 'g',
    8: 'h',
    9: 'i',
    10: 'j'
}
