import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from classes.exp_model import CModel
from utils import near_2pow


class Accf():
    """docstring for Accf"""
    def __init__(self):
        self.accfs = None
        self.std = None
        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None

#########################################################

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

#########################################################

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