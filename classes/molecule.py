import time

import multiprocessing.dummy as mp
import numpy as np
from multiprocessing import cpu_count

NCPU = cpu_count()

from tqdm import tqdm
from error import NotFullError, Error
from classes.exp_model import CModel
from classes.relaxation_groups import Relaxation_Group_Types_Names

class Molecule:
    """docstring for molecule"""

    def __init__(self, molecule):
        self.mol = molecule
        self.seq = ''
        self.start = 0
        self.end = 0
        self.residues = []
        self.natoms = 0
        # self.lines = None

        self.trace = None
        # self.xyz = None

    def get_res(self, resnr, res=None):
        if resnr > self.end or resnr < self.start:
            raise Error
        res = self.residues[resnr - self.start]
        return res

    def res_in_mol(self, resnr, res=None):
        if resnr > self.end or resnr < self.start:
            return False
        elif res and self.residues[resnr - self.start].res != res:
            return False
        return True

    def check_mol(self, file):
        for elem in vars(self).values():
            if elem is None:
                print("ERROR!! In file {}".format(file))
            if type(elem) is list:
                for res in elem:
                    res.check_res(file)

    def put_coords(self, md_num):
        rline = 0
        gline = 0
        aline = 0
        for res in self.residues:
            for group in res.groups.values():
                group.trace.coords[md_num] = self.trace.coords[md_num][gline: aline, :, :]
                # group.xyz.xyz[md_num] = self.xyz.xyz[md_num][gline: aline, :, :]
                gline = aline
            res.trace.coords[md_num] = self.trace.coords[md_num][rline: aline, :, :]
            # res.xyz.xyz[md_num] = self.xyz.xyz[md_num][rline: aline, :, :]
            rline = aline

    # def prepare_coords(self, md_num):
    #     aline = 0
    #     for atom in self.atom_iter():
    #         atom.xyz.xyz[md_num] = self.xyz.xyz[md_num][aline]
    #         aline += 1

    def calculate_accf(self, groups):
        """"
        Calculate acf or ccf
        groups - groups which need to be calculated.(NH, NH2, CHAl, etc)
        """
        start = time.time()
        print('ACCF calculating start')

        length = len(groups)
        pool = mp.Pool()

        arg_list = [groups[i * length // NCPU: (i + 1) * length // NCPU]
                    for i in range(NCPU)]

        respar = pool.map_async(par_accf, arg_list)

        pool.close()
        pool.join()
        respar.wait()

        end = time.time()
        print('ACCF calculating finish {:.3} s'.format(end - start))

    def fit(self, groups, nexp):
        """"
        Fitting acf or ccf
        groups - groups which need to be fitted.(NH, NH2, CHAl, etc)
        """
        print('Fitting')
        start = time.time()

        pool = mp.Pool()
        length = len(groups)

        arg_list = [(groups[i * length // NCPU: (i + 1) * length // NCPU], nexp)
                    for i in range(NCPU)]

        respar = pool.map_async(par_fit, arg_list)

        pool.close()
        try:
            pool.join()
        except:
            raise
        respar.wait()

        end = time.time()
        print('Fitting finish {:.3} s'.format(end - start))

    def output(self, groups, fd):
        """"
        Calculate acf or ccf
        groups - groups which need to be calculated.(NH, NH2, CHAl, etc)
        """

        for group in groups:
            # print(group.trace.get_result())
            fd.write(str(group.res.resnr) + ' ' + group.group_id)
            fd.write('\n')
            fd.write(group.trace.get_result()[0])
            fd.write('\n')
            fd.write(str(group.trace.get_result()[1]))
            fd.write('\n')

    def get_group_types():
        types = set()
        for group in self.group_iter():
            types.add(group.simple_type)
        return list(types).sorted()

    def get_calc_groups(self, groups=Relaxation_Group_Types_Names):
        to_calc = []
        for group in self.group_iter(groups):
            to_calc.append(group)
        return to_calc

    def save_accfs(self, group_objects):
        fd = open('./accfs.npy', 'wb')
        for group in group_objects:
            np.save(fd, group.trace.accf_mean)
            np.save(fd, group.trace.std)
        fd.close()

    def save_accfs_to_files(self, group=Relaxation_Group_Types_Names):
        pass
        # for  in :
            # print("Saving group: ")


    def load_accfs(self, group_objects):
        fd = open('./accfs.npy', 'rb')
        for group in group_objects:
            group.trace.accf_mean = np.load(fd)
            group.trace.std = np.load(fd)
        fd.close()

    def save_fits(self, group_objects):
        fd = open('./fits.npy', 'wb')
        for group in group_objects:
            np.save(fd, np.array(group.trace.fit.report()))
        fd.close()

    def load_fits(self, group_objects):
        fd = open('./fits.npy', 'rb')
        for group in group_objects:
            group.trace.coefs = np.load(fd)
        fd.close()

    def dump(self, res=None, group=None, atom=None):
        result = "MOLECULE * * * * * * * * * * * * * * * \n"
        result += "Molecule Name: {}\n".format(self.mol)
        result += "Sequence:      {}\n".format(self.seq)
        result += "Start:         {}\n".format(self.start)
        result += "End:           {}\n".format(self.end)
        result += "Residues:\n"
        if res == 'all':
            for res in self.residues:
                result += res.dump(group, atom)
        elif res:
            result += res.dump(group, atom)
        else:
            result += '\tAre not None\n'
        result += " * * * * * * * * * * * * * * * * * * * \n"
        return result

    def res_iter(self):
        for res in self.residues:
            yield res

    def group_iter(self, group_type=Relaxation_Group_Types_Names):
        if type(group_type) is not list:
            group_type = list(group_type)
        for res in self.res_iter():
            for group in res.groups.values():
                if group.simple_type in group_type:
                    yield group

    def accf_iter(self, group_type=Relaxation_Group_Types_Names):
        for group in self.group_iter(group_type):
            for accf in group.accfs:
                yield accf

    def __lt__(self, other):
        if self.mol < other.mol:
            return self
        else:
            return other

    def __gt__(self, other):
        if self.mol > other.mol:
            return self
        else:
            return other


def par_accf(args):
    for group in args:
        group.trace.accf()
        group.trace.calc_mean()


def par_fit(args):
    groups, nexp = args
    try:
        for group in groups:
            group.trace.fmodel = CModel(nexp)
            group.trace.fitting()
            print('Huge success', group.trace.fit.model_fit.best_values)
    except Exception as e:
        print("Undefined exception while main: {} {}".format(type(e), e))
        raise
