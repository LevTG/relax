#! /usr/bin/python3 -u
import gc
import re
import os
import sys
import copy
import time
import subprocess

# не работает
# sys.path.insert(0, '/usr/local/gromacs/share/gromacs/top/')

import numpy as np
import _pickle as pickle

import error as err

from prettytable import PrettyTable
from tqdm import tqdm

from classes.molecule import Molecule
from classes.residue import Residue, get_letter_by_res
from classes.mdtrace import Mdtrace
from classes.config import Config
from classes.relaxation_groups import Relaxation_Group_Types_Names

from utils import split_comments, load_np

# line for includes in files
include_line_re = re.compile(r'^\#include "(?P<include_file>\S+)"$')

# line for comments
comment_line_re = re.compile(r'^\s*\;\s*(?P<comment>[\S\s]*)')

# line for name of group
sec_line_re = re.compile(r'^\[\s*(?P<title>[A-Za-z]+)\s+\d?\]$')

# line for moleule type
molecule_line_re = re.compile(r'^(?P<Name>[\S]+)\s+(?P<nrexcl>\d+)$')

# comment line for residue
res_line_re = re.compile(r'^residue\s+(?P<resnr>\d+)\s+(?P<residue>\w+)\s+')

# line pattern for atoms in topology
# for searching nr, type, residue, atom, cgnr
top_line_re = re.compile(r'(?P<nr>[0-9]+)\s+(?P<atom_type>[A-Z\d\*]+)\s+(?P<resnr>[0-9]+)\s+(?P<residue>\S{3,})\s+(?P<atom>[A-Z0-9]+)\s+(?P<cgnr>[0-9]+)')

##########################################################

class Topology:
    def __init__(self, logger=None):
        self.molecules = []
        self.logger = logger

        self.top_path = None
        self.config = None
        self.md_count = 0

    def mol_iter(self):
        for molecule in self.molecules:
            yield molecule

    def res_iter(self):
        for mol in self.mol_iter():
            for res in mol.res_iter():
                yield res

    def group_iter(self, group_type=Relaxation_Group_Types_Names):
        for res in self.res_iter():
            for group in res.groups.values():
                if group.simple_type in group_type:
                    yield group

    def atom_iter(self, group_type=Relaxation_Group_Types_Names):
        for group in self.group_iter(group_type):
            for atom in group.atoms_nn.values():
                yield atom

    def read_topology(self, read_config, file = './topol.top'):
        print("Reading topology starts...")
        if not read_config:
            print('Missing read config file. Please input config file and repeat.')
            return
        elif not os.path.exists(file):
            print('ERROR!! file {} not exist!'.format(file))
            return
        self.top_path = os.path.abspath(file)
        print(self.top_path)
        self.config = read_config

        topol = open(file, 'r')
        sec = ''
        cmol = None
        for self.line, self.line_number in split_comments(topol):
            line = self.line
            line_number = self.line_number
            try:
                if include_line_re.search(line):
                    file_name = include_line_re.search(line).group('include_file')
                    self.read_topology(read_config, './' + file_name)
                elif sec_line_re.search(line):
                    sec = sec_line_re.search(self.line).group('title')
                elif sec == 'moleculetype':
                    cmoln, cmol = self.read_molecule_line()
                elif sec == 'atoms':
                    self.read_atom_line(cmol)
                else:
                    break

            except err.TopologyError as e:
                print("ERROR: {}\nIn topology file: {}\n{} : {}\n".format(e.msg, self.top_path, self.line_number, self.line))
            except Exception as e:
                print("Undefined exception while reading: {} {}\n{} {}".format(type(e), e, line_number, line))
                raise e

        try:
            # print(self.molecules[0].dump(group='all', atom='all', res='all'))
            # self.check_topol(self.top_path)
            self.overall_table()

        except err.NotFullError as e:
            topol.close()
            print("ERROR: {}\n{}".format(e.msg, e.mol.dump(e.res, e.group, e.atom)))
        except Exception as e:
            print("Undefined exception while checking: ", e)
            raise e
        for atom in self.molecules[0].atom_iter():
            self.molecules[0].natoms += 1
        print("read_topology finish")

    def read_atom_line(self, cmol):
        line = self.line
        if top_line_re.search(line):
            found_group = top_line_re.search(line)
            results = found_group.groupdict()
            nr    = results['nr']
            resnr = int(results['resnr'])
            res   = results['residue'][:3]
            ann   = results['atom']

            cres = self.get_residue(resnr, res, cmol)
            relax_groups = self.config.get_group_by_res_atom(res, ann)
            if not relax_groups:
                return

            for relax_group in relax_groups:
                atom = Atom(ann)
                gtn = relax_group.group_id
                if gtn not in cres.groups.keys():
                    cres.groups[gtn] = copy.deepcopy(relax_group)
                    cres.groups[gtn].mol = cmol
                    cres.groups[gtn].res = cres

                cgroup = cres.groups[gtn]
                atn = cgroup.get_atn(ann)
                cgroup.set_atom(ann, atn, atom)
                cgroup.atoms_nn[ann].atom_updater(cmol, cres, cres.groups[gtn], nr)
        else:
            raise err.TopologyError('Wrong line in [ atoms ] section.')

    def get_residue(self, resnr, resn, cmol):
        if cmol.res_in_mol(resnr, resn):
            return cmol.get_res(resnr, resn)
        else:
            if cmol.start == 0:
                cmol.start = resnr
                cmol.end = resnr - 1

            cres = Residue(cmol, resn, resnr)
            cmol.residues.append(cres)
            cmol.seq += get_letter_by_res(resn)
            cmol.end += 1

            return cres

    def read_molecule_line(self):
        line = self.line

        if molecule_line_re.search(line):
            cmoln = molecule_line_re.search(line).group('Name')
            cmol = Molecule(cmoln)
            self.molecules.append(cmol)
        else:
            raise err.TopologyError('Wrong line in [ moleculetype ] section.')

        return cmoln, cmol

    def get_mol(self, molecule):
        return [el for el in self.molecules if el.mol == molecule][0]

    def check_topol(self, file):
        print("Checking topology starts...")
        for molecule in sorted(self.molecules):
            molecule.check_mol(file)
        print("Checking topology finish")

    def overall_table(self):
        table = ''
        tonum = lambda x: '' if x == 0 else x
        for molecule in self.molecules:
            # table += molecule.seq + '\n'
            x = PrettyTable()
            x.title = molecule.mol
            x.field_names = ['resnr', 'res', 'NH', 'NH2', 'CH Al', 'CH Ar', 'CH2', 'CH3', 'ERRORS']
            total = {'NH': 0, 'NH2': 0, 'CHAl': 0, 'CHAr': 0, 'CH2': 0, 'CH3': 0, 'ERRORS': 0}
            for res in molecule.residues:
                allg = {}
                for g in res.groups.values():
                    grt = g.residue_type
                    gst = g.simple_type
                    if grt not in allg:
                        allg[grt] = {'NH': 0, 'NH2': 0, 'CHAl': 0, 'CHAr': 0, 'CH2': 0, 'CH3': 0}
                    allg[grt][gst] += 1
                    total[gst] += 1
                for g in res.groups.values():
                    grt = g.residue_type
                    for gst in allg[grt]:
                        allg[grt][gst] = tonum(allg[grt][gst])
                resnr = res.resnr
                resn = res.res
                x.add_row([str(resnr), resn, allg[resn]['NH'], allg[resn]['NH2'], allg[resn]['CHAl'], allg[resn]['CHAr'], allg[resn]['CH2'], allg[resn]['CH3'], tonum(res.errors)])
                total['ERRORS'] += res.errors
            x.add_row(['', 'TOTAL', total['NH'], total['NH2'], total['CHAl'], total['CHAr'], total['CH2'], total['CH3'], total['ERRORS']])
            table += str(x) + '\n'
        print(table)

    def read_coords(self):
        pass

    def save(self, file):
        fd = open(file, 'wb')
        pickle.dump(self, fd, protocol=4)
        fd.close()

    def save_coords(self, dump, file):
        # header = np.array([[a.res.res, a.res.resnr, a.group.group_id, a.name, a.nr] for a in self.molecules[0]])
        fd = open(file, 'wb')
        np.save(fd, dump)
        tframes = np.array(self.molecules[0].trace.tframe)
        np.save(fd, tframes)
        for trace in self.molecules[0].trace.coords:
            np.save(fd, trace)
        fd.close()

    def load(self, file):
        fd = open(file, 'rb')

        tmp = np.load(fd)
        obj = pickle.loads(tmp)
        cmol = obj.molecules[0]
        tframes = np.load(fd)
        obj.md_count = len(tframes)
        obj.add_traces(obj.md_count)
        cmol.trace.tframe[:] = tframes
        for i, trace in zip(range(obj.md_count), load_np(obj.md_count, fd)):
            cmol.trace.coords[i] = trace
        for i in range(obj.md_count):
            cmol.put_coords(i)
        self.check_all_coords()
        fd.close()
        return obj


def main():
    test_cnf = Config()
    test_cnf.read_config()
    test = Topology()
    test.read_topology(test_cnf, sys.argv[1])
    # print(test.molecules[0].residues[0])


if __name__ == '__main__':
    main()
