import re
import copy

from error import DoubleElementError, NotFullError, EmptyElementError

import utils as u
import numpy as np


class group_type:
    """
    Abstract class for relaxation groups\n
    group_type_name     -   (AminesSecondary, etc)
    group_id                group users name (in config)
    molecule            -   molecule this group refer to
    residue_number      -   residue number this group refer to
    residue_type        -   residue type this group refer to (LIS, ALA, etc.)

    """

    def __init__(self):
        self.group_type_name = None
        self.group_id = None
        self.simple_type = None
        self.residue_type = None
        self.atom_type_names = []
        self.atom_native_names = []
        self.coord = {}
        self.lines = {'Des': [], 'Group': []}

        #####################################
        self.mol = None
        self.res = None
        self.trace = None
        # self.xyz = None
        self.atoms_nn = {}
        self.atoms_tn = {}
        self.acfs = []
        self.acfs_data =None
        self.ccfs = []
        self.ccfs_data = None
        #####################################

@property
    def natoms(self):
        return len(self.atom_type_names)

    def have_atom_type(self, atom_type):
        return atom_type in self.atom_type_names

    def is_atom_in_group(self, ann=None, atn=None):
        return ann in self.atoms_nn.keys() or atn in self.atoms_tn.keys()

    def set_empty_atom_names(self):
        self.atom_native_names = [None for _ in self.atom_type_names]

    def set_residue(self, res):
        self.residue_type = res

    def set_lines(self, line, line_number, type_):
        self.lines[type_].append((line_number, line))

    def set_atom_names(self, atn, ann):
        self.atom_type_names.append(atn)
        self.atom_native_names.append(ann)

    def set_initial_coordinates(self, atom_type_name, coordinates):
        self.coord[atom_type_name] = coordinates[:]

    def set_coordinates(self, coordinates, atom_tn=None, atom_nn=None):
        if atom_tn:
            self.atoms_tn[atom_tn].coords = coordinates
        if atom_nn:
            self.atoms_nn[atom_nn].coords = coordinates

    def get_atn(self, ann):
        idx = np.inf
        for i in range(len(self.atom_type_names)):
            if ann in self.atom_native_names[i]:
                idx = i
        return self.atom_type_names[idx]

    def get_ann(self, atn):
        idx = np.inf
        for i in range(len(self.atom_type_names)):
            if atn in self.atom_type_names[i]:
                idx = i
        return self.atom_native_names[idx]

    # def accf(self):
    #     self.prepare_vec()
    #     self.trace.acf()
    #     # except:
    #     #     self.trace.ccf()

    def dump(self):
        result = "= = = = = = = = = = = = = = = = = = = =\n"
        result += "Group Type Name: {}\n".format(self.group_type_name)
        # result += "Group Id:        {}\n".format(self.group_id)
        result += "Simple Type:     {}\n".format(self.simple_type)
        result += "Residue:         {}\n".format(self.residue_type)
        result += "Coordinates:     {}\n".format(self.coord)
        if 'atoms_nn' not in vars(self).keys():
            result += "= = = = = = = = = = = = = = = = = = = =\n"
            return result
        except Exception as e:
            raise e
        result += "= = = = = = = = = = = = = = = = = = = =\n"
        return result

    def check_elem(self, file=None):
        e = EmptyElementError(self)
        for key, elem in vars(self).items():
            if elem is None:
                raise e
            if type(elem) is list:
                if u.not_full_list(elem):
                    raise e
            elif type(elem) is dict and key == 'coord':
                if u.not_full_dict(elem):
                    raise e

    def full_check_elem(self):
        if 'atoms' in vars(self).keys():
            e = NotFullError(elem=self)
        else:
            e = EmptyElementError(self)
        for key, elem in vars(self).items():
            if elem is None:
                raise e
            if type(elem) is list:
                if u.not_full_list(elem):
                    raise e
            elif type(elem) is dict and key == 'coord':
                if u.not_full_dict(elem):
                    raise e


class group_AminesSecondary(group_type):
    def __init__(self):
        super(group_AminesSecondary, self).__init__()
        self.group_type_name = "AminesSecondary"
        self.simple_type = "NH"
        # self.atom_type_names = ["N", "H"]
        self.set_empty_atom_names()


class group_AminesPrimary(group_type):
    """docstring for group_AminesPrimary"""

    def __init__(self):
        super(group_AminesPrimary, self).__init__()
        self.group_type_name = "AminesPrimary"
        self.simple_type = "NH2"
        # self.atom_type_names = ["N", "H1", "H2"]
        self.set_empty_atom_names()


class group_MethinesAlif(group_type):
    """docstring for group_Methines"""

    def __init__(self):
        super(group_MethinesAlif, self).__init__()
        self.group_type_name = "MethinesAlif"
        self.simple_type = "CHAl"
        # self.atom_type_names = ["C", "H"]
        self.set_empty_atom_names()


class group_MethinesArom(group_type):
    """docstring for group_Methines"""

    def __init__(self):
        super(group_MethinesArom, self).__init__()
        self.group_type_name = "MethinesArom"
        self.simple_type = "CHAr"
        # self.atom_type_names = ["C", "H"]
        self.set_empty_atom_names()


class group_Methylenes(group_type):
    """docstring for group_Methylenes"""

    def __init__(self):
        super(group_Methylenes, self).__init__()
        self.group_type_name = "Methylenes"
        self.simple_type = "CH2"
        # self.atom_type_names = ["C", "H1", "H2"]
        self.set_empty_atom_names()


class group_Methyles(group_type):
    """docstring for group_Methyles"""

    def __init__(self):
        super(group_Methyles, self).__init__()
        self.group_type_name = "Methyles"
        self.simple_type = "CH3"
        # self.atom_type_names = ["C", "H1", "H2", "H3"]
        self.set_empty_atom_names()


Relaxation_Group_Types = {
    'NH': group_AminesSecondary(),
    'NH2': group_AminesPrimary(),
    'CHAl': group_MethinesAlif(),
    'CHAr': group_MethinesArom(),
    'CH2': group_Methylenes(),
    'CH3': group_Methyles()
}

Relaxation_Group_Types_Names = list(Relaxation_Group_Types.keys())

def des_group_constructor(simple_type_name, line, line_number):
    new_group = copy.deepcopy(Relaxation_Group_Types[simple_type_name])
    new_group.set_lines(line, line_number, 'Des')
    return new_group


def group_constructor(elem, simple_type_name, line, line_number, residue, group_id):
    new_group = copy.deepcopy(elem.groups_des[simple_type_name])
    new_group.set_residue(residue)
    new_group.group_id = group_id
    new_group.set_lines(line, line_number, 'Group')
    return new_group


def group_updater(group, atom_type_name, line, line_number, type_, atom_native_name=None, coord=None):
    group.set_lines(line, line_number, type_)

    if coord:
        group.set_initial_coordinates(atom_type_name, coord)
