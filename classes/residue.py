from error import NotFullError, WrongAtomError


class Residue:
    """docstring for residue"""
    def __init__(self, mol, res, res_num):
        self.mol = mol
        self.res = res
        self.resnr = res_num
        self.groups = {}
        self.errors = 0

    def group_iter(self):
        for group_name, group in self.groups.items():
            yield group_name, group

    def check_res(self, file):
        dels = []
        for elem in vars(self).values():
            if elem is None:
                return
        for group_name, group in self.group_iter():
            try:
                group.check_elem(file)
            except NotFullError as e:
                print('Removed {}'.format(e.elem.group_id))
                dels.append(group_name)
                self.errors += 1
        for del_name in dels:
            del self.groups[del_name]


    # def full_check_res(self):
    #     for elem in vars(self).values():
    #         if elem is None:
    #             raise NotFullError(res=self)
    #         if type(elem) is dict:
    #             for group in elem.values():
    #                 try:
    #                     if group:
    #                         group.check_elem(f)
    #                     else:
    #                         raise NotFullError(res=self)
    #                 except NotFullError as e:
    #                     e.res = self
    #                     raise e

    def __str__(self):
        result = "RESIDUE- - - - - - - - - - - - - - - - \n"
        result += "Residue name:   {}\n".format(self.res)
        result += "Molecule:       {}\n".format(self.mol.mol)
        result += "Residue number: {}\n".format(self.resnr)
        result += "Groups:\n"
        for group in self.groups.values():
            result += str(group)
        result += " - - - - - - - - - - - - - - - - - - - \n"
        return result

    def dump(self, group=None):
        result = "RESIDUE- - - - - - - - - - - - - - - - \n"
        result += "Residue name:   {}\n".format(self.res)
        result += "Molecule:       {}\n".format(self.mol.mol)
        result += "Residue number: {}\n".format(self.resnr)
        result += "Groups:\n"
        if group == 'all':
            for group in self.groups.values():
                result += group.dump()
        elif group:
            result += group.dump()
        else:
            result += '\tAre not None\n'
        result += " - - - - - - - - - - - - - - - - - - - \n"
        return result


letter_res = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'PYL': 'O',
    'SEC': 'U',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
}


def get_letter_by_res(residue):
    return letter_res[residue]
