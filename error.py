class Error(Exception):
    """Base class for exceptions in this module."""
    def __init__(self):
        self.lines = None


class ConfigFileError(Error):
    """Raised when a mistake in Config File occures
        Attributes:
        message -- explanation of why the specific transition is not allowed
    """

    def __init__(self, message):
        self.msg = message

    def __str__(self):
        result = 'Description:\n'
        for line_number, line in self.lines['Des']:
            result += str(line_number) + ' : ' + line + '\n'
        result += 'Group:\n'
        for line_number, line in self.lines['Group']:
            result += str(line_number) + ' : ' + line + '\n'
        return result


class DoubleElementError(ConfigFileError):
    """docstring for DoubleElementError"""

    def __init__(self, message, element):
        self.msg = message
        self.lines = element.lines


class EmptyElementError(ConfigFileError):
    """docstring for EmptyElementError"""

    def __init__(self, element):
        self.msg = 'Not full element.'
        self.lines = element.lines


class GroupElementError(ConfigFileError):
    """docstring for GroupElementError"""

    def __init__(self, message, element):
        self.msg = message
        self.lines = element.lines


class CoordError(ConfigFileError):
    """docstring for CoordError"""

    def __init__(self, element):
        self.msg = 'Missing coordinates in group'
        self.lines = element.lines


class TopologyError(Error):
    def __init__(self, msg):
        self.msg = msg


class NotFullError(Error):
    def __init__(self, elem):
        self.msg = 'Not full element'
        self.elem = elem


# class NotFullError(Error):
#       def __init__(self, mol=None, res=None, group=None, atom=None):
#               self.msg = 'Not full element'
#               self.mol = mol
#               self.res = res
#               self.group = group
#               self.atom = atom

class WrongAtomError(Error):
    def __init__(self, mol=None, res=None, group=None):
        self.msg = 'Wrong Atom'


class WrongResidueInPdb(Error):
    def __init__(self, mol, resno, res_topol, res_pdb):
        self.msg = 'Error: In molecute \'{}\' residue number ' + \
                   ' \'{}\' is \'{}\', but \'{}\' is expected.' + \
                   'Check topology/pdb sequence mismatch\n'
        self.msg = self.msg.format(mol.mol, resno, res_topol, res_pdb)


class MissingFrameError(Error):
    """docstring for MissingFrameError"""

    def __init__(self, nframe, model_nr):
        super(MissingFrameError, self).__init__()
        self.msg = 'Missing Frame'
        self.nframe = nframe
        self.model_nr = model_nr
