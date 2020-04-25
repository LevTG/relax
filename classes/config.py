import re
import os.path
import logging

from prettytable import PrettyTable
# from collections import namedtuple

from utils import split_comments

import classes.relaxation_groups as rg
import error as err
# line for name of group
name_line_re = re.compile(r'^\[\s*(?P<title>[A-Za-z]+)\s*\]$')


des_title_line_re = re.compile(r'^\[\s*(?P<title>[A-Za-z]+)\s+(?P<simple_type>[\w\d]+)\s*\]$')


des_group_line_re = re.compile(r'^ATOM\s+(?P<atom_type_name>[A-Z]+\d*)\s+(?P<coordinates>(?:[-+]?\d+\.*\d*\s*){3})?$')

# line for config file
# residue_group_type (Gly_NH, etc)
group_info_line_re = re.compile(r'^\[\s*(?P<group_type_name>[\w\d]+)\s+\:\s+(?P<simple_type>[\w\d]+)\s+(?P<residue>[A-Z]{3,})\s*\]$')

# line for group in config file
group_line_re = re.compile(r'^ATOM\s+(?P<residue>[A-Z]{3,})\s+(?P<atom_type_name>[A-Z]+\d*)\s+(?P<atom_native_name>(?:[\w]+\d?\,?){1,})\s*(?P<coordinates>(?:[-+]?\d+\.*\d*\s*){3})?$')

# line for acf in config file
acf_line_re = re.compile(r'^ACF\s+(?P<atom1>[\w\d]+)-(?P<atom2>[\w\d]+)$')

logger = logging.getLogger()


class Config:
    def __init__(self, logger=None):
        self.group = {}
        self.logger = logger
        self.path = None
        self.groups_des = {}

    def read_config(self, file = './src/md2nmr.cfg'):
        print("Read config starts...")
        if not os.path.exists(file):
            print("ERROR! File {} not exist!".format(file))
            print('Read config finish')
            return

        self.path = os.path.abspath(file)
        config_file = open(file, 'r')
        part = ''
        has_coord = None
        read_group = ''
        stn = ''

        for self.line, self.line_number in split_comments(config_file):
            try:
                if name_line_re.search(self.line):
                    part = name_line_re.search(self.line).group('title')
                elif part == 'Descriptions':
                    stn = self.read_des_line(self.groups_des, stn)
                elif part == 'Groups':
                    read_group, has_coord = self.read_group_line(read_group, has_coord)
                elif part == 'Formulas':
                    pass

            except err.DoubleElementError as e:
                print("WARNING: {}\nIn config file: {}\n{}".format(e.msg, self.path, e))
            except err.GroupElementError as e:
                print("ERROR: {}\nIn config file: {}\n{}".format(e.msg, self.path, e))
                raise
            except KeyError as e:
                print("ERROR: Absense description for this group.\n{} : {}".format(self.line_number, self.line))
                raise
            except err.CoordError as e:
                print("ERROR: {}\nIn config file: {}\n{}".format(e.msg, self.path, e))
                raise
            except err.ConfigFileError as e:
                print("ERROR: {}\nIn config file: {}\n{} : {}\n".format(e.msg, self.path, self.line_number, self.line))
                raise
            except Exception as e:
                print("ERROR: Undefined exception: {} {}\nIn config file: {}\n{} : {}\n".format(type(e), e, self.path, self.line_number, self.line))
                raise
        # try:
            # self.check_config()
        # except err.EmptyElementError as e:
            # print("ERROR: {}\nIn config file: {}\n{}".format(e.msg, self.path, e))
            # logger.exception()
            # return
        config_file.close()
        self.overall_table()
        print("Read config finish")

    def read_des_line(self, groups_des, stn):
        line = self.line
        line_number = self.line_number

        if des_title_line_re.search(line):
            results = des_title_line_re.search(line).groupdict()
            stn = results['simple_type']

            if stn not in groups_des:
                groups_des[stn] = rg.des_group_constructor(stn, line, line_number)
            return stn

        elif des_group_line_re.search(line):
            results = des_group_line_re.search(line).groupdict()
            coord = list(map(int, results['coordinates'].split()))
            atn = results['atom_type_name']

            if atn[0] not in stn:
                raise err.GroupElementError('Wrong atom in definitions.', groups_des[stn])

            rg.group_updater(groups_des[stn], atn, line, line_number, 'Des', coord=coord)
            return stn

        else :
            raise err.ConfigFileError('Wrong line.')

    def read_group_line(self, read_group, has_coord):
        line = self.line
        line_number = self.line_number

        if name_line_re.search(line):
            return None, has_coord
        elif group_info_line_re.search(line):
            read_group = self.get_group_info()
            return read_group, None
        else:
            found_group = group_line_re.search(line)
            if found_group == None:
                raise err.ConfigFileError('Wrong line.')

            group_dict = found_group.groupdict()
            res = group_dict['residue']
            atn = group_dict['atom_type_name']
            ann = group_dict['atom_native_name']
            coord = group_dict['coordinates']

            ann = self.get_ann(ann)
            coord, has_coord = self.get_coord(coord, has_coord, read_group)

            if read_group.residue_type != res:
                raise err.GroupElementError('Wrong residue in line.', read_group)
            # elif atn[0] not in read_group.simple_type:
            #     read_group.add_lines(line, line_number, 'Group')
            #     raise GroupElementError('Wrong atom in line.', read_group)

            rg.group_updater(read_group, atn, line, line_number, 'Group', ann, coord)
            return read_group, has_coord

    def get_group_info(self):
        line = self.line
        line_number = self.line_number

        found_group = group_info_line_re.search(line)
        group_dict = found_group.groupdict()
        gtn = group_dict['group_type_name']
        res = group_dict['residue']
        stn = group_dict['simple_type']

        if stn not in rg.Relaxation_Group_Types.keys():
            raise err.ConfigFileError('Wrong group name.')

        if gtn not in self.group:
            self.group[gtn] = rg.group_constructor(self, stn, line, line_number, res, gtn)

        return self.group[gtn]

    def get_coord(self, coord, has_coord, read_group):
        if coord:
            coord = list(map(float, coord.split()))
            if has_coord == None:
                has_coord = True
            elif has_coord != True:
                read_group.set_lines(self.line, self.line_number, 'Group')
                raise err.CoordError(read_group)
        elif has_coord == None:
            has_coord = False
            # return (None, has_coord)
        elif has_coord != None and has_coord != False:
            read_group.set_lines(self.line, self.line_number, 'Group')
            raise err.CoordError(read_group)
        return (coord, has_coord)

    def get_ann(self, ann):
        ann = tuple(ann.split(','))
        return ann

    # FICHA_1
    def write_config(self, file='./config_new'):
        print("------------------------------------------------------------------------------")
        print("write_config starts...")
        out_file = open(file, 'w')
        line_number = 0
        gtn = ''

        for elem in self.group.values():
            # перебор по группам (AminesSecondary, etc)
            print(elem[gtn],file=out_file)
            for el in elem.values():
                print(el, file=out_file)
        print("write_config finished")
        print("------------------------------------------------------------------------------")
        # DOIT

    def check_config(self):
        print("Check config starts...")
        for group_type_name, group in sorted(self.group.items()):
            try:
                group.check_elem()
            except err.EmptyElementError as e:
                raise e
        print("Check config finish")

    def overall_table(self):
        allg = {}
        total = {'NH': 0, 'NH2': 0, 'CHAl': 0, 'CHAr': 0, 'CH2': 0, 'CH3': 0}

        for g in self.group.values():
            grt = g.residue_type
            gst = g.simple_type
            if grt not in allg:
                allg[grt] = {'NH': 0, 'NH2': 0, 'CHAl': 0, 'CHAr': 0, 'CH2': 0, 'CH3': 0}
            allg[grt][gst] += 1
            total[gst] += 1

        for g in self.group.values():
            grt = g.residue_type
            for gst in allg[grt]:
                if allg[grt][gst] == 0:
                    allg[grt][gst] = ''

        x = PrettyTable()
        x.field_names = ['', 'NH', 'NH2', 'CH Al', 'CH Ar', 'CH2', 'CH3']

        for res in sorted(allg.keys()):
            x.add_row([res, allg[res]['NH'], allg[res]['NH2'], allg[res]['CHAl'], allg[res]['CHAr'], allg[res]['CH2'], allg[res]['CH3']])

        x.add_row(['TOTAL', total['NH'], total['NH2'], total['CHAl'], total['CHAr'], total['CH2'], total['CH3']])

        print(x)

    def get_group_by_type(self, simple_type):
        return [el for el in self.group.values() if el.simple_type == simple_type]

    def get_group_by_res_atom(self, residue, cann):
        elems = [el for el in self.group.values() if residue == el.residue_type]
        groups = []
        for el in elems:
            for ann in el.atom_native_names:
                if cann in ann:
                    groups.append(el)
        return groups


def main():
    test = Config()

    test.read_config()
    # test.overall_table()
    # print(test.get_group_by_type('NH'))
    # print(test.get_group_by_res_atom('ALA', 'CHAl'))
    # print(test.group['MET_E_CH3'])


if __name__ == '__main__':
    main()
