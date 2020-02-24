import sys

import numpy as np

from h5py import File



class hdfAPI(File):
    """docstring for Hdf_API"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    ### Iterators

    def trj_iter(self):
        for key, item in self['/'].items():
            yield key, item

    # use 'acf' or 'ccf' to iterate trough
    # tcf groups in file
    def tcf_iter(self, tcf=''):
        if tcf:
            for _, item in self.trj_iter():
                yield item[tcf]
        else:
            for _, item in self.trj_iter():
                for tcf in ['acf', 'ccf']:
                    return item[tcf]

    def group_iter(self, tcf='', gname='', all=False):
        if not all:
            trj = self._get_zeroTrj()
            if gname in trj[tcf].keys():
                yield trj[tcf][gname]
                return

        groups = [gname] if gname else self.get_groupList()
        for gname in groups:
            if tcf:
                for item in self.tcf_iter(tcf):
                    if gname in item.keys():
                        g = item[gname]
                        yield g
            else:
                for item in self.tcf_iter():
                    if gname in item.keys():
                        g = item[gname]
                        yield g

    ### Getters

    def get_time(self):
        self._get_zeroTrj
        item = self._get_zeroTrj()
        return np.array(item['acf']['t'])

    def get_names(self, tcf, gname):
        items = list(self['/'].values())
        names = items[0][tcf][gname]['names'][()]
        return names

    def get_trjCount(self):
        return len(self.get_trjList())

    def get_groupList(self, tcf='', traj_name=''):
        if traj_name and tcf:
            groups = list(self['/'][traj_name][tcf].keys())
            groups.remove('t')
        elif tcf:
            items = list(self['/'].values())
            groups = list(items[0][tcf].keys())
            groups.remove('t')
        else:
            groups = set(self.get_groupList('acf'))
            groups.update(self.get_groupList('ccf'))
            groups = list(groups)
        return sorted(groups)

    def get_trjList(self):
        return self['/'].keys()

    def get_tcfList(self):
        return self._get_zeroTrj().keys()

    def get_TrjInfo(self):
        pass

    def _get_zeroGroup(self, tcf, gname):
        trj = self._get_zeroTrj()
        if gname in trj[tcf].keys():
            return trj[tcf][gname]
        else:
            print("ERROR!! {} not in {}".format(gname, tcf))

    def _get_zeroTrj(self):
        return list(self['/'].values())[0]

    ### Checkers

    def cmp_t(self):
        arr = self.get_time()
        for item in self.tcf_iter():
            if not np.array_equal(arr, item['t']):
                print("ERROR!! time arrays are not equal!")
                return False
        return True

    def cmp_groups(self):
        try:
            for tcf in ['acf', 'ccf']:
                groups = self.get_groupList(tcf)
                for key, item in self.trj_iter():
                    glist = self.get_groupList(tcf, key)
                    if not groups == glist:
                        print("ERROR!! group names are not equal")
                for gname in groups:
                    # get atoms, cf, gs, names, smarts
                    g = self._get_zeroGroup(tcf, gname)
                    atoms = g['atoms'][()]
                    cf = g['cf']
                    gs = g['group_size']
                    names = g['names'][()]
                    smarts = g['smarts'][()]
                    for group in self.group_iter(tcf, gname):
                        if atoms != group['atoms'][()]:
                            raise Exception("ERROR!! atoms are not equal")
                        if cf.shape != group['cf'].shape:
                            raise Exception("ERROR!! cf are not equal")
                        if gs != group['group_size']:
                            raise Exception("ERROR!! group size are not equal")
                        if names != group['names'][()]:
                            raise Exception("ERROR!! names are not equal")
                        if smarts != group['smarts'][()]:
                            raise Exception("ERROR!! smarts are not equal")
                return True
        except Exception as e:
            return False

    def check_file(self):
        if self.cmp_t() and self.cmp_groups():
            return True
        else:
            return False

    ### Manipulation

    def array_tcf(self, tcf, gname):
        g = self._get_zeroGroup(tcf, gname)
        trjCount = self.get_trjCount()
        arr = np.zeros((trjCount, *g['cf'].shape))
        for i, group in zip(range(trjCount), self.group_iter(tcf, gname, True)):
            arr[i] = group['cf']
        return arr

    def mean_tcf(self, tcf='', gname=''):
        g = self._get_zeroGroup(tcf, gname)
        mean = np.zeros(g['cf'].shape)
        partial_sum_squares = np.zeros(g['cf'].shape)
        trjCount = self.get_trjCount()
        if tcf and gname:
            arr = self.array_tcf(tcf, gname)
        return np.mean(arr, axis=0), np.std(arr, axis=0)

