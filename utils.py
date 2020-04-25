import re
import copy
import numpy as np
import math as m

# line for comments
comment_line_re = re.compile(r'^[\s\w\W]*\;(?P<comment>[\w\d\s]*)$')


def split_comments(file):
    i = 0
    for line in file:
        line = re.sub(r'\;[\w\W\s]*$', '', line).strip()
        if line == '':
            continue
        yield line, i + 1
        i += 1


def not_full_list(list_):
    for elem in list_:
        if elem is None:
            return True
    return False


def not_full_dict(dict_):
    for el in dict_.values():
        if el is None:
            return True
    return False


def join_lists(*lists):
    add = lambda x, y: x if x is not None else y
    all_lists = [None] * len(lists)
    for clist in lists:
        all_lists = list(map(add, all_lists, clist))
    return all_lists


def load_np(md_count, fd):
    for i in range(md_count):
        yield np.load(fd)


def get_ncopy(obj, n):
    for _ in range(n):
        yield copy.deepcopy(obj)


def near_2pow(n):
    power = m.ceil(m.log(n, 2))
    return 2 ** power
