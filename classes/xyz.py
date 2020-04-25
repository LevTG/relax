import time as t
import numpy as np
import matplotlib
matplotlib.use('Agg')
import multiprocessing.dummy as mp
import matplotlib.pyplot as plt

from numpy.fft import fft, ifft
from utils import near_2pow


class XYZ:
    def __init__(self, md_count):
        self.xyz = [None]*md_count
        self.md_count = md_count