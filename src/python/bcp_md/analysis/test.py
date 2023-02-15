
import numpy as np
import lmfit
import re
import os
import subprocess

__all__ = ['Dump', 'Thermo']

def sigmoid(v, x):
    return v['c']/(1.0 + np.exp(-v['a']*(x - v['b'])))


class Dump:

    def __init__(self, dump_files):
        self.fnames = dump_files
        self.groups = {}
        self.lx = None
        self.ly = None
        self.bins = None

    def set_groups(self, **kwargs):
        self.groups.update(kwargs)

    def set_lx(self, lx):
        self.lx = lx

    def set_ly(self, ly):
        self.ly = ly

    def set_bins(self, bins):
        self.bins = bins
        
    def compute_density(self):

        if self.lx is None: raise Exception("set lx")
        if self.ly is None: raise Exception("set ly")
        if self.bins is None: raise Exception("set bins")

        self._compute_density()

    def _compute_density(self):

        def func2minimize(pars, x, y):
            v = pars.valuesdict()
            m = sigmoid(v, x)

            return m - y 

        for f_ind, f_name in enumerate(self.fnames):
            if f_ind == 0:
                trj = np.loadtxt(f_name, skiprows=9)
            else:
                trj = np.row_stack((trj, np.loadtxt(f_name, skiprows=9)))

        trjs = {}
        for group in self.groups:
            trjs.update({group: np.concatenate([trj[trj[:, 2] == i] for i in self.groups[group]])})

        density = np.zeros([len(self.bins) - 1, len(self.groups) + 1])

        for group_ind, group in enumerate(self.groups, 1):
            hist, bin_edges = np.histogram(trjs[group][:,5], bins=self.bins)
            density[:, 0] = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
            density[:, group_ind] = hist/(self.lx * self.ly * (bin_edges[1] - bin_edges[0]))/len(self.fnames)

        # Fit a sigmoid curve to group1 density
            if group_ind == 1:
                pars = lmfit.Parameters()
                pars.add_many(('a', -0.5), ('b', np.average(density[:, 0])), ('c', 0.8))

                mi = lmfit.minimize(func2minimize, pars, args=(density[3:, 0], density[3:, group_ind]))
                #lmfit.printfuncs.report_fit(mi.params)
                popt = mi.params
                print (popt)
                density_fit = sigmoid(popt, density[:, 0])

        self.density = density


        
        
class Thermo:

    def __init__(self, log_files=None):

        if log_files is None: raise Exception("specify log files")

        self.fnames = log_files
        self.keys = None
        self.available_keys = []

        self._initialize()
    
    def _initialize(self):
        
        with open(self.fnames[0]) as fin:
            for line_ind, line in enumerate(fin.readlines(), 1):
                match_re = re.compile('^Step')
                m = match_re.match(line)
                if m:
                    self.available_keys = line.split(" ")[:-1]

    def set_keys(self, keys=None):
        if keys is not None:
            self.keys = keys
    
    def readin(self):

        if self.keys is None: raise Exception("set thermo keys")

        self._readin()

    def _readin(self):
        outputs = {}

        for f_ind, f_name in enumerate(self.fnames, 1):
            starts = []
            ends = []
            with open(f_name) as fin:
                for line_ind, line in enumerate(fin.readlines(), 1):
                    match_re = re.compile('^Step')
                    m = match_re.match(line)
                    if m:
                        start = line_ind
                        starts.append(start)
                        #print(start)
                    match_re = re.compile('^Loop')
                    m = match_re.match(line)
                    if m:
                        end = line_ind
                        ends.append(end)
                        #print(end)

            for ind, (start, end) in enumerate(zip(starts, ends)):
                subprocess.Popen(f"sed -n '{start},{end-1}p' {f_name} > temp_{ind}", shell=True).wait()
                file = open(f"temp_{ind}", 'r')

                if ind == 0:
                    thermos = file.readline().split(" ")[:-1]

                file.close()

            self.available_keys = thermos

            # Include by default column 0, which corresponds to Step
            columns = [0]
            # Iterate over thermo keys specified beforehand
            for thermo in self.keys:
                column = thermos.index(thermo)
                #print(column)
                columns.append(int(column))

            # Load output files w/o the first row
            print(f"* Read in: {f_name}")
            for ind in range(len(starts)):
                if ind == 0:
                    output = np.loadtxt(f"temp_{ind}", skiprows=1)[:, columns]
                else:
                    output = np.row_stack((output, np.loadtxt(f"temp_{ind}", skiprows=1)[1:, columns]))
                os.remove(f"temp_{ind}")

            if f_ind > 1:
                if output[0, 0] > 0:
                    output[:, 0] -= output[0, 0] 

            # Concatenate output files in the order of increasing Step
            if f_ind == 1:
                outputs = output.copy()
            else:
                output[:, 0] += outputs[-1, 0]
                outputs = np.concatenate((outputs, output[1:, :]))

        self.outputs = outputs


        
        
