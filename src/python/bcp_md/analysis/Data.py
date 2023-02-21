import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

__all__ = ['Data2D', 'Data_TimeSeries1D', 'Data_TimeSeries2D', 'Data_TimeSeries3D']

class Data2D:

    def __init__(self, z=None, name=None):

        self.z = z
        self.name = name

        self.plot_args = {}

    def plot(self, save=None, show=False, plot_args={}, **kwargs):

        self.plot_args.update(plot_args)
        print(plot_args)
        self._plot(save=save, show=show, **kwargs)

    def _plot(self, save=None, show=False, figsize=(4,3), **kwargs):

        self.fig, self.ax = plt.subplots( figsize=figsize )

        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)

        zmin = plot_args['zmin']
        zmax = plot_args['zmax']
        cmap = plot_args['cmap']
        cbarticks = plot_args['cbarticks']
        cbarlabel = plot_args['cbarlabel']
        xticks = plot_args['xticks']
        xticklabels = plot_args['xticklabels']
        xlabel = plot_args['xlabel']
        yticks = plot_args['yticks']
        yticklabels = plot_args['yticklabels']
        ylabel = plot_args['ylabel']

        self.im = self.ax.imshow(self.z, vmin=zmin, vmax=zmax, cmap=cmap, origin='lower')
        self.ax.set_aspect('auto')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = self.fig.colorbar(self.im, cax=cax)
        cbar.set_ticks(cbarticks)
        cbar.set_label(cbarlabel)
        self.ax.set_xticks(xticks)
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_xlabel(xlabel)
        self.ax.set_yticks(yticks)
        self.ax.set_ylim(yticks[0] - 0.5, yticks[-1] + 0.5)
        self.ax.set_yticklabels(yticklabels)
        self.ax.set_ylabel(ylabel)
        plt.tight_layout(pad=1, h_pad=None, w_pad=None, rect=None)
        
        if save:
            plt.savefig(save, dpi=1000)
        
        if show:
            plt.show()
        
        plt.close(self.fig.number)

class Data_TimeSeries1D:
    def __init__(self, data, Nrepeat, std_true=True):
        
        ndim = data.ndim
        
        if ndim == 1:
            self.Nrepeat = Nrepeat
            self.std_true = std_true
            self.run(data)
        else:
            print("Given array is not 1D")

    def run(self, data):
        data = data.reshape(len(data), -1)
        binsize = 1

        self.initial = data[0]
        self.rest = data[1:].reshape(-1, self.Nrepeat, binsize)
        
        self.mean = self.rest.mean(axis=1)
        self.mean = np.vstack((self.initial, self.mean))
        self.mean = self.mean.reshape(-1)
        
        if self.std_true:
            self.std = self.rest.std(axis=1)        
            self.std = np.vstack((np.zeros(binsize), self.std))
            self.std = self.std.reshape(-1)
            
class Data_TimeSeries2D:
    def __init__(self, data, Nrepeat, t_axis=0, std_true=True):
        
        ndim = data.ndim
        
        if ndim == 2:
            self.Nrepeat = Nrepeat
            self.t_axis = t_axis
            self.std_true = std_true
            self.run(data)
        else:
            print("Given array is not 2D")
        
    def run(self, data):
        binsize = data.shape[1]
        axis = 1
        if self.t_axis == 1:
            binsize = data.shape[0]
            axis = 0

        self.initial = data[0]
        self.rest = data[1:].reshape(-1, self.Nrepeat, binsize)
        
        self.mean = self.rest.mean(axis=axis)
        self.mean = np.vstack((self.initial, self.mean))
        
        if self.std_true:
            self.std = self.rest.std(axis=axis)        
            self.std = np.vstack((np.zeros(binsize), self.std))
        
class Data_TimeSeries3D:
    def __init__(self, data, Nrepeat, t_axis=0, axis=1, bin_axis=2, std_true=True):
        
        ndim = data.ndim
        
        if ndim == 3:
            self.Nrepeat = Nrepeat
            self.t_axis = t_axis
            self.axis = axis
            self.bin_axis = bin_axis
            self.std_true = std_true
            self.run(data)
        else:
            print("Given array is not 3D")

    def run(self, data):
        size = data.shape[self.axis]
        binsize = data.shape[self.bin_axis]
        
        self.initial = data[0]
        self.rest = data[1:].reshape(-1, self.Nrepeat, size, binsize)
        
        self.mean = self.rest.mean(axis=self.axis)
        self.mean = np.vstack((self.initial[None, :, :], self.mean))
        
        if self.std_true:
            self.std = self.rest.std(axis=self.axis)
            self.std = np.vstack((np.zeros((size, binsize))[None, :, :], self.std))
