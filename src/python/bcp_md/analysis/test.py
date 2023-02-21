import numpy as np
import lmfit
import re
import os
import math
import sys
import glob
import subprocess
import matplotlib.pyplot as plt
from scipy import stats
sys.setrecursionlimit(100000)

from . import Data
from . import Colormaps

__all__ = ['Dump', 'DumpProgress', 'Thermo']

class Dump:
    """
    Attributes:
        dump_files (list)
    """

    def __init__(self, dump_files):
        """
        Dump constructor
        """

        if isinstance(dump_files, str):
            self.fnames = [dump_files]
        elif isinstance(dump_files, list):
            self.fnames = dump_files


    def compute_density(self, groups, bins):
        """
        compute density through the film thickness
        Arguments:
        ----------
        groups (dict): types of beads for each group whose density is computed
        bins (list): bins to be used to compute density

        Returns:
        --------
        N/A
        """
        self.groups = groups
        self.bins = bins

        self._compute_density()


    def _compute_density(self):
        """
        Given groups and bins, compute density
        """

        def func2minimize(pars, x, y):
            v = pars.valuesdict()
            m = sigmoid(v, x)

            return m - y 

        for f_ind, f_name in enumerate(self.fnames):
            if f_ind == 0:
                trj, box = loadtrj(f_name)
                self.lx = box[0][1] - box[0][0]
                self.ly = box[1][1] - box[1][0]
            else:
                trj_, _ = loadtrj(f_name)
                trj = np.row_stack((trj, trj_))

        trjs = {}
        for group_name, group_types in self.groups.items():
            trjs.update({group_name: np.concatenate([trj[trj[:, 2] == i] for i in group_types])})

        density = np.zeros([len(self.bins) - 1, len(self.groups) + 1])

        for group_ind, group_name in enumerate(self.groups, 1):
            hist, bin_edges = np.histogram(trjs.get(group_name)[:,5], bins=self.bins)
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
        

    def compute_pattern(self, types, delta=[0.5, 0.5, 0.5], fft_resolution=2048, zoom_mag=4):
        """
        determine pattern made by beads of given types and compute repeat spacing
        Arguments:
        ----------
        types (list): types of beads to be used to make blobs/morphologies and pattern
        delta (list): binning size for each axis (default: [0.5, 0.5, 0.5])
        fft_resolution (int): 2Dfft (fft_resolution by fft_resolution)
        zoom_mag (int): zoom ratio

        Returns:
        --------
        N/A
        """
        distance_in_pixels = 4
        known_distance = 1
        distance_unit = 'sigma'

        scale = distance_in_pixels / known_distance #pixel/unit
        gridsize = 1.0 / scale #unit/pixel
        delta = [gridsize, gridsize, gridsize]
        
        for f_ind, f_name in enumerate(self.fnames):
            trj, box = loadtrj(f_name)
            
            if f_ind == 0:
                binary_2D = binarize_2D(trj, box, [0.5, 0.5], types) # Work as a mask to sellectively remove dots/lines
                binary_3D = binarize_3D(trj, box, delta, types)
                topview = view_xy(binary_3D, delta, method='topview')#'xray_average')
            else:
                binary_2D += binarize_2D(trj, box, [0.5, 0.5], types) # Work as a mask to sellectively remove dots/lines
                binary_3D = binarize_3D(trj, box, delta, types)
                topview = topview + view_xy(binary_3D, delta, method='topview')#'xray_average')
        
        topview /= len(self.fnames) # Averge over 
        topview /= np.max(topview)

        from scipy.ndimage import gaussian_filter

        binary_2D[binary_2D > 0] = 1
        binary_2D = gaussian_filter(binary_2D, 1)
        binary_2D[binary_2D > 0.5] = 1
        binary_2D[binary_2D < 0.5] = 0

        h, w = binary_2D.shape

        blobs = group2blob(binary_2D)

        lam=[]
        cyl=[]
        for blob_ind, blob in enumerate(blobs):
            #print (len(blob))
            if len(blob) < 10:
                for j in range(len(blob)):
                    binary_2D[blob[j][0]%h, blob[j][1]%w]=0
            else:
                if (len(blob) > 10) and (len(blob) < 650):
                    cyl.append(blob)
                    print (len(blob))
                    #for j in range(len(blob)):
                        #binary_2D[blob[j][0]%h, blob[j][1]%w]=1
                        #topview[blob[j][0]%h*2, blob[j][1]%w*2]=0
                        #topview[blob[j][0]%h*2+1, blob[j][1]%w*2]=0
                        #topview[blob[j][0]%h*2, blob[j][1]%w*2+1]=0
                        #topview[blob[j][0]%h*2+1, blob[j][1]%w*2+1]=0

                if len(blob) > 650:
                    lam.append(blob)
                    print (len(blob))
                    #for j in range(len(blob)):
                        #binary_2D[blob[j][0]%h, blob[j][1]%w]=0
                        #topview[blob[j][0]%h*2, blob[j][1]%w*2]=0
                        #topview[blob[j][0]%h*2+1, blob[j][1]%w*2]=0
                        #topview[blob[j][0]%h*2, blob[j][1]%w*2+1]=0
                        #topview[blob[j][0]%h*2+1, blob[j][1]%w*2+1]=0
        
        N=np.array(topview.shape)
        print ('imarray shape = [' + ', '.join("%d" % iind for iind in N) + ']')

        if N[0] > N[1]:
            im_h = 1
            im_w = N[1]/N[0]*im_h
            dpi = N[0]
        else:
            im_w = 1
            im_h = N[0]/N[1]*im_w
            dpi = N[1]

        #fig = plt.figure(frameon=False)
        #fig.set_size_inches(im_w, im_h)
        #ax = plt.Axes(fig,[0., 0., 1., 1.])
        #ax.set_axis_off()
        #fig.add_axes(ax)
        #ax.imshow(binary_2D, origin='lower')
        #plt.savefig('binary_2D.png', dpi=dpi)
        #plt.show()

        x = {
            'ticks': np.linspace(0,N[1],4),
            'ticklabels': np.linspace(0, N[1]*delta[1], 4, dtype=int),
            'label': r'$x$ ($\sigma$)',
        }

        y = {
            'ticks': np.linspace(0,N[0],4),
            'ticklabels': np.linspace(0, N[0]*delta[0], 4, dtype=int),
            'label': r'$y$ ($\sigma$)',
        }

        topview_image = Data2DImage(topview)
        topview_image.plot(x, y)

        q_del_px = 2*math.pi/fft_resolution # Delta q in pixel
        q_del = q_del_px*scale # Delta q in distance unit

        q = np.linspace(-1.5, 1.5, 7)
        qq_x = fft_resolution//2+q/q_del
        qq_y = fft_resolution//2+q/q_del

        " Apply 2D fft, shift spectrum, and take absolute values "
        fft = np.fft.fft2(topview)
        fshift = np.fft.fftshift(fft)
        magnitude_spectrum_shift = np.abs(fshift)

        from scipy import interpolate
        " Interpolate spectrum over a finer 2D grid (fft_resolution) "
        x = np.linspace(0, N[1]-2, N[1]-1)
        y = np.linspace(0, N[0]-2, N[0]-1)
        z = magnitude_spectrum_shift[1:N[0],1:N[1]]
        f = interpolate.interp2d(x, y, z, kind='linear')

        xnew = np.linspace(0,N[1]-2,fft_resolution+1)
        ynew = np.linspace(0,N[0]-2,fft_resolution+1)
        znew = f(xnew, ynew)

        x = {
            'ticks': qq_x,
            'ticklabels': q,
            'label': r'$q_{x}$ (1/$\sigma$)',
        }

        y = {
            'ticks': qq_y,
            'ticklabels': q,
            'label': r'$q_{y}$ (1/$\sigma$)',
        }

        znew_image = Data2DImage(znew)
        znew_image.plot(x, y, cmap=plt.get_cmap('jet'))


        " Calculate average in the radial direction "
        c_x = znew.shape[1]//(2*zoom_mag)
        c_y = znew.shape[0]//(2*zoom_mag)

        x, y = np.meshgrid(np.arange(fft_resolution//(1*zoom_mag)), np.arange(fft_resolution//(1*zoom_mag)))
        R = np.sqrt((x-c_x)**2 + (y-c_y)**2)

        ll = fft_resolution//2 - fft_resolution//(2*zoom_mag)
        ul = fft_resolution//2 + fft_resolution//(2*zoom_mag)
        znew_zoomed = znew[ll:ul, ll:ul]

        f = lambda r : znew_zoomed[(R >= r-.5) & (R < r+.5)].mean()
        r = np.linspace(1, fft_resolution//(2*zoom_mag), fft_resolution//(2*zoom_mag))
        mean = np.vectorize(f)(r)
        mean *= r**3

        " Fit a gaussian plus line curve to I-q "
        left_ind = 10#30
        right_ind = 60#120
        x = r[left_ind:right_ind]
        y = mean[left_ind:right_ind]

        pars = lmfit.Parameters()
        pars.add_many(('amp', max(y)), ('cen', np.average(x)), ('wid', np.std(x)), ('slope', 0), ('intercept', 0))

        def gaussian_plus_line(v, x):
            """ line + 1-d gaussian """

            gauss = v['amp'] * np.exp(-(x-v['cen'])**2 / v['wid'])
            line = v['slope']*x + v['intercept']
            return gauss + line

        def func2minimize(pars, x, y):

            v = pars.valuesdict()
            m = gaussian_plus_line(v, x)

            return m - y

        mi = lmfit.minimize(func2minimize, pars, args=(x, y))
        lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)

        print(lmfit.fit_report(mi))

        q0 = mi.params['cen'] * q_del_px * scale
        q0_plus_std = (mi.params['cen'] + mi.params['cen'].stderr) * q_del_px * scale
        l0 = 2 * math.pi / q0
        l0_minus_sigma = 2 * math.pi / q0_plus_std
        print ('L0 = %.4f' % l0)


    def compute_localfraction(self, types1, types2):
        """
        compute local fraction of a component in the blobs made by blocks of a certain type
        the complement set of types1 and types2 over types 2
        Arguments:
        ----------
        types1 (list): types of a component
        types2 (list): types of a type of block

        Returns:
        --------
        N/A
        """
        delta = 0.75

        complement = list(set(types1) & set(types2))
        if len(complement) == 0: raise Exception("set types right")

        for f_ind, f_name in enumerate(self.fnames):
            trj, box = loadtrj(f_name)

            if f_ind == 0:
                binary2D_1 = binarize_2D(trj, box, [delta, delta], complement)
                binary2D_2 = binarize_2D(trj, box, [delta, delta], types2)
            else:
                binary2D_1 += binarize_2D(trj, box, [delta, delta], complement)
                binary2D_2 += binarize_2D(trj, box, [delta, delta], types2)

        fraction = binary2D_1 / binary2D_2

        N = np.array(fraction.shape)

        x = {
            'ticks': np.linspace(0, N[1]-1, 4),
            'ticklabels': np.linspace(0, N[1]*delta, 4, dtype=int),
            'label': r'$x$ ($\sigma$)',
        }

        y = {
            'ticks': np.linspace(0, N[0]-1, 4),
            'ticklabels': np.linspace(0, N[0]*delta, 4, dtype=int),
            'label': r'$y$ ($\sigma$)',
        }

        from matplotlib.colors import ListedColormap
        cmaps = {}
        G2Y = np.dstack((np.linspace(0,1,256), np.linspace(100/255,1,256), np.linspace(0,0,256)))
        cmaps['G2Y'] = ListedColormap(G2Y[0], name='G2Y')

        fraction_image = Data2DImage(fraction)
        fraction_image.plot(x, y, cbar_params={'ticks': [0, 1], 'label': 'fraction'}, cmap=cmaps['G2Y'], facecolor='black')

        self.localfraction = fraction


    def compute_perobject(self, types_A, types_B, periodic=False, obj_flag=False):
        """
        1. identify every object (blob)
        2. identify chains associated with each object and compute how far each bead is
            located from the interface (in both normal and parallel directions)

        Only after every object is identified, change obj_flag to True to proceed
        For infinitely long object to be handled, change periodic to True

        Arguments:
        ----------
        types_A (list): types of beads constructing objects
        types_B (list): types of beads constructing matrix

        Returns:
        --------
        N/A
        """
        
        from scipy.ndimage import gaussian_filter
        from skimage import measure

        delta=0.75

        print("this uses only one dump file")
        f_name = self.fnames[-1]
        
        trj, box = loadtrj(f_name)

        trj9 = np.zeros([len(trj)*9,len(trj[0])])
        for i in range(9):
            trj9[len(trj)*i:len(trj)*(i+1),:] = trj
        trj9[:len(trj),3] = trj[:,3] - box[0][1] #
        trj9[:len(trj),4] = trj[:,4] + box[1][1]
        trj9[len(trj):len(trj)*2,3] = trj[:,3]
        trj9[len(trj):len(trj)*2,4] = trj[:,4] + box[1][1]
        trj9[len(trj)*2:len(trj)*3,3] = trj[:,3] + box[0][1]
        trj9[len(trj)*2:len(trj)*3,4] = trj[:,4] + box[1][1]
        trj9[len(trj)*3:len(trj)*4,3] = trj[:,3] - box[0][1] #
        trj9[len(trj)*3:len(trj)*4,4] = trj[:,4]
        trj9[len(trj)*4:len(trj)*5,3] = trj[:,3]
        trj9[len(trj)*4:len(trj)*5,4] = trj[:,4]
        trj9[len(trj)*5:len(trj)*6,3] = trj[:,3] + box[0][1]
        trj9[len(trj)*5:len(trj)*6,4] = trj[:,4]
        trj9[len(trj)*6:len(trj)*7,3] = trj[:,3] - box[0][1] #
        trj9[len(trj)*6:len(trj)*7,4] = trj[:,4] - box[1][1]
        trj9[len(trj)*7:len(trj)*8,3] = trj[:,3]
        trj9[len(trj)*7:len(trj)*8,4] = trj[:,4] - box[1][1]
        trj9[len(trj)*8:len(trj)*9,3] = trj[:,3] + box[0][1]
        trj9[len(trj)*8:len(trj)*9,4] = trj[:,4] - box[1][1]

        for ind, tp in enumerate(types_A):
            if ind == 0:
                logic = trj[:, 2] == tp
            else:
                logic = np.logical_or(logic, trj[:, 2] == tp)
        trj_filter = trj[logic]

        binary2D = binarize_2D(trj, box, [delta, delta], types_A)
        binary2D[binary2D > 0] = 1

        img = Data2DImage(binary2D)
        img.plot()

        h, w = binary2D.shape

        blobs, blob_tmp = group2blob(binary2D)

        lam = []
        cyl = []
        for i in range(len(blobs)):
            if len(blobs[i]) < 10:
                for j in range(len(blobs[i])):
                    binary2D[blobs[i][j][0]%h, blobs[i][j][1]%w] = 0
            else:
                if (len(blobs[i]) > 37.5) and (len(blobs[i]) < 260):
                    print ('cyl', len(blobs[i]))
                    cyl.append(blobs[i])

                if len(blobs[i]) > 260:
                    print ('lam', len(blobs[i]))
                    lam.append(blobs[i])
        cyl_and_lam = cyl + lam

        mols = np.unique(trj[:,1]) # list of all mol-IDs in the given trajectory
        mols_morph = [] # list of the lists of molecules constructing each object
        for iind in range(len(cyl_and_lam)):
            mols_morph_temp = []

            for jind in range(len(mols)):
                kind = 0
                flag = 0

                blockA = trj_filter[trj_filter[:,1] == mols[jind]] # trj of beads contructing block A of each chain

                while (flag == 0 and kind < len(blockA)):
                    bins = trj2bin(blockA[kind], box, [delta, delta])
                    if [int(round(bins[0])), int(round(bins[1]))] in cyl_and_lam[iind]:
                        flag = 1
                        mols_morph_temp.append(blockA[kind, 1])

                    kind = kind + 1

            mols_morph.append(mols_morph_temp)
        
        plt.figure(num=2)
        plt.imshow(binary2D,cmap=plt.get_cmap('gray'),vmin=0,vmax=1,origin='lower')
        plt.gca().set_aspect('equal',adjustable='box')
        plt.show()

        N = np.array(binary2D.shape)

        if N[0] > N[1]:
            im_h = 1
            im_w = N[1]/N[0]*im_h
            dpi = N[0]
        else:
            im_w = 1
            im_h = N[0]/N[1]*im_w
            dpi = N[1]

        #fig = plt.figure(num=5)
        #fig.set_size_inches(im_w, im_h)
        #ax = plt.Axes(fig,[0.,0.,1.,1.])
        #ax.set_axis_off()
        #fig.add_axes(ax)
        fig, ax = plt.subplots()
        ax.imshow(binary2D, cmap=plt.get_cmap('gray'), vmin=0, vmax=1, origin='lower')

        '''
        binary reproduced with contours
        majority (matrix) = 0, minority (morph) = 1 while majority has 1 in binary
        '''
        reproduced_binary = np.ones(binary2D.shape)

        moi = cyl_and_lam #morphology of interest
        contours = []
        for moi_id in range(len(moi)):
        #for moi_id in [0]:
            moi_array = np.array(moi[moi_id])

            moi_lx = abs(max(moi_array[:,0]) - min(moi_array[:,0])) + 1
            moi_ly = abs(max(moi_array[:,1]) - min(moi_array[:,1])) + 1
            #margin_x=int(moi_lx*0.1) #margin proportional to lx
            margin_x = 5 #margin fixed regardless of lx
            #margin_y=int(moi_ly*0.1) #margin proportional to ly
            margin_y = 5 #margin fixed regardless of ly

            #a window for each moi to draw a contour over pbc
            moi_window = np.zeros([moi_lx + margin_x*2, moi_ly + margin_y*2])
            vec_trs = np.array([margin_x - min(moi_array[:,0]), margin_y - min(moi_array[:,1])]) #translation vector

            for iind in range(len(moi_array)):
                moi_window[moi_array[iind,0] + vec_trs[0], moi_array[iind,1] + vec_trs[1]] = 1

            if periodic == False:
                # when every object if finite

                moi_gaussian = gaussian_filter(moi_window, sigma=1.1)
                contour = measure.find_contours(moi_gaussian, 0.5) # 0.6 by default, the smaller the bigger contour
                contour = contour[0] - vec_trs
                contour = np.round(contour).astype(int)
                contours.append(contour)

                ax.plot(contour[:,1]%w, contour[:,0]%h,' rs', ms=1) #draw contour discontinued over pbc
                #ax.plot(contour[:,1],contour[:,0],'rs',ms=1) #draw contour continued over pbc

                reproduced_binary[contour[:,0]%h, contour[:,1]%w] = 0.5

                minority_filter = trj_filter[trj_filter[:,1] == mols_morph[moi_id][0]]
                bins = trj2bin(minority_filter[1], box, [delta,delta]) #change index if filling does not work properly

                blob_tmp = []
                floodfill_pbc(reproduced_binary, int(round(bins[0])), int(round(bins[1])), blob_tmp) #floodfill within a contour line, which is a close loop

            else:
                # when an object is infinitely long in either direction

                minority_filter = trj_filter[trj_filter[:,1] == mols_morph[moi_id][0]]
                bins = trj2bin(minority_filter[1], box, [delta,delta]) #change index if filling does not work properly

                test_contour = []
                contour_finder(binary2D, int(round(bins[0])), int(round(bins[1])), test_contour)
                test_contour = np.asarray(test_contour)

                ax.plot(test_contour[:,1]%w, test_contour[:,0]%h, 'mo', ms=3)

                contours.append(test_contour)
                reproduced_binary[test_contour[:,0]%h, test_contour[:,1]%w] = 0.5

                floodfill_pbc(reproduced_binary, int(round(bins[0])), int(round(bins[1])), blob_tmp) #floodfill within a contour line, which is a close loop
                ax.plot(int(round(bins[1])), int(round(bins[0])), 'bd', ms=4)

        #plt.gca().set_aspect('equal',adjustable='box')
        #plt.savefig('MD_wo.png',dpi=dpi)
        plt.show()

        plt.figure(num=5)
        plt.imshow(reproduced_binary, cmap=plt.get_cmap('gray'), vmin=0, vmax=1, origin='lower')
        plt.gca().set_aspect('equal',adjustable='box')
        plt.show()

        if obj_flag == True:
            binary9 = np.zeros([3*h,3*w]) #replicate 3 times in x and y directions 
            binary9[0*h:1*h,0*w:1*w] = reproduced_binary
            binary9[1*h:2*h,0*w:1*w] = reproduced_binary
            binary9[2*h:3*h,0*w:1*w] = reproduced_binary
            binary9[0*h:1*h,1*w:2*w] = reproduced_binary
            binary9[1*h:2*h,1*w:2*w] = reproduced_binary
            binary9[2*h:3*h,1*w:2*w] = reproduced_binary
            binary9[0*h:1*h,2*w:3*w] = reproduced_binary
            binary9[1*h:2*h,2*w:3*w] = reproduced_binary
            binary9[2*h:3*h,2*w:3*w] = reproduced_binary

            gaussian_binary = gaussian_filter(binary9, sigma=0.9) # 1.1 / 0.9
            del binary9

            dif_x = np.diff(gaussian_binary, axis=1)
            line = np.zeros([gaussian_binary.shape[0], 1])
            dif_x = np.concatenate((dif_x, line), axis=1)

            dif_y = np.diff(gaussian_binary, axis=0)
            line = np.zeros([1, gaussian_binary.shape[1]])
            dif_y = np.concatenate((dif_y, line), axis=0)

            numerator = 2.0*np.multiply(dif_x,dif_y)
            denominator = np.power(dif_x,2) - np.power(dif_y,2)

            angles = 0.5*np.arctan2(numerator,denominator)
            angles = angles[1*h:2*h, 1*w:2*w]

            plt.figure(num=7)
            plt.imshow(reproduced_binary, cmap=plt.get_cmap('gray'), vmin=0, vmax=1, origin='lower')

            distances_final=[]
            rotated_final = []

            # for each morphology
            for moi_id in range(len(moi)):
            #for moi_id in [1]:

                contour = contours[moi_id]
                mol_track = []
                distances = []
                rotated = []

                for probeind in range(len(contour)):
                    angle = angles[contour[probeind,0]%h, contour[probeind,1]%w]
                    ray = np.array([math.cos(angle), math.sin(angle)])
                    unitray = ray/np.linalg.norm(ray)

                    pt_on_trj=np.array([(contour[probeind,1]%w)*delta, (contour[probeind,0]%h)*delta])

                    # change the magnitude in order to make vectors toward majority
                    mag = 4.5
                    if reproduced_binary[int((contour[probeind,0]+unitray[1]*mag)%h), int((contour[probeind,1]+unitray[0]*mag)%w)] != 1:
                        # does not fall within majority, take the opposite direction
                        unitray = unitray*(-1)

                    plt.plot([int((contour[probeind,1]+unitray[0]*mag)%w)], [int((contour[probeind,0]+unitray[1]*mag)%h)], 'rd', ms=3, alpha=0.25)

                    d1 = 2.0
                    d2i = 5.0
                    d2 = 5.0
                    d_parallel = np.cross(unitray, trj9[:,3:5] - pt_on_trj)
                    trj9_by_d1 = trj9[np.abs(d_parallel) < d1]
                    d_normal = np.dot(trj9_by_d1[:,3:5] - pt_on_trj, unitray)

                    mol_id = {}
                    mol_id_temp=[]
                    for iind in range(len(trj9_by_d1)):
                        if d_normal[iind] > -d2i/2.0 and d_normal[iind] < d2i/2.0:

                            for tp_ind, tp in enumerate(types_B):
                                if tp_ind == 0:
                                    logic_B = trj9_by_d1[iind, 2] == tp
                                else:
                                    logic_B = np.logical_or(logic_B, trj9_by_d1[iind, 2] == tp)

                            for tp_ind, tp in enumerate(types_A):
                                if tp_ind == 0:
                                    logic_A = trj[int(trj9_by_d1[iind, 0])-2, 2] == tp
                                else:
                                    logic_A = np.logical_or(logic_A, trj[int(trj9_by_d1[iind, 0])-2, 2] == tp)

                            if logic_B and logic_A:
                                if trj9_by_d1[iind,1] not in mol_track:
                                    mol_id_temp.append(trj9_by_d1[iind,1])
                                    mol_track.append(trj9_by_d1[iind,1])                                

                    mol_id.update({'interest': np.unique(np.array(mol_id_temp))})

                    for iind in range(len(mol_id['interest'])):

                        # reconstruct a chain in case it is over boundaries
                        chain_rep = trj[trj[:,1] == mol_id['interest'][iind]]
                        chain = trj[trj[:,1] == mol_id['interest'][iind]]
                        for jind in range(1,len(chain)):
                            chain[jind,3:5] = chain[jind-1,3:5] + d_pbc2D(chain_rep[jind,3:5], chain_rep[jind-1,3:5], [box[0][1]-box[0][0], box[1][1]-box[1][0]])

                        types = np.unique(chain[:,2])
                        idx_bead2 = sum(chain[:,2] == types[0])
                        idx_bead1 = idx_bead2 - 1

                        'use a vector connecting two blocks'
                        vector_inter = chain[idx_bead2, 3:6] - chain[idx_bead1, 3:6]
                        unitvector_inter = vector_inter[:2]/np.linalg.norm(vector_inter[:2])

                        center = chain[idx_bead1,3:6] + vector_inter*0.5

                        'use a normal vector at the interface'
                        angle_n = angles[int((center[1]/delta)%h), int((center[0]/delta)%w)]
                        unitvector_n = np.array([math.cos(angle_n), math.sin(angle_n)])

                        if reproduced_binary[int((center[1]/delta + unitvector_n[1]*mag)%h), int((center[0]/delta + unitvector_n[0]*mag)%w)] != 1:
                            unitvector_n = unitvector_n*(-1)

                        y_new = unitvector_n
                        z_new = np.array([0, 0, 1])
                        x_new = np.cross(y_new, z_new)

                        x = np.array([1, 0, 0])
                        y = np.array([0, 1, 0])

                        theta = np.arccos(np.dot(x, x_new))
                        sign = 1 if np.cross(x, x_new)[2] > 0 else -1
                        theta *= sign

                        '''
                        distance array
                        0: atom-ID, 1: atom-type, 2: mol-ID, 3: distance normal to interface, 4: distance parallel to interface
                        '''
                        distance = np.zeros([len(chain),5])
                        distance[:,:3] = chain[:,:3]
                        distance[:,3] = np.dot((chain[:,3:5] - center[:2]), unitvector_n)
                        distance[:,4] = np.cross((chain[:,3:5] - center[:2]), unitvector_n)

                        distances.append(distance)

                        '''
                        rotated distance array
                        as if the normal vector is y-axis
                        '''
                        rtt = np.zeros([len(chain),6])
                        rtt[:,:3] = chain[:,:3]
                        rtt_matrix = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
                        rtt[:,3:5] = np.array([np.matmul(rtt_matrix, chain[i,3:5] - center[:2]) for i in range(len(chain))])
                        rtt[:,5] = chain[:,5] - center[2]

                        rotated.append(rtt)

                        #plt.plot([center[0]/delta, center[0]/delta + unitvector_n[0]*mag], [center[1]/delta, center[1]/delta + unitvector_n[1]*mag], 'm*-', ms=4)
                        #plt.plot(center[0]/delta, center[1]/delta, 'g*', ms=4)
                        plt.arrow(center[0]/delta, center[1]/delta, unitvector_n[0]*mag, unitvector_n[1]*mag, head_width=1)

                distances_final.append(distances)
                rotated_final.append(rotated)

            plt.plot(chain[:,3]/delta,chain[:,4]/delta,'bo',ms=3)
            plt.plot(center[0]/delta, center[1]/delta, 'r*', ms=7)
            plt.plot(int((center[0]/delta)%w), int((center[1]/delta)%h), 'ms', ms=4)
            print (unitvector_n)
            #plt.plot([center[0]/delta, center[0]/delta + unitvector_n[0]*5], [center[1]/delta, center[1]/delta + unitvector_n[1]*5], 'm*-', ms=4)
            plt.arrow(center[0]/delta, center[1]/delta, unitvector_n[0]*mag, unitvector_n[0]*mag, head_width=2, facecolor='red')
            plt.gca().set_aspect('equal',adjustable='box')
            plt.show()

            distance = {
                    'dot': [],
                    'line': [],
                    }

            for ind, val in enumerate(distances_final):
                if ind < len(cyl):
                    distance['dot'].append(val)
                else:
                    distance['line'].append(val)

            np.save(f"distance_normal_for_{f_name}", distance)


    def convert_to_lmpdata(self, data_file=None):
        """
        Extract chain topology from an lmp data file and write a new data file using a dump file
        Arguments:
        ----------
        data_file (str): name of data file from which chain topology is extracted

        Returns:
        --------
        N/A
        """
        if len(self.fnames) != 1: raise Exception("use only one dump file")
        if data_file is None: raise Exception("specify lmp data to extract chain topology")

        dump_file = self.fnames[0]

        box = []
        f = open(self.fnames, 'r')
        for iind in range(5):
            f.readline()
        for iind in range(3):
            line = f.readline().split()
            box.append([float(line[0]), float(line[1])])
        f.close()

        atoms = np.loadtxt(dump_file, skiprows=9)
        atoms = atoms[np.argsort(atoms[:,0])]
        atoms = atoms[atoms[:, 2] < 5]

        flag = 0 

        f = open(data_file, 'r')
        while flag == 0:
            line = f.readline()

            match = re.search('^\d+ atoms', line)
            if match:
                n_atoms = int(line.split()[0])

            match = re.search('^\d+ atom types', line)
            if match:
                atom_types = int(line.split()[0])
                masses = np.empty((atom_types, 2), dtype=int)

            match = re.search('^\d+ bonds', line)
            if match:
                n_bonds = int(line.split()[0])

            match = re.search('^\d+ bond types', line)
            if match:
                bond_types = int(line.split()[0])

            match = re.search('^Masses', line)
            if match:
                line = f.readline()
                for atom_type in range(atom_types):
                    line = f.readline()
                    masses[atom_type, 0] = int(line.split()[0])
                    masses[atom_type, 1] = int(line.split()[1])

            match = re.search('^Atoms', line)
            if match:
                flag = 1

        f.close()

        print('* # of bonds = {}'.format(n_bonds))

        subprocess.Popen('grep -A {} "Atoms" {} > temp.txt'.format(n_atoms + 1, data_file), shell=True).wait()
        atoms_tmp = np.loadtxt('temp.txt', skiprows=2)[:,:6]
        os.remove('temp.txt')
        substrate = atoms_tmp[atoms_tmp[:, 2] > np.max(atoms[:,2])]
        atoms = np.row_stack((atoms, substrate))
        subprocess.Popen('grep -A {} "Bonds" {} > temp.txt'.format(n_bonds + 1, data_file), shell=True).wait()
        bonds = np.loadtxt('temp.txt', skiprows=2)
        os.remove('temp.txt')

        f = open('data_from_dump.txt', 'w')
        f.write('# data generated from dump name: {}\n\n'.format(dump_file))
        f.write('%d atoms\n' % atoms.shape[0])
        f.write('%d atom types\n' % atom_types)
        f.write('%d bonds\n' % bonds.shape[0])
        f.write('%d bond types\n' % bond_types)
        f.write('\n')
        f.write('%s %s xlo xhi\n' % (box[0][0], box[0][1]))
        f.write('%s %s ylo yhi\n' % (box[1][0], box[1][1]))
        f.write('%s %s zlo zhi\n' % (box[2][0], box[2][1]))
        f.write('\n')
        f.write('Masses\n\n')
        for mass in masses:
            f.write('%d %d\n' % (mass[0], mass[1]))
        f.write('\n')
        f.write('Atoms\n\n')
        for atom in atoms:
            f.write('%d\t%d\t%d\t%s\t%s\t%s\n' % (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5]))
        f.write('\n')
        f.write('Bonds\n\n')
        for bond in bonds:
            f.write('%d\t%d\t%d\t%d\n' % (bond[0], bond[1], bond[2], bond[3]))
        f.close()


class Data2DImage:

    def __init__(self, data):
        self.data = data

    def plot(self, x=None, y=None, cbar_params=None, cmap=None, show=False, save=False, name='result.png', dpi=200, **kwargs):

        if cmap is not None:
            cmap = cmap
        else:
            cmap = plt.get_cmap('gray')

        fig, ax = plt.subplots(figsize=(4,3), dpi=200)
        ax.set_aspect('equal', adjustable='box')

        if 'facecolor' in kwargs:
            ax.set_facecolor(kwargs.get('facecolor'))

        im = ax.imshow(self.data, cmap=cmap, origin='lower')

        if cbar_params:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_ticks(cbar_params.get('ticks'))
            cbar.set_label(cbar_params.get('label'))

        if x is not None:
            xticks = x['ticks']
            xticklabels = x['ticklabels']
            xlabel = x['label']
            ax.set_xlim(xticks[0], xticks[-1])
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_xlabel(xlabel)

        if y is not None:
            yticks = y['ticks']
            yticklabels = y['ticklabels']
            ylabel = y['label']
            ax.set_ylim(yticks[0], yticks[-1])
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            ax.set_ylabel(ylabel)

        plt.tight_layout(pad=0.5, h_pad=None, w_pad=None, rect=None)

        if save:
            plt.savefig(name, dpi=dpi)
        if show:
            plt.show()

    #def fft(self, fft_resolution=2048, zoom_mag=4):


class DumpProgress:
    """
    Attributes:
        path: path to the directory of simulation results
        dirname_pattern: (default: 'equal_{}')
        dir_start: from this equil number
        dir_end: to this equil number
        fname_pattern: (default: 'dump*.{}')
        Nevery: use input values every this many steps
        Nrepeat: number of steps to use input values for computing mean values
        Nfreq: compute mean values every this many steps
        end: from 0 to this
        timestep: timestep used for the simulation (default: 0.006tau)
        comm: mpi communicator
    """
    def __init__(self, path=None, dirname_pattern='equil_{}', dir_start=0, dir_end=0, fname_pattern='dump*.{}', Nevery=10000, Nrepeat=1, Nfreq=10000, end=10000, timestep=0.006, comm=None):
        """
        DumpProgress constructor
        """
        if path == None:
            print('path was not given and set to current working directory.')
            path = os.getcwd()

        if dir_start > dir_end:
            print('dir_start cannot be larger than dir_end. dir_start is set to 0.')
            dir_start = 0

        if Nevery > Nfreq:
            print('Nevery cannot be larger than Nfreq. Nevery is set to Nfreq.')
            Nevery = Nfreq

        if Nevery > end: 
            print('Nevery cannot be larger than end. Nevery is set to end.')
            Nevery = end

        if Nfreq > end:
            print('Nfreq cannot larger than end. Nfreq is set to end.')
            Nfreq = end

        self.path = path
        self.dirname_pattern = dirname_pattern
        self.dir_start = dir_start
        self.dir_end = dir_end
        self.fname_pattern = fname_pattern
        self.Nevery = Nevery
        self.Nrepeat = Nrepeat
        self.Nfreq = Nfreq
        self.end = end
        self.comm = comm

        self.timestep = timestep
        self.results = {
                'timestep': timestep
                }

        self.source_dirs = []
        for iind in range(dir_start, dir_end + 1):
            self.source_dirs.append(dirname_pattern.format(iind))

        self.freqs = np.linspace(Nfreq, end, int(end/Nfreq), dtype=int)

        fnames = []
        if path is not None:
            fnames.append(glob.glob(os.path.join(path, self.source_dirs[0], self.fname_pattern.format(0)))[0])
        else:
            path = os.getcwd()
            fnames.append(glob.glob(os.path.join(path, self.source_dirs[0], self.fname_pattern.format(0)))[0])
        buff_path = path

        steps = [0]
        for dir_ind, source_dir in enumerate(self.source_dirs):

            for freq_ind, freq in enumerate(self.freqs):

                patterns = [fname_pattern.format(i) for i in np.linspace(freq - self.Nevery*(self.Nrepeat -1), freq, self.Nrepeat, dtype=int)]

                for pattern_ind, pattern in enumerate(patterns):
                    if path is not None:
                        fnames.append(glob.glob(os.path.join(path, source_dir, pattern))[0])
                    else:
                        fnames.append(glob.glob(os.path.join(source_dir, pattern))[0])

                steps.append((dir_start + dir_ind)*end + freq)

        self.fnames = fnames
        self.steps = np.asarray(steps)
        self.results.update({'steps': np.asarray(steps)})

        rank = self.comm.Get_rank()

        if rank == 0:
            print("The number of files to be processed = {}".format(len(fnames)))

            f = open(self.fnames[0], 'r')
            box = []
            for i in range(5):
                f.readline()
            for i in range(3):
                line = f.readline().split()
                box.append([float(line[0]), float(line[1])])
            f.close()

            lx = box[0][1] - box[0][0]
            ly = box[1][1] - box[1][0]

            trj = np.loadtxt(self.fnames[-1], skiprows=9)
            trj = trj[trj[:,2] <= np.unique(trj[:,2])[:-2][-1]]
            bins = np.linspace(0, 80, 161)
            hist, bin_edges = np.histogram(trj[:,5], bins=bins)
            z = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
            rho = hist/(lx * ly * (bin_edges[1] - bin_edges[0]))

            def func2minimize(pars, x, y):
                v = pars.valuesdict()
                m = sigmoid(v, x)

                return m - y

            pars = lmfit.Parameters()
            pars.add_many(('a', -0.5), ('b', np.average(z)), ('c', 0.8))
            mi = lmfit.minimize(func2minimize, pars, args=(z, rho))
            popt = mi.params
            h = mi.params['b']

            self.box = box
            self.lx = lx
            self.ly = ly
            self.h = h
            self.results.update({"height": h.value})

        else:
            self.box = None
            self.lx = None
            self.ly = None

        self.box = comm.bcast(self.box, root=0)
        self.lx = comm.bcast(self.lx, root=0)
        self.ly = comm.bcast(self.ly, root=0)        


    def save_results(self, path=None):
        """
        save the results dict in the given path

        Arguments:
        ----------
        path: path to directory in which the results dict is saved

        Returns:
        ----------
        N/A
        """
        
        if self.comm.rank == 0:
            if path is not None:
                self.path = path

            if os.path.exists(os.path.join(self.path, 'results.npy')):
                results = np.load(os.path.join(self.path, 'results.npy'), allow_pickle=True).item()
            else:
                results = {}

            keys = self.results.keys()
            keys_old = results.keys()
            for key in keys:
                if 'timestep' == key:
                    value = self.results.get('timestep')
                    results.update({'timestep': value})

                elif 'steps' == key:
                    value = self.results.get('steps')
                    if 'steps' in keys_old:
                        value_old = results.get('steps')

                        if value.shape[0] > value_old.shape[0]:
                            results.update({'steps': value})
                    else:
                        results.update({'steps': value})

                else:
                    results.update({key: self.results.get(key)})

            np.save(os.path.join(self.path, 'results.npy'), results)
        

    def specify_timesteps(self, timesteps):
        """
        specify timesteps by manually providing an array

        Arguments:
        ----------
        timesteps (array): use these timesteps to compute observables

        Returns:
        ----------
        N/A
        """

        self.Nrepeat = 1
        fnames = []
        for dir_ind, source_dir in enumerate(self.source_dirs):

            patterns = [self.fname_pattern.format(int(timestep)) for timestep in timesteps]

            for pattern_ind, pattern in enumerate(patterns):
                if self.path is not None:
                    fnames.append(glob.glob(os.path.join(self.path, source_dir, pattern))[0])
                else:
                    fnames.append(glob.glob(os.path.join(source_dir, pattern))[0])

        self.fnames = fnames
        print("The number of files to be processed = {}".format(len(fnames)))    


    def fix_blankdump(self, f_ind):
        """
        fix blank dump files by replacing it by an adjacent one

        Arguments:
        ----------
        f_ind (int): index of self.fnames

        Returns:
        ----------
        trj (array): trajectory array that replaces a broken dump file
        """
        try:

            trj = np.loadtxt(self.fnames[f_ind], skiprows=9)

            f = open(self.fnames[f_ind],'r')
            for i in range(3):
                f.readline()
            n_atoms = int(f.readline().split()[0])
            f.close()

            flag = 0
            shift = 1
            while (flag == 0):
                if trj.shape[0] != n_atoms:
                    trj = np.loadtxt(self.fnames[f_ind - shift], skiprows=9)

                    f = open(self.fnames[f_ind - shift],'r')
                    for i in range(3):
                        f.readline()
                    n_atoms = int(f.readline().split()[0])
                    f.close()

                    shift +=1

                else:
                    flag = 1

        except:
            'empty dump files can be created due to storage limit'

            flag = 0
            shift = 1
            while (flag == 0):

                try:
                    trj = np.loadtxt(self.fnames[f_ind - shift], skiprows=9)

                    f = open(self.fnames[f_ind - shift],'r')
                    for i in range(3):
                        f.readline()
                    n_atoms = int(f.readline().split()[0])
                    f.close()

                    if trj.shape[0] == n_atoms:
                        flag = 1
                    else:
                        shift += 1

                except:
                    shift += 1

        return trj


    def compute_density(self, types, zlo, zhi, n_bins):
        """
        compute density through the film thickness

        Arguments:
        ----------
        types: types of beads to be used to compute density
        zlo: lower limit
        zhi: upper limit
        n_bins: number of bins

        Returns:
        ----------
        N/A
        """

        self.zlo = zlo
        self.zhi = zhi
        self.bins = np.linspace(zlo, zhi, n_bins)

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        avg_rows_per_process = int(len(self.fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1: 
            end_row = len(self.fnames)

        density_tmp = np.empty([len(self.fnames), n_bins - 1])
        for iind in range(start_row, end_row):

            trj = self.fix_blankdump(iind)

            density_tmp[iind, :] = self._computeDensity(trj, types)

        del trj

        if rank == 0:

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1: 
                    end_row = len(self.fnames)

                recv = np.empty([len(self.fnames), n_bins - 1])
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                density_tmp[start_row:end_row] = recv[start_row:end_row]
        else:
            send = density_tmp
            req = self.comm.Isend(send, dest=0)
            req.Wait()

        if rank == 0:

            result = {}

            rho = Data.Data_TimeSeries2D(density_tmp, self.Nrepeat)

            density_final = rho.mean
            density_final = density_final.transpose()

            result.update({'mean': rho.mean.transpose()})
            result.update({'std': rho.std.transpose()})

            res = Data.Data2D(z=density_final)
            xticks = np.linspace(0,len(self.source_dirs)*len(self.freqs), 4)
            xticklabels = [format(i*self.Nfreq*0.006*pow(10,-6), '.4f') for i in xticks]
            plot_args = {
                'zmin': 0,
                'zmax': 0.8,
                'cmap': Colormaps.cmaps['W2B_8'],
                'cbarticks': [0, 0.8],
                'cbarlabel': r'density ($m/\sigma^{3}$)',
                'xticks': xticks,
                'xticklabels': xticklabels,
                'xlabel': r'time ($\times 10^{6} \tau$)',
                'yticks': np.linspace(0, n_bins-1, 5),
                'yticklabels': np.linspace(self.bins[0], self.bins[-1], 5),
                'ylabel': r'$z$ ($\sigma$)',
            }

            res.plot(save='density_evol.png', show=False, plot_args=plot_args)

            t = np.linspace(0, density_final.shape[1]-1, density_final.shape[1])
            tt, zz = np.meshgrid(t - 0.5, np.linspace(0, n_bins-2, n_bins-1))

            result.update({'tt': tt})
            result.update({'zz': zz})
            result.update({'xticks': xticks})
            result.update({'xticklabels': xticklabels})

            self.results.update({'density': result})


    def _computeDensity(self, trj, types):
        """
        Given types and bins, compute density
        """
        for ind, tp in enumerate(types):
            if ind == 0:
                logic = trj[:,2] == tp
            else:
                logic = np.logical_or(logic, trj[:,2] == tp)

        trj_filter = trj[logic]

        hist, bin_edges = np.histogram(trj_filter[:,5], bins=self.bins)
        density = hist/(self.lx*self.ly*(bin_edges[1] - bin_edges[0]))

        return density


    def compute_localfraction(self, types1, types2, zlo, zhi, n_bins):
        """
        compute local fraction of a component in the blobs made by blocks of a certain type
        the complement set of types1 and types2 over types 2

        Arguments:
        ----------
        types1 (list): types of a component
        types2 (list): types of a type of block
        zlo: lower limit
        zhi: upper limit
        n_bins: number of bins

        Returns:
        ----------
        N/A
        """

        self.bins = np.linspace(zlo, zhi, n_bins)

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        avg_rows_per_process = int(len(self.fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1: 
            end_row = len(self.fnames)

        fraction_tmp = np.empty([len(self.fnames), n_bins - 1])
        for iind in range(start_row, end_row):

            trj = self.fix_blankdump(iind)

            fraction_tmp[iind, :] = self._compute_localfraction(trj, types1, types2)

        del trj    

        if rank == 0:

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1: 
                    end_row = len(self.fnames)

                recv = np.empty([len(self.fnames), self.n_bins - 1])
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                fraction_tmp[start_row:end_row] = recv[start_row:end_row]
        else:
            send = fraction_tmp
            req = self.comm.Isend(send, dest=0)
            req.Wait()

        if rank == 0:

            result = {}

            fraction = Data.Data_TimeSeries2D(fraction_tmp, self.Nrepeat)

            fraction_final = fraction.mean
            fraction_final = fraction_final.transpose()

            result.update({'mean': fraction.mean.transpose()})
            result.update({'std': fraction.std.transpose()})

            res = Data.Data2D(z=fraction_final)
            xticks = np.linspace(0,len(self.source_dirs)*len(self.freqs),4)
            xticklabels = [format(i*self.Nfreq*0.006*pow(10,-6), '.4f') for i in xticks]
            plot_args = {
                'zmin': 0,
                'zmax': 1.0,
                'cmap': Colormaps.cmaps['G2B'],
                'cbarticks': [0, 1.0],
                'cbarlabel': r'$f_{\mathrm{C}}$',
                'xticks': xticks,
                'xticklabels': xticklabels,
                'xlabel': r'time ($\times 10^{6} \tau$)',
                'yticks': np.linspace(0, n_bins-1, 5),
                'yticklabels': np.linspace(self.bins[0], self.bins[-1], 5),
                'ylabel': r'$z$ ($\sigma$)',
            }

            res.plot(save='fraction_evol.png', show=False, plot_args=plot_args)

            self.results.update({'fraction': result})


    def _compute_localfraction(self, trj, types1, types2):

        complement = list(set(types1) & set(types2))
        if len(complement) == 0: raise Exception("set types right")

        for ind, tp in enumerate(complement):
            if ind == 0:
                logic_complement = trj[:,2] == tp
            else:
                logic_complement = np.logical_or(logic_complement, trj[:,2] == tp)

        for ind, tp in enumerate(types2):
            if ind == 0:
                logic = trj[:,2] == tp
            else:
                logic = np.logical_or(logic, trj[:,2] == tp)

        trj_complement = trj[logic_complement]
        trj_filter = trj[logic]

        hist_complement, bin_edges = np.histogram(trj_complement[:,5], bins=self.bins)
        hist, bin_edges = np.histogram(trj_filter[:,5], bins=self.bins)

        fraction = hist_complement / hist

        return fraction


    def compute_orientation(self, bond_types, zlo, zhi, width, buff, director):
        """
        compute orientation of morphology (Hermanns order parameter, S at every height) w.r.t director

        Arguments:
        ----------
        bond_types (list): types of bonds that connect A and B blocks (AB)
        zlo: lower limit
        zhi: uppder limit
        width: binning size
        director (list): director vertor

        Returns:
        --------
        N/A
        """

        self.zlo = zlo
        self.zhi = zhi
        self.width = width
        self.buff = buff 
        self.n_bins = int((zhi - zlo)/width)
        self.director = np.asarray(director)

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        if rank == 0:
            # number of bonds varies
            path = os.path.join(self.path, 'equil_0')
            files = [f for f in os.listdir(path) if 'data' in f]
            path_to_file = os.path.join(path, files[0])

            print('* Chain geometry info extracted')
            print('* path_to_file: {}'.format(path_to_file))

            f = open(path_to_file, 'r') 

            ind = 0
            flag = 0

            while flag == 0:
                line = f.readline()
                match = re.search('^\d+ bonds', line)
                if match:
                    n_bonds = int(line.split()[0])

                match = re.search('^Atoms', line)
                if match:
                    flag = 1

                ind += 1

            f.close()
            print('* # of bonds = {}'.format(n_bonds))

            subprocess.Popen('grep -A {} "Bonds" {} > temp.txt'.format(n_bonds + 1, path_to_file), shell=True).wait()
            bonds = np.loadtxt('temp.txt', skiprows=2)
            os.remove('temp.txt')
            #bonds = np.loadtxt(os.path.join(self.path, 'equil_0', 'pseudo_bonds.txt'))
            # choose bridge bonds between blocks
            for ind, bond_type in enumerate(bond_types):
                if ind == 0:
                    logic = bonds[:,1] == bond_type
                else:
                    logic = np.logical_or(logic, bonds[:,1] == bond_type)
            bonds = bonds[logic].astype(int)
            self.bonds = bonds

        else:
            bonds = None

        bonds = self.comm.bcast(bonds, root=0)
        self.bonds = bonds

        avg_rows_per_process = int(len(self.fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1:
            end_row = len(self.fnames)

        # 0: time, 1: bin along z
        S_tmp = np.empty([len(self.fnames), self.n_bins])
        for iind in range(start_row, end_row):

            trj = self.fix_blankdump(iind)

            trj = trj[np.argsort(trj[:,0])]

            cos_tmp = self._compute_cos(trj, iind)
            tmp = self._compute_S(cos_tmp)
            S_tmp[iind, :] = tmp

        if rank == 0:
            S_1D = np.empty([len(self.fnames), self.n_bins])
            S_1D[start_row:end_row, :] = S_tmp[start_row:end_row]

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1:
                    end_row = len(self.fnames)

                recv = np.empty([len(self.fnames), self.n_bins])
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                S_1D[start_row:end_row, :] = recv[start_row:end_row]
        else:
            send = S_tmp
            req = self.comm.Isend(send, dest=0)
            req.Wait()
    
        if rank == 0:

            result = {}
            op = Data.Data_TimeSeries2D(S_1D, self.Nrepeat)

            S_final = op.mean
            S_final = S_final.transpose()

            print(np.max(S_final), np.min(S_final))

            result.update({'mean': op.mean})
            result.update({'std': op.std})

            xticks = np.linspace(0,len(self.source_dirs)*len(self.freqs),4)

            from mpl_toolkits.axes_grid1 import make_axes_locatable
            fig, ax = plt.subplots(figsize=(4,3))
            im = ax.imshow(S_final, vmin=-0.5, vmax=1.0, cmap=Colormaps.cmaps['BWR'], origin='lower')
            ax.plot([xticks[0], xticks[-1]], [self.h, self.h], 'k--', lw=0.5)
            ax.set_aspect('auto')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_ticks([-0.5, 0, 1.0])
            cbar.set_label(r'$S$')
            ax.set_xlim(xticks[0] - 0.5, xticks[-1] + 0.5)
            ax.set_xticks(xticks)
            ax.set_xticklabels([format(i*self.Nfreq*0.006*pow(10,-6), '.4f') for i in xticks])
            ax.set_xlabel(r'time ($\times 10^{6} \tau$)')
            ax.set_yticks(np.linspace(0, self.n_bins-1,5))
            ax.set_ylim(0 - 0.5, self.n_bins-1 + 0.5)
            ax.set_yticklabels(np.linspace(self.zlo, self.zhi, 5))
            ax.set_ylabel(r'$z$ ($\sigma$)')
            plt.tight_layout(pad=1,h_pad=None,w_pad=None,rect=None)
            plt.savefig('angle_map.png', dpi=1000)

            self.results.update({'angle_map': result})


    def compute_chainalignment(self, bond_types, zlo, zhi, width, buff, director):
        """
        compute angles of chains at every height w.r.t director

        Arguments:
        ----------
        bond_types (list): types of bonds that connect A and B blocks (AB)
        zlo: lower limit
        zhi: uppder limit
        width: binning size
        director (list): director vertor

        Returns:
        --------
        N/A
        """

        self.zlo = zlo
        self.zhi = zhi
        self.width = width
        self.buff = buff
        self.n_bins = int((zhi - zlo)/width)
        self.director = np.asarray(director)

        bins = np.linspace(0, 180, 72+1) # delta ranging from 0 to 180, delta_phi = 2.5
        #bins = np.linspace(0, 360, 144+1) # dleta ranging from 0 to 360, delta_phi = 2.5

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        if rank == 0:
            if zhi == None:
                h = self.h
            else:
                h = zhi
        else:
            h = None

        h = self.comm.bcast(h, root=0)

        ll_max = zhi - width # lower limit max

        if rank == 0:
            # number of bonds varies
            path = os.path.join(self.path, 'equil_0')
            files = [f for f in os.listdir(path) if 'data' in f]
            path_to_file = os.path.join(path, files[0])

            print('* Chain geometry info extracted')
            print('* path_to_file: {}'.format(path_to_file))

            f = open(path_to_file, 'r')

            ind = 0
            flag = 0

            while flag == 0:
                line = f.readline()
                match = re.search('^\d+ bonds', line)
                if match:
                    n_bonds = int(line.split()[0])

                match = re.search('^Atoms', line)
                if match:
                    flag = 1

                ind += 1

            f.close()
            print('* # of bonds = {}'.format(n_bonds))

            subprocess.Popen('grep -A {} "Bonds" {} > temp.txt'.format(n_bonds + 1, path_to_file), shell=True).wait()
            bonds = np.loadtxt('temp.txt', skiprows=2)
            os.remove('temp.txt')
            #bonds = np.loadtxt(os.path.join(self.path, 'equil_0', 'pseudo_bonds.txt'))
            # choose bridge bonds between blocks
            for ind, bond_type in enumerate(bond_types):
                if ind == 0:
                    logic = bonds[:,1] == bond_type
                else:
                    logic = np.logical_or(logic, bonds[:,1] == bond_type)
            bonds = bonds[logic].astype(int)
            self.bonds = bonds

        else:
            bonds = None

        bonds = self.comm.bcast(bonds, root=0)
        self.bonds = bonds

        avg_rows_per_process = int(len(self.fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1:
            end_row = len(self.fnames)

        # column 0: time, 1: height, 2: angle_bin
        angle_hist_tmp = np.empty((len(self.fnames), self.n_bins*3, len(bins) - 1))
        for iind in range(start_row, end_row):

            trj = np.loadtxt(self.fnames[iind], skiprows=9)
            trj = trj[np.argsort(trj[:,0])]

            cos_tmp = self._compute_cos(trj, iind)

            for ind, ll in enumerate(np.linspace(0, self.zhi - self.width, self.n_bins)):
                ul = ll + self.width

                cos_filtered = cos_tmp[np.logical_and(cos_tmp[:,0] > ll - self.buff, cos_tmp[:,0] < ul + self.buff)]

                for bond_ind, bond_type in enumerate(bond_types, 1):
                    angles = np.degrees(np.arccos(cos_filtered[cos_filtered[:,2] == bond_type, 1]))

                    hist, bin_edges = np.histogram(angles, bins)
                    angle_hist_tmp[iind, ind + bond_ind*self.n_bins, :] = hist

                angles = np.degrees(np.arccos(cos_filtered[:,1]))

                hist, bin_edges = np.histogram(angles, bins)
                angle_hist_tmp[iind, ind, :] = hist

            bin_edges = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2

        del trj

        if rank == 0:
            angle_hist = np.empty((len(self.fnames), self.n_bins*3, len(bin_edges)))
            angle_hist[start_row:end_row, :, :] = angle_hist_tmp[start_row:end_row, :, :]

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1:
                    end_row = len(self.fnames)

                recv = np.empty((len(self.fnames), self.n_bins*3, len(bin_edges)))
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                angle_hist[start_row:end_row, :, :] = recv[start_row:end_row, :, :]
        else:
            send = angle_hist_tmp
            req = self.comm.Isend(send, dest=0)
            req.Wait()

        if rank == 0:

            result = {}
            result.update({'z': np.linspace(0, self.zhi - self.width, self.n_bins) + self.width/2})
            result.update({'phi': bin_edges})

            align = Data.Data_TimeSeries3D(angle_hist, self.Nrepeat)

            result.update({'hist_mean': align.mean})
            result.update({'hist_std': align.std})

            self.results.update({'chain_alignment': result})


    def _compute_cos(self, trj, iind):
        '''
        compute the position of the mid point of A-b-B chain
        and the angle made by the AB vector and the director
        '''

        ind_A = self.bonds[:,2] - 1
        ind_B = self.bonds[:,3] - 1

        #####
        a0 = trj[ind_A - 4]
        next_pt = trj[ind_A - 3]
        a1 = a0.copy()
        a1[:,3:] = self.d_pbc(next_pt, trj[ind_A - 4]) + trj[ind_A - 4,3:]
        next_pt = trj[ind_A - 2]
        a2 = a0.copy()
        a2[:,3:] = self.d_pbc(next_pt, trj[ind_A - 3]) + trj[ind_A - 3,3:]
        next_pt = trj[ind_A - 1]
        a3 = a0.copy()
        a3[:,3:] = self.d_pbc(next_pt, trj[ind_A - 2]) + trj[ind_A - 2,3:]
        next_pt = trj[ind_A - 0]
        a4 = a0.copy()
        a4[:,3:] = self.d_pbc(next_pt, trj[ind_A - 1]) + trj[ind_A - 1,3:]

        A = (a0 + a1 + a2 + a3 + a4)/5

        next_pt = trj[ind_B + 1]
        b0 = a0.copy()
        b0[:,3:] = self.d_pbc(next_pt, trj[ind_B - 0]) + trj[ind_B - 0,3:]
        next_pt = trj[ind_B + 2]
        b1 = a0.copy()
        b1[:,3:] = self.d_pbc(next_pt, trj[ind_B + 1]) + trj[ind_B + 1,3:]
        next_pt = trj[ind_B + 3]
        b2 = a0.copy()
        b2[:,3:] = self.d_pbc(next_pt, trj[ind_B + 2]) + trj[ind_B + 2,3:]
        next_pt = trj[ind_B + 4]
        b3 = a0.copy()
        b3[:,3:] = self.d_pbc(next_pt, trj[ind_B + 3]) + trj[ind_B + 3,3:]
        next_pt = trj[ind_B + 5]
        b4 = a0.copy()
        b4[:,3:] = self.d_pbc(next_pt, trj[ind_B + 4]) + trj[ind_B + 4,3:]

        B = (b0 + b1 + b2 + b3 + b4)/5

        #A = trj[ind_A - 3]
        #B = trj[ind_B + 3]

        d = self.d_pbc(A, B)
        z_avg = (trj[ind_A,5] + trj[ind_B,5])/2

        cos = np.empty((len(self.bonds), 4))
        #cos = np.empty((A.shape[0], 2)) #

        cos[:,0] = z_avg
        cos[:,1] = np.dot(d, self.director) / np.linalg.norm(d, axis=1)
        cos[:,2] = self.bonds[:,1]
        cos[:,3] = np.zeros(len(cos))
        cos[d[:,1] < 0, 3] = 1

        return cos
    

    def _compute_S(self, cos):
        '''
        compute the Hermanns order parameter
        '''

        # DEFAULT: width = 1.0, buffer = 0.5 (both above and below)
        bin_starts = np.linspace(self.zlo, self.zhi - self.width, int((self.zhi-self.zlo)/self.width))
        bin_ends = bin_starts + self.width
        cos_binned = [cos[np.logical_and(cos[:,0] > bin_starts[i] - self.buff, cos[:,0] < bin_ends[i] + self.buff), 1] for i in range(len(bin_starts))]

        S = np.empty(len(bin_starts))

        for bin_ind in range(len(cos_binned)):
            if len(cos_binned[bin_ind]) == 0:
                S[bin_ind] = np.nan
            else:
                S[bin_ind] = (3*np.mean(cos_binned[bin_ind]**2) - 1)/2

        return S


    def d_pbc(self, vector1, vector2):
        """
        reconstruct beads considering pbc
        """
        boxlength = [self.lx, self.ly]

        l_ref = [self.lx/2.0, self.ly/2.0]
        vector = vector1[:,3:] - vector2[:,3:]

        for i in [0, 1]:
            vector[vector[:,i] < (-1)*l_ref[i], i] += boxlength[i]
            vector[vector[:,i] > l_ref[i], i] -= boxlength[i]

        return vector


    def compute_lambda(self, types):
        """
        compute stratification of beads of given types comparing the top and bottom halves

        Arguments:
        ----------
        types (list): types of beads of interesting

        Returns:
        --------
        N/A
        """

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        avg_rows_per_process = int(len(self.fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1: 
            end_row = len(self.fnames)

        if rank == 0:
            h = self.h
        else:
            h = None 

        h = self.comm.bcast(h, root=0)

        lambda_tmp = np.empty(len(self.fnames))
        for iind in range(start_row, end_row):

            trj = self.fix_blankdump(iind)

            for tp_ind, tp in enumerate(types):
                if tp_ind == 0:
                    logic = trj[:,2] == tp
                else:
                    logic = np.logical_or(logic, trj[:,2] == tp)

            logic_top = trj[:,5] > h/2

            numerator = len(trj[np.logical_and(logic, logic_top)]) - len(trj[np.logical_and(logic, ~logic_top)])
            denominator = len(trj[logic])

            lambda_tmp[iind] = numerator/denominator

        if rank == 0:
            lambda_res = np.empty(len(self.fnames))
            lambda_res[start_row:end_row] = lambda_tmp[start_row:end_row]

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1:
                    end_row = len(self.fnames)

                recv = np.empty(len(self.fnames))
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                lambda_res[start_row:end_row] = recv[start_row:end_row]
        else:
            send = lambda_tmp
            req = self.comm.Isend(send, dest=0)
            req.Wait()

        if rank == 0:
            self.results.update({'lambda': lambda_res})


    def compute_objsize(self, types, delta=[0.75, 0.75, 0.5]):
        '''
        carry out film resconstruction by looking at every height
        with the matrix flood-filled, calculate the size of each morphological object
        '''

        from scipy.ndimage import gaussian_filter

        sys.setrecursionlimit(100000)

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()

        fnames = self.fnames.copy()
        fnames.pop(0) #not to compute for the initial configuration

        fnames = [fnames[self.Nrepeat*ind: self.Nrepeat*(ind + 1)] for ind in range(int(len(fnames)/self.Nrepeat))]
     
        avg_rows_per_process = int(len(fnames)/size)

        start_row = rank * avg_rows_per_process
        end_row = start_row + avg_rows_per_process
        if rank == size - 1:
            end_row = len(fnames)

        if rank == 0:
            h = self.h
        else:
            h = None

        h = self.comm.bcast(h, root=0)

        delta_h = 2
        ll_max = int(h - delta_h)

        res = np.empty((len(fnames), ll_max + 1, 3))
        for row_ind in range(start_row, end_row):

            for f_ind, fname in enumerate(fnames[row_ind]):
                trj, box = loadtrj(fname)

                if f_ind == 0:
                    binary3D = binarize_3D(trj, box, delta, types)
                    binary3D[binary3D > 0] = 1
                else:
                    new_binary3D = binarize_3D(trj, box, delta, types)
                    new_binary3D[new_binary3D > 0] = 1
                    dim_z = min(binary3D.shape[2], new_binary3D.shape[2])
                    binary3D = binary3D[:, :, :dim_z] + new_binary3D[:, :, :dim_z]

            del trj

            for ind, ll in enumerate(np.linspace(0, ll_max, ll_max + 1, dtype=int)):
                ul = ll + delta_h

                binary2D = binary3D[:, :, ll*2: ul*2].sum(axis=2)

                binary2D[binary2D > 0] = 1
                binary2D = gaussian_filter(binary2D, 1)
                binary2D[binary2D > 0.5] = 1
                binary2D[binary2D < 0.5] = 0

                'film reconstruction snapshots'
                #plt.figure()
                #plt.imshow(binary2D, origin='lower')
                #plt.savefig('test_{:02d}_{:02d}.png'.format(row_ind, ind))
                #plt.close()

                object_size = []
                line = []
                dot = []
                try:
                    blobs = group2blob(binary2D)

                    h, w = binary2D.shape

                    for blob_ind, blob in enumerate(blobs):
                        #print (len(blob))
                        if len(blob) < 10:
                            for j in range(len(blob)):
                                binary2D[blob[j][0]%h, blob[j][1]%w]=0
                        elif (len(blob) > 37.5) and (len(blob) < 300):
                            dot.append(len(blob))
                        elif (len(blob) > 300):
                            line.append(len(blob))
                        object_size.append(len(blob))

                except:
                    object_size.append(0)

                res[row_ind, ind, 0] = (ll + ul)/2

                if (len(line) + len(dot)) > 0:
                    ## simply counting the relative number of dots
                    res[row_ind, ind, 1] = len(dot) / (len(line) + len(dot))

                    ## comparing the areas covered by dots and lines
                    res[row_ind, ind, 2] = sum(dot) / sum(object_size)
                else:
                    res[row_ind, ind, 1] = -1

                    res[row_ind, ind, 2] = -1

        if rank == 0:

            for iind in range(1, size):
                start_row = iind*avg_rows_per_process
                end_row = start_row + avg_rows_per_process
                if iind == size - 1:
                    end_row = len(fnames)

                recv = np.empty((len(fnames), ll_max + 1 , 3))
                req = self.comm.Irecv(recv, source=iind)
                req.Wait()

                res[start_row:end_row] = recv[start_row:end_row]
        else:
            send = res
            req = self.comm.Isend(send, dest=0)
            req.Wait()

        if rank == 0:

            self.results.update({'f_dot': res})
            #for ind, val in enumerate(res): 
            #    plt.figure()
            #    plt.plot(val[:,0], val[:,1])
            #    plt.plot(val[:,0], val[:,2])
            #    plt.ylim(0, 1)
            #    plt.savefig('test2_{:02d}.png'.format(ind))
            #    plt.close()


class Thermo:
    """
    Attributes:
        log_files (list): log files
    """

    def __init__(self, log_files=None):
        """
        Thermo constructor
        """

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


def sigmoid(v, x):
    return v['c']/(1.0 + np.exp(-v['a']*(x - v['b'])))


def loadtrj(f_name, skiprows=9):

    f = open(f_name, 'r')
    box = []
    for i in range(5):
        f.readline()
    for i in range(3):
        line = f.readline().split()
        box.append([float(line[0]), float(line[1])])
    f.close()

    trj = np.loadtxt(f_name, skiprows=skiprows)
    trj = trj[np.argsort(trj[:,0])]
    trj[:,5] = trj[:,5] - box[2][0]

    return trj, box


def binarize_3D(trj, box, delta, types):
    """
    Arguments:
    ----------
    trj: MD trajectories, N by 6 2D array (atom-id, mol-id, atom-type, x, y, z)
    box: lower and upper limits of simulation box in each axis
        box[first_index][second_index]
        first_index = 0, 1, 2 for x, y, z
        second_index = 0, 1 for lower, upper limits
    delta: bin width 
        delta[index]
        index = 0, 1, 2, for x, y, z
    types: types used to construct blobs

    Returns:
    --------
    binary: A 3D array. Each element corresponds to number of atoms for a corresdponding bin.

    """

    binx = np.linspace(box[0][0], box[0][1], int((box[0][1]-box[0][0])/delta[0])+1)
    biny = np.linspace(box[1][0], box[1][1], int((box[1][1]-box[1][0])/delta[1])+1)
    binz = np.linspace(0, math.ceil(max(trj[:,5])), int(math.ceil(max(trj[:,5]))/delta[2])+1)

    for ind, tp in enumerate(types):
        if ind == 0:
            logic = trj[:,2] == tp
        else:
            logic = np.logical_or(logic, trj[:,2] == tp)
    
    trj_filter = trj[logic]

    binary = stats.binned_statistic_dd(trj_filter[:,[4, 3, 5]], None, statistic='count', bins=[biny, binx, binz]).statistic

    return binary


def binarize_2D(trj, box, delta, types):
    """
    Discretize simulation box into voxels and identifies a voxel within which each atom falls

    Arguments:
    ----------
    trj: MD trajectories given as a 2D array (N by 6)
        atom-id, mol-id, atom-type, x, y, z
    box: lower and upper limits of simulation box in each axis
        box[first_index][second_index]
        first_index = 0, 1, 2 for x, y, z
        second_index = 0, 1 for lower, upper limits
    delta: bin width 
        delta[index]
        index = 0, 1, 2, for x, y, z
    types: types used to construct blobs

    Returns:
    --------
    binary: A 2D array. Each element corresponds to number of atoms for a corresdponding bin.

    """
    
    binx = np.linspace(box[0][0], box[0][1], int((box[0][1]-box[0][0])/delta[0])+1)
    biny = np.linspace(box[1][0], box[1][1], int((box[1][1]-box[1][0])/delta[1])+1)

    for ind, tp in enumerate(types):
        if ind == 0:
            logic = trj[:,2] == tp
        else:
            logic = np.logical_or(logic, trj[:,2] == tp)
    
    trj_filter = trj[logic]

    binary = stats.binned_statistic_dd(trj_filter[:,[4, 3]], None, statistic='count', bins=[biny, binx]).statistic

    return binary
   

def view_xy(binary_3D, delta, method='topview'):
    """
    Generate an orthographic, top-view image by determining the highest material point along z axis at each grid point on the xy plane

    Arguments:
    ----------
    binary_3D: A 3D array (binx by biny by binz) whose elements are atom counts
    delta: bin width 
        delta[index]
        index = 0, 1, 2, for x, y, z

    Returns:
    --------
    topview: A 2D array (binx by biny). Simply put, binary_3D projected on xy plane

    """

    array = binary_3D

    if method == 'topview':
        'weight voxel - sem like topview'
        topview = np.zeros([array.shape[0], array.shape[1]])
        for iind in range(topview.shape[0]):
            for jind in range(topview.shape[1]):
                kind = array.shape[2]
                flag = 0
                while (flag == 0):
                    kind = kind - 1
                    if array[iind,jind,kind] != 0:
                        flag = 1
                    if kind == 0:
                        flag = 1
                topview[iind,jind] = array[iind,jind,kind]*(kind*delta[2])**1 #power of 1: linear & 2: quadratic

    elif method == 'xray_sum':
        'sum up over thickness - xray'
        topview = np.sum(binary_3D, axis=2)

    elif method == 'xray_average':
        'averagve over thickness - xray'
        topview = np.average(binary_3D, axis=2)

    return topview        


def floodfill_pbc(binary_2D, x, y, blob, repeat=4):
    """
    Floodfill over periodic boundary conditions

    Arguments:
    ----------
    binary_2D: A 2D array (binx by biny) whose elements represent atoms counts or weighted intensities
    x, y: indeces for binary_2D[x][y]
    blob: sets of [i, j] for each object determined by floodfill
    repeat: floodfill is applied over this many periods

    Returns:
    --------

    """
    array = binary_2D

    if array[x%len(array), y%len(array[x%len(array)])] == 1:
        array[x%len(array), y%len(array[x%len(array)])] = 0.5
        blob.append([x,y])

        if x > len(array)*(-repeat):
            floodfill_pbc(array, x-1, y, blob)
        if x < repeat*len(array)-1:
            floodfill_pbc(array, x+1, y, blob)
        if y > len(array[x%len(array)])*(-repeat):
            floodfill_pbc(array, x, y-1, blob)
        if y < repeat*len(array[x%len(array)])-1:
            floodfill_pbc(array, x, y+1, blob)


def group2blob(binary2D):

    h, w = binary2D.shape

    blobs = []
    flag=0
    while (flag == 0):
        candidate = []
        blob_tmp = []
        for i in range(h):
            for j in range(w):
                if binary2D[i,j] == 1:
                    candidate.append([i,j])
        if len(candidate) > 0:
            floodfill_pbc(binary2D, candidate[0][0], candidate[0][1], blob_tmp)
            blobs.append(blob_tmp)
        else:
            flag=1

    return blobs, blob_tmp


def trj2bin(coords, box, delta):
    ''' 
    maps a point in coordinate system to a point in binary image
    '''
    bin_y = coords[4]%(box[1][1]-box[1][0])//delta[1]
    bin_x = coords[3]%(box[0][1]-box[0][0])//delta[0]

    return [bin_y, bin_x]


def d_pbc2D(vector1, vector2, boxlength):
    l_ref = [boxlength[0]/2.0, boxlength[1]/2.0]
    vector = vector1 - vector2
    
    if vector[0] < (-1)*l_ref[0]:
        vector[0] = vector[0] + boxlength[0]
    elif vector[0] > l_ref[0]:
        vector[0] = vector[0] - boxlength[0]
    else:
        vector[0] = vector[0]

    if vector[1] < (-1)*l_ref[1]:
        vector[1] = vector[1] + boxlength[1]
    elif vector[1] > l_ref[1]:
        vector[1] = vector[1] - boxlength[1]
    else:
        vector[1] = vector[1]

    return vector


def contour_finder(matrix, x, y, contour, repeat=4):
    'floodfills over periodic boundary conditions'

    if matrix[x%len(matrix), y%len(matrix[x%len(matrix)])] == 0.5:
        matrix[x%len(matrix), y%len(matrix[x%len(matrix)])] = 0.4

        if x > len(matrix)*(-repeat):
            if matrix[(x-1)%len(matrix), y%len(matrix[x%len(matrix)])] == 0:
                contour.append([x,y])
            contour_finder(matrix, x-1, y, contour)

        if x < repeat*len(matrix) - 1:
            if matrix[(x+1)%len(matrix), y%len(matrix[x%len(matrix)])] == 0:
                contour.append([x,y])
            contour_finder(matrix, x+1, y, contour)
             
        if y > len(matrix[x%len(matrix)])*(-repeat):
            if matrix[x%len(matrix), (y-1)%len(matrix[x%len(matrix)])] == 0:
                contour.append([x,y])
            contour_finder(matrix, x, y-1 ,contour)

        if y < repeat*len(matrix[x%len(matrix)]) - 1:
            if matrix[x%len(matrix), (y+1)%len(matrix[x%len(matrix)])] == 0:
                contour.append([x,y])
            contour_finder(matrix, x, y+1, contour)


