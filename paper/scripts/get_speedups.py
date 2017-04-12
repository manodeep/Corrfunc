#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import cPickle as pickle
except ImportError:
    import pickle
        
from Corrfunc.io import read_catalog
import multiprocessing
max_threads = multiprocessing.cpu_count()

    
def read_file(filename):
    """
    Reads in the file I created manually (by recompiling and adding timers)
    
    Not used any more but left for historical reasons (the first 'speedup'
    plots were generated with this function)
    """
    dtype = np.dtype([('same_cell', np.int32),
                      ('N1', np.int),
                      ('N2', np.int),
                      ('time', np.float)
                      ])
    if pd is not None:
        timings = pd.read_csv(filename, header=None,
                              engine="c",
                              dtype={'same_cell': np.int32,
                                     'N1': np.int,
                                     'N2': np.int,
                                     'time': np.float},
                              index_col=None,
                              names=['same_cell', 'N1', 'N2', 'time'],
                              delim_whitespace=True)
    else:
        timings = np.loadtxt(filename, dtype=dtype)
    return timings


class nf(float):
    def __repr__(self):
        str_repr = '%.1f' % (self.__float__(),)
        if str_repr[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


def run_wp(boxsize, x, y, z, pimax, nthreads=max_threads, isa=None):
    import Corrfunc
    from Corrfunc.theory import wp
    from os.path import dirname, abspath, join as pjoin
    #binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    #                "../theory/tests/", "bins")
    binfile = './bins'
    _, cell_time = wp(boxsize, pimax, nthreads, binfile,
                      x, y, z, c_cell_timer=True, isa=isa,
                      verbose=True)
    
    return cell_time

        
def main():
    import sys
    if len(sys.argv) == 1:
        print("Running cell timers for wp")
        all_isa = ['avx', 'sse42', 'fallback']
        #x, y, z = read_catalog()
        boxsize = 1100.0
        pimax = 45.0
        
        points = np.loadtxt('halos_emulator_1100box_Neff3_00.txt')
        numpart = int(1.*len(points))
        assert (points >= 0).all() and (points < 1100.).all()
        dtype = points.dtype  # float64
        points = points.reshape(-1).view(dtype=[('x',dtype,3)])
        subsample = np.random.choice(points, numpart, replace=False)
        subsample = subsample.view(dtype=dtype).reshape(-1,3)
        x, y, z = subsample.T
        
        cell_timings = dict()
        serial_timings = dict()
        for isa in all_isa:
            cell_timings[isa] = dict()
            serial_timings[isa] = dict()

            # First run the serial (single threaded)
            timings = run_wp(boxsize, x, y, z, pimax,
                             nthreads=1, isa=isa)
            (serial_timings[isa])[1] = timings

            # then the one with all threads
            # max_threads is not required but being explicit
            timings = run_wp(boxsize, x, y, z, pimax,
                             nthreads=max_threads, isa=isa)
            (cell_timings[isa])[max_threads] = timings

        with open('wp_cell_timers.pkl', 'wb') as outfile:
            pickle.dump([all_isa, cell_timings, serial_timings], outfile,
                        protocol=pickle.HIGHEST_PROTOCOL)

    else:
        
        timings_file = sys.argv[1]
        print("Loading benchmarks from file = {0}".format(timings_file))
        with open(timings_file, 'rb') as pkl_file:
            all_isa, _, serial_timings = pickle.load(pkl_file)
            
        legend = ['AVX', 'SSE4.2', 'Fallback']
        base_string = 'wp'

        all_speedup = []
        base_timing = (serial_timings['fallback'])[1]['time_in_ns']
        N1_parts = (serial_timings['fallback'])[1]['N1']
        N2_parts = (serial_timings['fallback'])[1]['N2']
        gridsize = 40
        cb_range = [0.0, 3.0e-5]
        contour_nlevels = 4
        xlimits = [0, 40]
        ylimits = xlimits
        xlabel = 'Number of points in a cell'
        ylabel = xlabel

        
        cb_diff = (cb_range[1] - cb_range[0])
        '''
        positive_Ncolors = int((cb_range[1] - 1.0) / cb_diff * 256)
        negative_Ncolors = 256 - positive_Ncolors
        colors1 = cm.OrRd(np.linspace(0.0, 1.0, negative_Ncolors))
        colors2 = cm.viridis(np.linspace(0.0, 1.0, positive_Ncolors))
        # combine them and build a new colormap
        colors = np.vstack((colors1, colors2))
        mycmap = mcolors.LinearSegmentedColormap.from_list('my_colormap',
                                                           colors)
        '''
        mycmap = 'viridis'
        matplotlib.style.use('default')
        # Label levels with specially formatted floats
        if plt.rcParams["text.usetex"]:
            cntr_fmt = r'%r\%%'
        else:
            cntr_fmt = '%r%%'

        # Want fallback to appear first
        all_isa.reverse()
        legend.reverse()

        for ii, isa in enumerate(all_isa):
            if ii == 0:
                continue
            
            this_timing = (serial_timings[isa])[1]['time_in_ns']
            ind = (np.where((this_timing > 0.0) & (base_timing > 0.0)))[0]
            speedup = base_timing[ind] / this_timing[ind]
            all_speedup.append(speedup)
            print("Min speedup = {0}. Max = {1}".format(
                min(speedup), max(speedup)))
            bad = (np.where(speedup <= 1.0))[0]
            bad_timings_base = np.sum(base_timing[ind[bad]])
            bad_timings = np.sum(this_timing[ind[bad]])
            print("Cells with slowdown  {3}({4:4.3f}%): Base takes - {0:8.3f} "
                  "sec while {1} takes {2:8.3f} seconds".format(
                      bad_timings_base/1e9,
                      legend[ii],
                      bad_timings/1e9,
                      len(bad),
                      100.0 * len(bad) / len(ind)))

            good = (np.where(speedup > 1.0))[0]
            good_timings_base = np.sum(base_timing[ind[good]])
            good_timings = np.sum(this_timing[ind[good]])
            print("Cells with speedup {3}({4:4.3f}%): Base takes - {0:8.3f} sec "
                  "while {1} takes {2:8.3f} seconds".format(
                      good_timings_base/1e9,
                      legend[ii],
                      good_timings/1e9,
                      len(good),
                      100.0 * len(good) / len(ind)))

            fig = plt.figure(1, figsize=(8, 8))
            figsize = 0.6
            left = 0.1
            bottom = 0.1
            top_aspect = 0.15
            hist_area = [left, bottom + figsize, figsize, figsize * top_aspect]
            axhist = plt.axes(hist_area)
            axhist.autoscale(enable=True, axis="y")
            axhist.set_xlim(xlimits)
            plt.setp(axhist.get_xticklabels(), visible=False)
            axhist.axis('off')
            axhist.hist(N1_parts[ind], gridsize, range=xlimits,
                        color='0.5')

            hist_time_area = [left + figsize, bottom,
                              figsize*top_aspect, figsize]
            ax_time = plt.axes(hist_time_area)
            ax_time.autoscale(enable=True, axis="x")
            ax_time.set_ylim(ylimits)
            plt.setp(ax_time.get_yticklabels(), visible=False)
            plt.setp(ax_time.get_xticklabels(), visible=False)
            ax_time.axis('off')
            ax_time.hist(N1_parts[ind], gridsize, weights=this_timing[ind],
                         range=xlimits, orientation="horizontal",
                         color='0.5')

            im_area = [left, bottom, figsize, figsize]
            ax = plt.axes(im_area)
            ax.set_autoscale_on(False)
            ax.set_xlim(xlimits)
            ax.set_ylim(ylimits)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            xedges = np.linspace(xlimits[0], xlimits[1], gridsize)
            yedges = np.linspace(ylimits[0], ylimits[1], gridsize)
            cell_time, xedges, yedges = np.histogram2d(
                N1_parts, N2_parts, (xedges, yedges),
                weights=base_timing, normed=False)

            cell_time /= np.sum(cell_time)
            cell_time *= 100.0
            cell_time_1d = cell_time.flatten()
            sorted_ind = np.argsort(cell_time_1d)
            cum_sorted_time = np.cumsum(cell_time_1d[sorted_ind])
            correct_order_cum_time = np.empty_like(cum_sorted_time)
            for kk, ct in zip(sorted_ind, cum_sorted_time):
                correct_order_cum_time[kk] = ct

            correct_order_cum_time = correct_order_cum_time.reshape(
                cell_time.shape)
            extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
            xarr, yarr = np.meshgrid(xedges[0:-1], yedges[0:-1])
            contours = ax.contour(xarr, yarr,
                                  correct_order_cum_time, contour_nlevels,
                                  linewidths=3.0,
                                  extent=extent,
                                  cmap=cm.Greys)

            # Recast levels to new class
            # Reverse the levels to show that the contours represent
            # enclosed fraction of time spent
            contours.levels = [nf(val) for val in contours.levels[::-1]]
            ax.clabel(contours, contours.levels, fmt=cntr_fmt,
                      inline=True, fontsize=10)

            # Now plot the image for the speedup
            normalized_this_timing = this_timing/this_timing.sum()
            im = ax.hexbin(N1_parts[ind], N2_parts[ind],
                           #C=speedup[ind],
                           C=normalized_this_timing[ind],
                           vmin=cb_range[0], vmax=cb_range[1],
                           cmap=mycmap, gridsize=gridsize)
            plt.figtext(left + figsize - 0.03, bottom + figsize - 0.05,
                        '{0}'.format(legend[ii]), fontsize=16, ha='right')
            cbar_offset = 0.08
            cbar_width = 0.03
            cbar_ax = fig.add_axes([left + figsize + figsize*top_aspect +
                                    cbar_offset, bottom,
                                    cbar_width, figsize])
            cb = fig.colorbar(im, extend='both', format="%.1f",
                              ticks=np.linspace(cb_range[0], cb_range[1],
                                                cb_diff + 1.0),
                              cax=cbar_ax)
            cb.set_label('Speedup rel. to non-vectorized code')
            if 'laptop' in timings_file:
                exts = 'laptop_'
            elif 'stampede' in timings_file:
                exts = 'stampede_'
            elif 'bender' in timings_file:
                exts = 'bender_'
            else:
                exts = ''
            
            plt.savefig('{1}_{2}Speedup_{0}.png'.format(legend[ii],
                                                        base_string,
                                                        exts),
                        dpi=400)
            
            plt.savefig('{1}_{2}Speedup_{0}.pdf'.format(legend[ii],
                                                        base_string,
                                                        exts),
                        dpi=400)
            fig.clear()
            ax.clear()
            axhist.clear()
            ax_time.clear()
            plt.close(fig)

if __name__ == '__main__':
    main()
