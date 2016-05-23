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


def read_file(filename):
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
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


def main():
    base_dir = '../xi_theory/wp/'
    base_string = 'wp'
    files = ['timings_naive', 'timings_sse', 'timings_avx']
    files = [base_dir + f for f in files]
    legend = ['Naive', 'SSE4.2', 'AVX']
    numfiles = len(files)
    all_timings = []
    for filename in files:
        timings = read_file(filename)
        all_timings.append(timings)

    all_speedup = []
    base_timing = (all_timings[0])['time']
    N1_parts = (all_timings[0])['N1']
    N2_parts = (all_timings[0])['N2']
    gridsize = 40
    cb_range = [0.0, 5.0]
    contour_nlevels = 4
    xlimits = [0, 1000]
    ylimits = xlimits
    xlabel = 'Number of points in a cell'
    ylabel = xlabel

    cb_diff = (cb_range[1] - cb_range[0])
    positive_Ncolors = int((cb_range[1] - 1.0) / cb_diff * 256)
    negative_Ncolors = 256 - positive_Ncolors
    colors1 = cm.OrRd(np.linspace(0.0, 1.0, negative_Ncolors))
    colors2 = cm.viridis(np.linspace(0.0, 1.0, positive_Ncolors))
    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mycmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    matplotlib.style.use('default')
    # Label levels with specially formatted floats
    if plt.rcParams["text.usetex"]:
        cntr_fmt = r'%r\%%'
    else:
        cntr_fmt = '%r%%'

    for i in xrange(numfiles):
        if i == 0:
            continue
        this_timing = (all_timings[i])['time']
        ind = (np.where((this_timing > 0.0) & (base_timing > 0.0)))[0]
        speedup = base_timing[ind] / this_timing[ind]
        all_speedup.append(speedup)
        print("Min speedup = {0}. Max = {1}".format(
            min(speedup), max(speedup)))
        bad = (np.where(speedup <= 1.0))[0]
        bad_timings_base = np.sum(base_timing[ind[bad]])
        bad_timings = np.sum(this_timing[ind[bad]])
        print("Cells with slowdown  {3}({4:4.3f}%): Base takes - {0:8.3f} sec "
              "while {1} takes {2:8.3f} seconds".format(
                  bad_timings_base,
                  legend[i],
                  bad_timings,
                  len(bad),
                  100.0 * len(bad) / len(ind)))

        good = (np.where(speedup > 1.0))[0]
        good_timings_base = np.sum(base_timing[ind[good]])
        good_timings = np.sum(this_timing[ind[good]])
        print("Cells with speedup {3}({4:4.3f}%): Base takes - {0:8.3f} sec "
              "while {1} takes {2:8.3f} seconds".format(
                  good_timings_base,
                  legend[i],
                  good_timings,
                  len(good),
                  100.0 * len(good) / len(ind)))

        fig = plt.figure(1, figsize=(8, 8))
        figsize = 0.7
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
        im = ax.hexbin(N1_parts[ind], N2_parts[ind], C=speedup[ind],
                       vmin=cb_range[0], vmax=cb_range[1],
                       cmap=mycmap, gridsize=gridsize)
        plt.figtext(left + figsize - 0.03, bottom + figsize - 0.05,
                    '{0}'.format(legend[i]), fontsize=16, ha='right')
        cbar_offset = 0.05
        cbar_width = 0.03
        cbar_ax = fig.add_axes([left + figsize + cbar_offset, bottom,
                                cbar_width, figsize])
        cb = fig.colorbar(im, extend='both', format="%.1f",
                          ticks=np.linspace(cb_range[0], cb_range[1],
                                            cb_diff + 1.0),
                          cax=cbar_ax)
        cb.set_label('Speedup rel. to non-vectorized code')
        plt.savefig('{1}_Speedup_{0}.png'.format(legend[i], base_string),
                    dpi=400)
        plt.savefig('{1}_Speedup_{0}.pdf'.format(legend[i], base_string),
                    dpi=400)
        fig.clear()
        ax.clear()
        axhist.clear()
        plt.close(fig)

if __name__ == '__main__':
    main()
