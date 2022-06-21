"""
Example python code to call the 3 correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
try:
    from future.utils import bytes_to_native_str
except ImportError:
    print("\n\tPlease run python -m pip install Corrfunc before using "
          "the 'Corrfunc' package\n")
    raise

from os.path import dirname, abspath, exists, splitext, join as pjoin
import time
import numpy as np

import _countpairs
from _countpairs import \
    countpairs as DD,\
    countpairs_rp_pi as DDrppi,\
    countpairs_wp as wp,\
    countpairs_xi as xi,\
    countspheres_vpf as vpf,\
    countpairs_s_mu as DDsmu


def read_text_file(filename, encoding="utf-8"):
    """
    Reads a file under python3 with encoding (default UTF-8).
    Also works under python2, without encoding.
    Uses the EAFP (https://docs.python.org/2/glossary.html#term-eafp)
    principle.
    """
    try:
        with open(filename, 'r', encoding) as f:
            r = f.read()
    except TypeError:
        with open(filename, 'r') as f:
            r = f.read()
    return r


def read_catalog(filebase=None):
    """
    Reads a galaxy/randoms catalog.

    :param filebase: (optional)
        The fully qualified path to the file. If omitted, reads the
        theory galaxy catalog under ../tests/data/

    Returns:
    * ``x y z`` - Unpacked numpy arrays compatible with the installed
    version of ``Corrfunc``.

    **Note** If the filename is omitted, then first the fast-food file
    is searched for, and then the ascii file. End-users should always
    supply the full filename.
    """

    def read_ascii_catalog(filename, return_dtype=None):
        if return_dtype is None:
            msg = 'Return data-type must be set and a valid numpy data-type'
            raise ValueError(msg)

        # check if pandas is available - much faster to read in the data
        # using pandas
        print("Reading in the data...")
        try:
            import pandas as pd
        except ImportError:
            pd = None

        if pd is not None:
            df = pd.read_csv(filename, header=None,
                             engine="c",
                             dtype={"x": return_dtype,
                                    "y": return_dtype,
                                    "z": return_dtype},
                             delim_whitespace=True)
            x = np.asarray(df[0], dtype=return_dtype)
            y = np.asarray(df[1], dtype=return_dtype)
            z = np.asarray(df[2], dtype=return_dtype)

        else:
            x, y, z, _ = np.genfromtxt(filename, dtype=return_dtype,
                                       unpack=True)

        return x, y, z

    def read_fastfood_catalog(filename, return_dtype=None, need_header=None):
        if return_dtype is None:
            msg = "Return data-type must be set and a valid numpy data-type"
            raise ValueError(msg)

        import struct
        with open(filename, "rb") as f:
            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            idat = struct.unpack(bytes_to_native_str(b'@iiiii'),
                                 f.read(20))[0:5]
            skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == 20 and skip2 == 20,\
                "fast-food file seems to be incorrect (reading idat)"
            ngal = idat[1]

            if need_header is not None:
                # now read fdat
                skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                fdat = struct.unpack(bytes_to_native_str(b'@fffffffff'),
                                     f.read(36))[0:9]
                skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                assert skip1 == 36 and skip2 == 36,\
                    "fast-food file seems to be incorrect (reading fdat )"

                skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                znow = struct.unpack(bytes_to_native_str(b'@f'), f.read(4))[0]
                skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                assert skip1 == 4 and skip2 == 4,\
                    "fast-food file seems to be incorrect (reading redshift)"
            else:
                fdat_bytes = 4 + 36 + 4
                znow_bytes = 4 + 4 + 4
                # seek over the fdat + znow fields + padding bytes
                # from current position
                f.seek(fdat_bytes + znow_bytes, 1)

            # read the padding bytes for the x-positions
            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == ngal * 4 or skip1 == ngal * 8, \
                "fast-food file seems to be corrupt (padding bytes)"

            # seek back 4 bytes from current position
            f.seek(-4, 1)
            pos = {}
            for field in 'xyz':
                skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                assert skip1 == ngal * 4 or skip1 == ngal * 8, \
                    "fast-food file seems to be corrupt (padding bytes a)"
                # the next division must be the integer division
                input_dtype = np.float32 if skip1 // ngal == 4 else np.float64
                array = np.fromfile(f, input_dtype, ngal)
                skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
                pos[field] = array if return_dtype == input_dtype \
                    else [return_dtype(a) for a in array]

        x = pos['x']
        y = pos['y']
        z = pos['z']

        if need_header is not None:
            return idat, fdat, znow, x, y, z
        else:
            return x, y, z

    if filebase is None:
        filename = pjoin(dirname(abspath(_countpairs.__file__)),
                         "../tests/data/", "gals_Mr19")
        dtype = np.float32
        allowed_exts = {'.ff': read_fastfood_catalog,
                        '.txt': read_ascii_catalog,
                        '.dat': read_ascii_catalog,
                        '.csv': read_ascii_catalog
                        }

        for e in allowed_exts:
            if exists(filename + e):
                f = allowed_exts[e]
                x, y, z = f(filename + e, dtype)
                return x, y, z
        raise IOError("Could not locate {0} with any of these extensions \
        = {1}".format(filename, allowed_exts.keys()))
    else:
        # Likely an user-supplied value
        if exists(filebase):
            extension = splitext(filebase)[1]
            f = read_fastfood_catalog if '.ff' in extension else read_ascii_catalog

            # default return is double
            x, y, z = f(filebase, np.float64)
            return x, y, z

        raise IOError("Could not locate file {0}", filebase)


def main():
    tstart = time.time()
    t0 = tstart
    x, y, z = read_catalog()
    w = np.ones((1,len(x)), dtype=x.dtype)
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")
    boxsize = 420.0
    nthreads = 4
    pimax = 40.0
    binfile = pjoin(dirname(abspath(_countpairs.__file__)),
                    "../tests/", "bins")
    autocorr = 1
    periodic = 1
    numbins_to_print = 5

    print("Running 3-D correlation function DD(r)")
    results_DD, _ = DD(autocorr=autocorr,
                       nthreads=nthreads,
                       binfile=binfile,
                       X1=x, Y1=y, Z1=z, weights1=w,
                       weight_type='pair_product',
                       periodic=periodic,
                       boxsize=(boxsize, boxsize, boxsize),
                       verbose=True,
                       output_ravg=True,
                       isa=-1)

    print("\n#      **** DD(r): first {0} bins  *******       "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg       npairs       weight_avg")
    print("##################################################################")
    for ibin in range(numbins_to_print):
        items = results_DD[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d} {4:13.4f}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("------------------------------------------------------------------")


    print("Running 3-D correlation function DD(r) with different bin refinement")
    results_DD, _ = DD(autocorr=autocorr,
                       nthreads=nthreads,
                       binfile=binfile,
                       X1=x, Y1=y, Z1=z,
                       periodic=periodic,
                       boxsize=(boxsize, boxsize, boxsize),
                       verbose=True,
                       output_ravg=True,
                       xbin_refine_factor=3,
                       ybin_refine_factor=1,
                       zbin_refine_factor=2)

    print("\n#      **** DD(r): first {0} bins  *******       "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg       npairs")
    print("################################################")
    for ibin in range(numbins_to_print):
        items = results_DD[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d}"
              .format(items[0], items[1], items[2], items[3]))
    print("------------------------------------------------")


    print("\nRunning 2-D correlation function DD(rp,pi)")
    results_DDrppi, _ = DDrppi(autocorr=autocorr,
                               nthreads=nthreads,
                               pimax=pimax,
                               binfile=binfile,
                               X1=x, Y1=y, Z1=z, weights1=w,
                               X2=x, Y2=y, Z2=z, weights2=w,
                               periodic=periodic,
                               boxsize=(boxsize, boxsize, boxsize),
                               verbose=True,
                               output_rpavg=True,
                               weight_type='pair_product')

    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs    weight_avg")
    print("#########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-------------------------------------------------------------------------")

    mu_max = 0.5
    nmu_bins = 10

    print("\nRunning 2-D correlation function DD(s,mu)")
    results_DDsmu, _ = DDsmu(autocorr=autocorr,
                             nthreads=nthreads,
                             binfile=binfile,
                             mu_max=mu_max,
                             nmu_bins=nmu_bins,
                             X1=x,
                             Y1=y,
                             Z1=z,
                             weights1=w,
                             weight_type='pair_product',
                             verbose=True,
                             periodic=periodic,
                             boxsize=(boxsize, boxsize, boxsize),
                             output_savg=True)
    print("\n#            ****** DD(s,mu): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      smin        smax       savg     mu_max     npairs    weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDsmu[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")



    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp, _, _ = wp(boxsize=boxsize, pimax=pimax, nthreads=nthreads,
                          binfile=binfile, X=x, Y=y, Z=z,
                          weights=w, weight_type='pair_product',
                          verbose=True, output_rpavg=True)
    print("\n#            ******    wp: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        wp       npairs     weight_avg")
    print("#########################################################################")
    for ibin in range(numbins_to_print):
        items = results_wp[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-------------------------------------------------------------------------")


    print("\nRunning 3-D auto-correlation function xi(r)")
    results_xi, _ = xi(boxsize=boxsize, nthreads=nthreads, binfile=binfile,
                       X=x, Y=y, Z=z, weights=w,
                       weight_type='pair_product',
                       verbose=True, output_ravg=True)
    print("\n#            ******    xi: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        xi       npairs    weight_avg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_xi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")
    print("Done with all four correlation calculations.")


    print("\nRunning VPF pN(r)")
    rmax = 10.0
    nbin = 10
    nspheres = 10000
    num_pN = 8
    seed = -1
    results_vpf, _ = vpf(rmax=rmax, nbins=nbin, nspheres=nspheres,
                         num_pN=num_pN, seed=seed, X=x, Y=y, Z=z, verbose=True,
                         periodic=periodic,
                         boxsize=(boxsize, boxsize, boxsize),
                         )

    print("\n#            ******    pN: first {0} bins  *******         "
          .format(numbins_to_print))
    print('#       r    ', end="")

    for ipn in range(num_pN):
        print('        p{0:0d}      '.format(ipn), end="")

    print("")

    print("###########", end="")
    for ipn in range(num_pN):
        print('################', end="")
    print("")

    for ibin in range(numbins_to_print):
        items = results_vpf[ibin]
        print('{0:10.2f} '.format(items[0]), end="")
        for ipn in range(num_pN):
            print(' {0:15.4e}'.format(items[ipn + 1]), end="")
        print("")

    print("-----------------------------------------------------------")

    tend = time.time()
    print("Done with all functions. Total time taken = {0:10.1f} seconds. \
    Read-in time = {1:10.1f} seconds.".format(tend - tstart, t1 - t0))

if __name__ == "__main__":
    main()
