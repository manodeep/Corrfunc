"""
Example python code to call the 3 correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys
import re
import numpy as np
import time

# Import from current directory first,
# and then from the package. 
try: 
    import _countpairs
    if sys.version_info[0] >= 3:
        def rd(filename):
            with open(filename, encoding="utf-8") as f:
                r = f.read()
                
            return r
    else:
        def rd(filename):
            with open(filename) as f:
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


        import numpy as np
        def read_ascii(filename,return_dtype=None):
            if return_dtype is None:
                raise ValueError("Return data-type must be set and a valid numpy data-type")

            ### check if pandas is available - much faster to read in the data through pandas
            print("Reading in the data...")
            try:
                import pandas as pd
                df = pd.read_csv(file,header=None,engine="c",dtype={"x":return_dtype,"y":return_dtype,"z":return_dtype},delim_whitespace=True)
                x = np.asarray(df[0],dtype=return_dtype)
                y = np.asarray(df[1],dtype=return_dtype)
                z = np.asarray(df[2],dtype=return_dtype)
            except ImportError:
                print("Warning: Could not read in data with pandas -- due to error : {}. Falling back to slower numpy.".format(sys.exc_info()[0]))
                x,y,z = np.genfromtxt(file,dtype=return_dtype,unpack=True)

            return x,y,z

        def read_fastfood(filename,return_dtype=None):
            if return_dtype is None:
                raise ValueError("Return data-type must be set and a valid numpy data-type")

            import struct
            with open(filename, "rb") as f:
                skip1 = struct.unpack('@i',f.read(4))[0]
                idat  = struct.unpack('@iiiii',f.read(20))[0:5]
                skip2 = struct.unpack('@i',f.read(4))[0]
                assert skip1  == 20 and skip2 == 20,"fast-food file seems to be incorrect (reading idat)"
                ngal = idat[1]
                ## now read fdat
                skip1 = struct.unpack('@i',f.read(4))[0]
                fdat  = struct.unpack('@fffffffff',f.read(36))[0:9]
                skip2 = struct.unpack('@i',f.read(4))[0]
                assert skip1  == 36 and skip2 == 36,"fast-food file seems to be incorrect (reading fdat )"

                skip1 = struct.unpack('@i',f.read(4))[0]
                znow  = struct.unpack('@f',f.read(4))[0]
                skip2 = struct.unpack('@i',f.read(4))[0]
                assert skip1  == 4 and skip2 == 4,"fast-food file seems to be incorrect (reading redshift)"

                ### Now read the first 4 bytes for x positions
                skip1 = struct.unpack('@i',f.read(4))[0]
                assert skip1 == ngal*4 or skip1 == ngal*8, "fast-food file seems to be corrupt (padding bytes x)"

                ### seek back 4 bytes from current position
                f.seek(-4, 1)
                pos = {}
                for field in 'xyz':
                    skip1 = struct.unpack('@i',f.read(4))[0]
                    assert skip1 == ngal*4 or skip1 == ngal*8, "fast-food file seems to be corrupt (padding bytes for {}, pad = {} sizeof = {})".format(field, skip1, skip1//ngal)
                    input_dtype = np.float32 if skip1//ngal == 4 else np.float
                    array = np.fromfile(f, input_dtype, ngal)
                    skip2 = struct.unpack('@i',f.read(4))[0]
                    pos[field] = array if dtype is None else dtype(array)

            x = pos['x']
            y = pos['y']
            z = pos['z']

            return x,y,z

        if filebase is None:
            filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../tests/data/","gals_Mr19")
            ## Figure out the datatype, use the header file in the include directory
            ## because that is most likely correct (common.mk might have been modified
            ## but not recompiled)
            include_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "../../include/", "countpairs.h")
            includes = rd(include_file)
            vector_type = re.search(r'(\w+)\s*\*\s*rupp\s*\;', includes, re.I).group(1)
            allowed_types = {"float":np.float32,"double":np.float}
            if vector_type not in list(allowed_types.keys()):
                print("Error: Unknown precision={} found in header file {}. Allowed types are `{}'".format(vector_type,include_file,allowed_types))
                sys.exit()

            dtype = allowed_types[vector_type]
            allowed_exts = {'.ff' :read_fastfood,
                            '.txt':read_ascii,
                            '.dat':read_ascii,
                            '.csv':read_ascii
                            }


            for e in allowed_exts:
                if os.path.exists(filename+e):
                    f = allowed_exts[e]
                    x,y,z = f(filename+e, dtype)
                    return x,y,z
            raise IOError("Could not locate {} with any of these extensions = {}".format(filename, exts.keys()))
        else:
            ### Likely an user-supplied value
            if os.path.exists(filebase):
                extension = os.path.splitext(filebase)[1][1:].strip().lower()
                f = read_fastfood if '.ff' in extension else read_ascii

                ### default return is double
                x,y,z = f(filebase, np.float)
                return x,y,z
            else:
                raise IOError("Could not locate file {}",filebase)



        
except ImportError:
    from Corrfunc import _countpairs, rd, utils
    from .utils import read_catalog
    
tstart=time.time()
t0=tstart
x,y,z = read_catalog()
t1=time.time()    
print("Done reading the data - time taken = {0:10.1f} seconds.\nBeginning Correlation functions calculations".format(t1-t0))

boxsize=420.0
nthreads=4
pimax=40.0
binfile=os.path.join(os.path.dirname(os.path.abspath(__file__)),"../tests/","bins")
autocorr=1
numbins_to_print=5

print("Running 3-D correlation function DD(r)")
results_DD = _countpairs.countpairs(autocorr,nthreads,binfile,x,y,z,x,y,z)
print("\n#      **** DD(r): first {} bins  *******       ".format(numbins_to_print))
print("#      rmin        rmax       rpavg       npairs")
print("################################################")
for ibin in range(numbins_to_print):
    items = results_DD[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d}".format(items[0],items[1],items[2],items[3]))
print("------------------------------------------------")


print("\nRunning 2-D correlation function DD(rp,pi)")
results_DDrppi = _countpairs.countpairs_rp_pi(autocorr,nthreads,pimax,binfile,x,y,z,x,y,z)
print("\n#            ****** DD(rp,pi): first {} bins  *******      ".format(numbins_to_print))
print("#      rmin        rmax       rpavg     pi_upper     npairs")
print("###########################################################")
for ibin in range(numbins_to_print):
    items = results_DDrppi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")


print("\nRunning 2-D projected correlation function wp(rp)")
results_wp = _countpairs.countpairs_wp(boxsize,pimax,nthreads,binfile,x,y,z)
print("\n#            ******    wp: first {} bins  *******         ".format(numbins_to_print))
print("#      rmin        rmax       rpavg        wp       npairs")
print("##########################################################")
for ibin in range(numbins_to_print):
    items = results_wp[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")

print("\nRunning 3-D auto-correlation function xi(r)")
results_xi = _countpairs.countpairs_xi(boxsize,nthreads,binfile,x,y,z)
print("\n#            ******    xi: first {} bins  *******         ".format(numbins_to_print))
print("#      rmin        rmax       rpavg        xi       npairs")
print("##########################################################")
for ibin in range(numbins_to_print):
    items = results_xi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")

print("Done with all four correlation calculations.")

print("\nRunning VPF pN(r)")
rmax=10.0
nbin=10
nspheres=10000
num_pN=3
seed=-1
results_vpf = _countpairs.countspheres_vpf(rmax,nbin,nspheres,num_pN,seed,x,y,z)

print("\n#            ******    pN: first {} bins  *******         ".format(numbins_to_print))
print('#       r    ',end="")

for ipn in range(num_pN):
    print('        p{:0d}      '.format(ipn),end="")

print("")

print("###########",end="")
for ipn in range(num_pN):
    print('################',end="")
print("")


for ibin in range(numbins_to_print):
    items = results_vpf[ibin]
    print('{0:10.2f} '.format(items[0]),end="")
    for ipn in range(num_pN):
        print(' {0:15.4e}'.format(items[ipn+1]),end="")
    print("")

print("-----------------------------------------------------------")

tend=time.time()
print("Done with all THEORY python functions. Total time taken = {0:10.1f} seconds. Read-in time = {1:10.1f} seconds.".format(tend-tstart,t1-t0))


