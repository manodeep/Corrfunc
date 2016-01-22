#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import (absolute_import, print_function)
import os
import sys
import re
from Corrfunc.utils import rd

## Make sure we are running on posix (Linux, Unix, MAC OSX)    
if os.name != 'posix':
    sys.exit("Sorry, Windows is not supported")


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

common = rd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "common.mk"))
name = re.search(r'DISTNAME\s*:*=\s*(\w+)', common).group(1)
major = re.search(r'MAJOR\s*:*=\s*(\d)',common).group(1)
minor = re.search(r'MINOR\s*:*=\s*(\d)',common).group(1)
patch = re.search(r'PATCHLEVEL\s*:*=\s*(\d)',common).group(1)
version = "{}.{}.{}".format(major,minor,patch)
min_py_major = int(re.search(r'MIN_PYTHON_MAJOR\s*:=\s*(\d)',common).group(1))
min_py_minor = int(re.search(r'MIN_PYTHON_MINOR\s*:=\s*(\d)',common).group(1))


min_numpy_major = int(re.search(r'MIN_NUMPY_MAJOR\s*:=\s*(\d)',common).group(1))
min_numpy_minor = int(re.search(r'MIN_NUMPY_MINOR\s*:=\s*(\d)',common).group(1))
   
## Enforce minimum python version
if sys.version_info[0] < min_py_major or (sys.version_info[0] == min_py_major  and sys.version_info[1] < min_py_minor):
    raise RuntimeError('Sorry. Found python {}.{} but minimum required python version is {}.{}'.format(sys.version_info[0],sys.version_info[1],min_py_major,min_py_minor))

## numpy 1.7 supports python 2.4-2.5; python 3.1-3.3. 

try:
    from setuptools import setup, Extension, find_packages
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup, Extension, find_packages
    from distutils.command.build_ext import build_ext

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()
    
        
            
class build_ext_subclass( build_ext ):
    def build_extensions(self):
        ### Everything has already been configured within Make - so just call make
        import subprocess
        for ext in self.extensions:
            sources = ext.sources
            
            ## Blatantly copied from the Android distutils setup.py
            if sources is None:
                raise Exception, \
                    ("in 'ext_modules' option (extension '%s')," +
                     "'sources must be present and must be '" + 
                     "a list of source filename") % ext.name
            
            if len(sources)==1:
                ext_dir = os.path.dirname(sources[0])
            else:
                ### not debugged - might be wrong
                ext_dir = os.path.commonprefix(sources)
                
            command = "cd {} && make -j2".format(ext_dir)
            # print("about to execute command `{}`. sources = {}".format(command,sources))
            proc = subprocess.Popen(command, stderr=subprocess.STDOUT, shell=True)
            output, stderr = proc.communicate(input)
            status = proc.wait()
            if status:
                raise Exception("command = {} failed with status {:d}".format(command,status),
                                output, status)

            import shutil
            import errno
            try:
                os.makedirs(os.path.dirname(self.get_ext_fullpath(ext.name)))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            shutil.copyfile('{}.so'.format(os.path.join(ext_dir,ext.name)),self.get_ext_fullpath(ext.name))
            # print("Made extension {}.so in path = {} ".format(ext.name,self.get_ext_fullpath(ext.name)))


python_dirs = ["xi_theory/python_bindings",
                "xi_mocks/python_bindings"]
extensions = []
for pdir in python_dirs:
    mk = rd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    pdir, "Makefile"))
    project       = re.search(r'PROJECT\s*:*=\s*(\w+)', mk).group(1)
    src_files     = re.findall(r'SOURCES\s*:*=\s*(\w+\.c)', mk)

    sources = [os.path.join(pdir,f) for f in src_files]
    # print("Found project = {} in dir = {} with sources = {}".format(project,pdir,sources))
    ext = Extension("{}".format(project),
                    sources=sources,
                    )
    
    extensions.append(ext)


### Taken from numpy setup.py    
def setup_packages():

    ## protect the user in case they run python setup.py not from root directory
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)
    
    metadata = dict(
        name=name,
        version=version,
        author='Manodeep Sinha',
        author_email='manodeep@gmail.com',
        maintainer='Manodeep Sinha',
        maintainer_email='manodeep@gmail.com',
        url='http://github.com/manodeep/Corrfunc',
        description='Blazing fast correlation functions on the CPU',
        long_description=rd('README.md'),
        classifiers = [
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            "Intended Audience :: Science/Research",
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: Linux, OSX',
            'Programming Language :: C',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            ],
        license='MIT',
        ### Solaris might work, Windows will almost certainly not work
        platforms = [ "Linux", "Mac OS-X", "Unix"],
        keywords=['correlation functions','simulations','surveys','galaxies'],
        packages=find_packages(),
        ext_package=name,
        ext_modules=extensions,
        scripts=['Corrfunc/call_correlation_functions.py','Corrfunc/call_correlation_functions_mocks.py'],
        include_package_data=True,
        package_data={
            '':['xi_theory/tests/bins',
                'xi_mocks/tests/bins','xi_mocks/tests/angular_bins',
                'xi_theory/tests/Mr19*','xi_theory/tests/cmass*','xi_theory/tests/data/*.ff','xi_theory/tests/data/*.txt',
                'xi_mocks/tests/Mr19*','xi_mocks/tests/data/*.dat','xi_mocks/tests/data/*.ff','xi_mocks/tests/data/*.txt',
                'xi_theory/xi_of_r/*.c','xi_theory/xi_of_r/*.h','xi_theory/xi_of_r/Makefile',
                'xi_theory/xi_rp_pi/*.c','xi_theory/xi_rp_pi/*.h','xi_theory/xi_rp_pi/Makefile',
                'xi_theory/wp/*.c','xi_theory/wp/*.h','xi_theory/wp/Makefile',
                'xi_theory/xi/*.c','xi_theory/xi/*.h','xi_theory/xi/Makefile',
                'xi_theory/vpf/*.c','xi_theory/vpf/*.h','xi_theory/vpf/Makefile',
                'xi_theory/tests/*.c','xi_theory/tests/*.h','xi_theory/tests/Makefile',
                'xi_theory/python_bindings/*.c','xi_theory/python_bindings/*.h','xi_theory/python_bindings/Makefile',
                ],
            },
        install_requires=['setuptools','numpy>={}.{}'.format(min_numpy_major,min_numpy_minor),'python>={}.{}'.format(min_py_major,min_py_minor)],
        zip_safe=False,
        cmdclass = {'build_ext': build_ext_subclass },
        )


    ### Now the actual setup
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_packages()
