#!/usr/bin/env python
from __future__ import print_function
import os,sys,re

min_py_major=2
min_py_minor=6

## Enforce minimum python version
if sys.version_info[0] < min_py_major or (sys.version_info[0] == min_py_major  and sys.version_info[1] < min_py_minor):
    sys.exit('Sorry, Python < {}.{} is not supported'.format(min_py_major,min_py_minor))

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

if sys.argv[-1] == "upload":
    os.system("python setup.py sdist upload")
    sys.exit()
    
        
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
name = re.search(r'DISTNAME\s*=\s*([a-z]+[A-z]+)', common, re.I).group(1)
major = re.search(r'MAJOR\s*=\s*(\d)',common, re.I).group(1)
minor = re.search(r'MINOR\s*=\s*(\d)',common, re.I).group(1)
patch = re.search(r'PATCHLEVEL\s*=\s*(\d)',common, re.I).group(1)
version = "{}.{}.{}".format(major,minor,patch)

setup(
    name=name,
    version=version,
    author='Manodeep Sinha',
    author_email='manodeep@gmail.com',
    maintainer='Manodeep Sinha',
    maintainer_email='manodeep@gmail.com',
    url='http://github.com/manodeep/Corrfunc',
    description='Blazing fast correlation functions on the CPU',
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: Linux, OSX',
        'Programming Language :: Python',
        ],
    license='MIT',
    keywords=['correlation functions','simulations','surveys','galaxies'],
#     platform=[
#     packages=['countpairs','countpairs_mocks'],
#     package_dir={'countpairs':'xi_theory/'}
    scripts=['xi_theory/python_bindings/call_correlation_functions.py','xi_mocks/python_bindings/call_correlation_functions_mocks.py'],
    include_package_data=True,
    package_data={
        '':['xi_theory/tests/bins',
            'xi_mocks/tests/bins','xi_mocks/tests/angular_bins',
            'xi_theory/tests/Mr19*','xi_theory/tests/cmass*','xi_theory/tests/data/*.ff','xi_theory/tests/data/*.txt',
            'xi_mocks/tests/Mr19*','xi_mocks/tests/data/*.dat','xi_mocks/tests/data/*.ff','xi_mocks/tests/data/*.txt',
            ],
        },
    install_requires=['setuptools','numpy>=1.9','python>={}.{}'.format(min_py_major,min_py_minor)],
    zip_safe=False,
    )

