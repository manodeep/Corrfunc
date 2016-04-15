#!/usr/bin/env python
# -*- encoding: utf-8 -*-


from __future__ import (division, print_function)

import os
import os.path as path
import fnmatch
import sys
from sys import version_info
import re

# partial import
import Corrfunc
from Corrfunc import rd

# Make sure we are running on posix (Linux, Unix, MAC OSX)
if os.name != 'posix':
    sys.exit("Sorry, Windows is not supported")


common = rd(path.join(path.dirname(path.abspath(__file__)),
                      "common.mk"))
name = re.search(r'DISTNAME\s*:*=\s*(\w+)', common).group(1)
major = re.search(r'MAJOR\s*:*=\s*(\d)', common).group(1)
minor = re.search(r'MINOR\s*:*=\s*(\d)', common).group(1)
patch = re.search(r'PATCHLEVEL\s*:*=\s*(\d)', common).group(1)
version = "{0}.{1}.{2}".format(major, minor, patch)

# Check that version matches
if Corrfunc.__version__ != version:
    msg = "ERROR: Version mis-match. Python version found = {0} \
    while C version claims {1}".format(Corrfunc.__version__, version)
    sys.exit(msg)

min_py_major = int(re.search(r'MIN_PYTHON_MAJOR\s*:=\s*(\d)', common).group(1))
min_py_minor = int(re.search(r'MIN_PYTHON_MINOR\s*:=\s*(\d)', common).group(1))
min_np_major = int(re.search(r'MIN_NUMPY_MAJOR\s*:=\s*(\d)', common).group(1))
min_np_minor = int(re.search(r'MIN_NUMPY_MINOR\s*:=\s*(\d)', common).group(1))

# Enforce minimum python version
if version_info[0] < min_py_major or \
   (version_info[0] == min_py_major and version_info[1] < min_py_minor):
    raise RuntimeError('Sorry. Found python {0}.{1} but minimum required \
    python version is {2}.{3}'.format(version_info[0],
                                      version_info[1],
                                      min_py_major, min_py_minor))
# numpy 1.7 supports python 2.4-2.5; python 3.1-3.3.
try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command.build_ext import build_ext

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


def run_command(command):
    import subprocess
    # print("about to execute command `{0}`. sources = {1}"
    # .format(command, sources))
    proc = subprocess.Popen(command, stderr=subprocess.STDOUT, shell=True)
    output, stderr = proc.communicate(input)
    status = proc.wait()
    if status:
        raise Exception("command = {0} failed with output = {1} status {2:d}\n"
                        .format(command, output, status))


class BuildExtSubclass(build_ext):
    def build_extensions(self):
        # Everything has already been configured - so just call make
        for ext in self.extensions:
            sources = ext.sources

            # Blatantly copied from the Android distutils setup.py
            # But then modified to make it python3 compatible!
            if sources is None:
                raise Exception("in 'ext_modules' option (extension '%s')," +
                                "'sources must be present and must be '" +
                                "a list of source filename") % ext.name

            if len(sources) == 1:
                ext_dir = path.dirname(sources[0])
            else:
                # not debugged - likely to be wrong
                ext_dir = path.commonprefix(sources)

            command = "cd {0} && make ".format(ext_dir)
            run_command(command)

            import shutil
            import errno
            try:
                os.makedirs(path.dirname(self.get_ext_fullpath(ext.name)))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            # full_name = '{0}.so.{1}'.format(path.join(ext_dir, ext.name),
            #                                 version)
            full_build_name = '{0}.{1}'.format(
                self.get_ext_fullpath(ext.name), version)

            full_name = '{0}.so'.format(path.join(ext_dir, ext.name))
            full_build_name = '{0}'.format(self.get_ext_fullpath(ext.name))
            pkg_sourcedir = '{0}'.format(path.join(ext_dir, '../../Corrfunc'))
            pkg_in_srcdir = '{0}/{1}.so'.format(pkg_sourcedir, ext.name)

            shutil.copyfile(full_name, full_build_name)

            # just copy the newly created library in the Corrfunc module directory.
            # Installed Corrfunc version will automatically get the extensions
            #os.remove(pkg_in_srcdir)
            #os.symlink('{0}'.format(path.join('../', full_name)),
            #           pkg_in_srcdir)
            shutil.copyfile(full_name, pkg_in_srcdir)

def generate_extensions(python_dirs):
    extensions = []
    for pdir in python_dirs:
        mk = rd(path.join(path.dirname(path.abspath(__file__)),
                          pdir, "Makefile"))
        project = re.search(r'PROJECT\s*:*=\s*(\w+)', mk).group(1)
        src_files = re.findall(r'SOURCES\s*:*=\s*(\w+\.c)', mk)

        sources = [path.join(pdir, f) for f in src_files]
        # print("Found project = {0} in dir = {1} with sources = {2}".format(
        #     project, pdir, sources))
        ext = Extension("{0}".format(project),
                        sources=sources,
                        )

        extensions.append(ext)
    return extensions


# Only python >= 3.5 supports the recursive glob, hence
# defining the function that works on all reasonable pythons
# http://stackoverflow.com/questions/2186525/use-a-glob-to-
# find-files-recursively-in-python
def recursive_glob(rootdir='.', patterns=['*']):
    return [path.join(looproot, filename)
            for looproot, _, filenames in os.walk(rootdir)
            for filename in filenames for p in patterns
            if fnmatch.fnmatch(filename, p)]


# Taken from numpy setup.py
def setup_packages():

    # protect the user in case they run python setup.py not from root directory
    src_path = path.dirname(path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # create a list of the python extensions
    python_dirs = ["xi_theory/python_bindings",
                   "xi_mocks/python_bindings"]
    extensions = generate_extensions(python_dirs)

    # only run this if not creating source dist
    if "sdist" not in sys.argv:
        command = "make install"
        run_command(command)

    # find all the data-files required.
    # Now the lib + associated header files have been generated
    # and put in lib/ and include/
    # This step must run after ``make install``
    dirs_patterns = {'xi_theory/tests/data': ['*.ff', '*.txt',
                                              '*.txt.gz', '*.dat'],
                     'xi_mocks/tests/data': ['*.ff', '*.txt',
                                             '*.txt.gz', '*.dat'],
                     'xi_theory/tests': ['Mr19*', 'bins*', 'cmass*'],
                     'xi_mocks/tests': ['Mr19*', 'bins*', 'angular_bins*'],
                     'include': ['count*.h'],
                     'lib': ['libcount*.a']
                     }
    data_files = []
    for d in dirs_patterns:
        patterns = dirs_patterns[d]
        f = recursive_glob(d, patterns)
        data_files.extend(f)

    # change them to be relative to package dir rather than root
    data_files = ["../{0}".format(d) for d in data_files]
    long_description = rd('README.rst')

    # All book-keeping is done.
    base_url = "https://github.com/manodeep/Corrfunc"
    classifiers = ['Development Status :: 4 - Beta',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: MIT License',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: C',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5']
    metadata = dict(
        name=name,
        version=version,
        author='Manodeep Sinha',
        author_email='manodeep@gmail.com',
        maintainer='Manodeep Sinha',
        maintainer_email='manodeep@gmail.com',
        url=base_url,
        download_url='{0}/archive/{1}-{2}.tar.gz'.format(
            base_url, name, version),
        description='Blazing fast correlation functions on the CPU',
        long_description=long_description,
        classifiers=classifiers,
        license='MIT',
        # Solaris might work, Windows will almost certainly not work
        platforms=["Linux", "Mac OSX", "Unix"],
        keywords=['correlation functions', 'simulations',
                  'surveys', 'galaxies'],
        packages=[name],
        ext_package=name,
        ext_modules=extensions,
        package_data={'': data_files},
        include_package_data=True,
        install_requires=['setuptools',
                          'numpy>={0}.{1}'.format(min_np_major, min_np_minor),
                          'future'],
        zip_safe=False,
        cmdclass={'build_ext': BuildExtSubclass})

    # Now the actual setup
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_packages()
