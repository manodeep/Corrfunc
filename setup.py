#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import print_function

import os
from os.path import join as pjoin, dirname, abspath, commonprefix
import sys
from sys import version_info
import re
import subprocess

if sys.version_info[0] < 3:
        import __builtin__ as builtins
else:
    import builtins

# Make sure we are running on posix (Linux, Unix, MAC OSX)
if os.name != 'posix':
    sys.exit("Sorry, Windows is not supported")

base_url = "https://github.com/manodeep/Corrfunc"
projectname = 'Corrfunc'

# global variables
version = ''
compiler = ''

builtins.__CORRFUNC_SETUP__ = True

# partial import
import Corrfunc
from Corrfunc import read_text_file, write_text_file, which

# numpy 1.7 supports python 2.4-2.5; python 3.1-3.3.
try:
    from setuptools import setup, Extension, find_packages
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup, Extension

    def find_packages(path='.'):
        ret = []
        for root, _, files in os.walk(path):
            if '__init__.py' in files:
                ret.append(re.sub('^[^A-z0-9_]+', '', root.replace('/', '.')))
                return ret
    from distutils.command.build_ext import build_ext


def strip_line(line, sep=os.linesep):
    """
    Removes occurrence of character (sep) from a line of text
    """

    try:
        return line.strip(sep)
    except TypeError:
        return line.decode('utf-8').strip(sep)


def run_command(command, **kwargs):
    proc = subprocess.Popen(command, shell=True,
                            **kwargs)
    stdout, stderr = proc.communicate(None)
    status = proc.wait()
    if status:
        msg = "command = {0} failed with stdout = {1} stderr = {2} "\
              "status {3:d}\n".format(command, stdout, stderr, status)
        raise RuntimeError(msg)
    return stdout, stderr


def get_dict_from_buffer(buf, keys=['DISTNAME', 'MAJOR',
                                    'MINOR', 'PATCHLEVEL',
                                    'PYTHON',
                                    'MIN_PYTHON_MAJOR',
                                    'MIN_PYTHON_MINOR',
                                    'MIN_NUMPY_MAJOR',
                                    'MIN_NUMPY_MINOR']):
    """
    Parses a string buffer for key-val pairs for the supplied keys.

    Returns: Python dictionary with all the keys (all keys in buffer
             if None is passed for keys) with the values being a list
             corresponding to each key.

    Note: Return dict will contain all keys supplied (if not None).
          If any key was not found in the buffer, then the value for
          that key will be [] such that dict[key] does not produce
          a KeyError.

    Slightly modified from:
    "http://stackoverflow.com/questions/5323703/regex-how-to-"\
    "match-sequence-of-key-value-pairs-at-end-of-string

    """


    pairs = dict()
    if keys is None:
        keys = "\S+"
        regex = re.compile(r'''
        \n            # all key-value pairs are on separate lines
        \s*           # there might be some leading spaces
        (             # start group to return
        (?:{0}\s*)    # placeholder for tags to detect '\S+' == all
        \s*:*=\s*     # optional spaces, optional colon, = , optional spaces
        .*            # the value
        )             # end group to return
        '''.format(keys), re.VERBOSE)
        validate = False
    else:
        keys = [k.strip() for k in keys]
        regex = re.compile(r'''
        \n            # all key-value pairs are on separate lines
        \s*           # there might be some leading spaces
        (             # start group to return
        (?:{0}\s*)    # placeholder for tags to detect '\S+' == all
        \s*:*=\s*     # optional spaces, optional colon, = , optional spaces
        .*            # the value
        )             # end group to return
        '''.format('|'.join(keys)), re.VERBOSE)
        validate = True
        for k in keys:
            pairs[k] = []

    matches = regex.findall(buf)
    for match in matches:
        key, val = match.split('=', 1)

        # remove colon and leading/trailing whitespace from key
        key = (strip_line(key, ':')).strip()

        # remove newline and leading/trailing whitespace from value
        val = (strip_line(val)).strip()
        if validate and key not in keys:
            msg = "regex produced incorrect match. regex pattern = {0} "\
                  "claims key = [{1}] while original set of search keys "\
                  "= {2}".format(regex.pattern, key, '|'.join(keys))
            raise AssertionError(msg)
        pairs.setdefault(key, []).append(val)

    return pairs


def replace_first_key_in_makefile(buf, key, replacement, outfile=None):
    '''
    Replaces first line in 'buf' matching 'key' with 'replacement'.
    Optionally, writes out this new buffer into 'outfile'.

    Returns: Buffer after replacement has been done
    '''

    regexp = re.compile(r'''
    \n\s*           # there might be some leading spaces
    (               # start group to return
    (?:{0}\s*)      # placeholder for tags to detect '\S+' == all
    \s*:*=\s*       # optional spaces, optional colon, = , optional spaces
    .*              # the value
    )               # end group to return
    '''.format(key), re.VERBOSE)
    matches = regexp.findall(buf)
    if matches is None:
        msg = "Could not find key = {0} in the provided buffer. "\
              "Pattern used = {1}".format(key, regexp.pattern)
        raise ValueError(msg)

    # Only replace the first occurence
    newbuf = regexp.sub(replacement, buf, count=1)
    if outfile is not None:
        write_text_file(outfile, newbuf)

    return newbuf


def requirements_check():
    common_mk_file = pjoin(dirname(abspath(__file__)),
                           "common.mk")
    common = read_text_file(common_mk_file)

    common_dict = get_dict_from_buffer(common)
    name = common_dict['DISTNAME'][0]
    major = common_dict['MAJOR'][0]
    minor = common_dict['MINOR'][0]
    patch = common_dict['PATCHLEVEL'][0]
    if name is None or major is None or minor is None or patch is None:
        msg = "ERROR: Did not find at least one of the keys "\
              "(DISTNAME, MAJOR, MINOR, PATCHLEVEL) in 'common.mk'.\n"\
              "Checks can not run - aborting installation. "\
              "projectname = {1} major = {2} minor = {3} patch = {4}\n\n"\
              "You can fix this by re-issuing git clone {0}".\
              format(base_url, projectname, major, minor, patch)
        raise AssertionError(msg)

    if projectname != name:
        msg = 'Mis-match between C project name and python project name'\
              'C claims project = {0} while python has {1}'.\
              format(name, projectname)
        raise AssertionError(msg)

    global version
    version = "{0}.{1}.{2}".format(major, minor, patch)
    # Check that version matches
    if Corrfunc.__version__ != version:
        msg = "ERROR: Version mis-match. Python version found = {0} \
        while C version claims {1}".format(Corrfunc.__version__, version)
        raise AssertionError(msg)

    # Okay common.mk has been updated to use current python
    # for building the extensions as required. Now check for
    # min. python version
    min_py_major = int(common_dict['MIN_PYTHON_MAJOR'][0])
    min_py_minor = int(common_dict['MIN_PYTHON_MINOR'][0])

    # Enforce minimum python version
    if version_info[0] < min_py_major or \
       (version_info[0] == min_py_major and version_info[1] < min_py_minor):
        msg = "Sorry. Found python {0}.{1} but minimum required "\
              "python version is {2}.{3}".format(version_info[0],
                                                 version_info[1],
                                                 min_py_major,
                                                 min_py_minor)
        raise AssertionError(msg)

    # Since arbitrary python can be used even within the Makefile
    # make sure that the current python executable is the same as the
    # one specified in common.mk. Easiest way is to replace,
    # but no need to do so if creating a source distribution
    if 'sdist' not in sys.argv:
        this_python = sys.executable
        key = "PYTHON"
        replacement = '\n{0}:={1}'.format(key, this_python)
        common = replace_first_key_in_makefile(common, key, replacement,
                                                   common_mk_file)

    # Check if CC is in argv:
    CC = "CC"
    for iarg, arg in enumerate(sys.argv):
        if CC not in arg:
            continue

        if '=' in arg:
            # user passed `CC=/path/to/compiler`
            # therefore, split on '=' to get the
            # compiler (i.e, 'CC') and the value
            # (i.e, the name/path of compiler)
            key, value = arg.strip().split('=')
        else:
            # Space-separated or spaces and an '=' sign
            key = arg.strip()
            if key != CC:
                msg = "Something strange has happened. Expected to find "\
                      "a custom compiler from the command-line but \n"\
                      "found command-line argument '{0}' (that matches "\
                      "pattern ['CC=/path/to/compiler']). Parsing "\
                      "produced CC={1}".format(arg, key)
                raise ValueError(msg)

            check_arg = iarg+1
            # Is there an "=" sign (i.e., `CC = /path/to/compiler`)
            # or did the user simply pass `CC /path/to/compiler`
            if check_arg >= len(sys.argv):
                msg = "Found compiler key = {} but could not locate "\
                      "compiler value - no further command-line "\
                      "parameters were passed.\nPlease pass the "\
                      "custom compiler name either `CC=compiler`"\
                      "or as `CC=/path/to/compiler`".format(key)
                raise ValueError(msg)

            # The user could have specified `CC =compiler` or
            # `CC = compiler`. The following 'if' condition checks
            # for the first case, the 'else' checks for the second
            # case (`CC = compiler`)
            if '=' in sys.argv[check_arg] and \
               sys.argv[check_arg].strip() != '=':
                _, value = sys.argv[check_arg].strip().split('=')
            else:
                # Otherwise, there was white-space separated '='
                # we can delete that command-line argument containing
                # just the '=' sign.
                del sys.argv[check_arg]
                # should be parsing the compiler value now
                if check_arg >= len(sys.argv):
                    msg = "Found compiler key = CC but could not locate "\
                          "compiler value (either as `CC=/path/to/CC` "\
                          "or as `CC /path/to/CC`"
                    raise ValueError(msg)

                value = sys.argv[check_arg].strip()

            # this deletes the argument containing the compiler name
            del sys.argv[check_arg]

        if key != CC or value == '':
            msg = "Something strange has happened. Expected to find a "\
                  "custom compiler from the command-line but found \n"\
                  "command-line argument '{0}' (that matches pattern "\
                  "['CC=/path/to/compiler']). Parsing produced CC={1} "\
                  "and $CC={2}".format(arg, key, value)

            raise ValueError(msg)

        # check if value is a valid compiler
        full_compiler = which(value)
        if full_compiler is None:
            msg = "Found compiler = '{0}' on the command-line but '{0}' "\
                  "can not be resolved from the shell.\n"\
                  "Please specify CC=/path/to/compiler in the "\
                  "python -m pip setup.py call.".format(value)
            raise ValueError(msg)

        replacement = '\n{0}:={1}'.format(CC, value)
        replace_first_key_in_makefile(common, CC,
                                      replacement, common_mk_file)

        global compiler
        compiler = value

        # Delete the 'CC' key, the compiler name and the '='
        # have already been deleted
        del sys.argv[iarg]
        break

    return common_dict


class BuildExtSubclass(build_ext):
    def build_extensions(self):
        # Everything has already been configured - so just call make
        for ext in self.extensions:
            sources = ext.sources

            # Blatantly copied from the Android distutils setup.py
            # But then modified to make it python3 compatible!
            if sources is None or not isinstance(sources, list):
                msg = "in 'ext_modules' option (extension {0}),"\
                      "sources={1} must be present and must be "\
                      "a list of source filename".\
                      format(ext.name, sources)
                raise AssertionError(msg)

            if len(sources) == 1:
                ext_dir = dirname(sources[0])
            else:
                # not debugged - could be wrong
                # (where sources from multiple directories
                # are specified simultaneously)
                ext_dir = commonprefix(sources)

            # global variable compiler is set if passed in
            # command-line
            extra_string = ''
            if compiler != '':
                extra_string = 'CC={0}'.format(compiler)
            command = "cd {0} && make {1}".format(ext_dir, extra_string)
            run_command(command)

            import shutil
            import errno
            try:
                os.makedirs(dirname(self.get_ext_fullpath(ext.name)))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            # full_name = '{0}.so.{1}'.format(pjoin(ext_dir, ext.name),
            #                                 version)
            full_build_name = '{0}.{1}'.format(
                self.get_ext_fullpath(ext.name), version)

            full_name = '{0}.so'.format(pjoin(ext_dir, ext.name))
            full_build_name = '{0}'.format(self.get_ext_fullpath(ext.name))
            pkg_sourcedir = '{0}'.format(pjoin(ext_dir, '../../Corrfunc'))
            pkg_in_srcdir = '{0}/{1}.so'.format(pkg_sourcedir, ext.name)

            shutil.copyfile(full_name, full_build_name)

            # just copy the newly created library in the Corrfunc module directory.
            # Installed Corrfunc version will automatically get the extensions
            # os.remove(pkg_in_srcdir)
            # os.symlink('{0}'.format(pjoin('../', full_name)),
            #           pkg_in_srcdir)
            shutil.copyfile(full_name, pkg_in_srcdir)


def generate_extensions(python_dirs):
    extensions = []
    for pdir in python_dirs:
        mk = read_text_file(pjoin(dirname(abspath(__file__)),
                                  pdir, "Makefile"))
        makefile_dict = get_dict_from_buffer(mk, ['PROJECT',
                                                  'SOURCES'])
        project = makefile_dict['PROJECT'][0].strip()
        src_files = makefile_dict['SOURCES']
        if project is None or src_files is None:
            msg = "In directory = {0}, can not locate either the project "\
                  "name or the list of source files."\
                  "Got project = {1} and sources = {2}."\
                  .format(pdir, project, src_files)
            raise AssertionError(msg)

        sources = [pjoin(pdir, f) for f in src_files]
        # print("Found project = {0} in dir = {1} with sources = {2}\n".
        #       format(project, pdir, sources))
        ext = Extension("{0}".format(project),
                        sources=sources)
        extensions.append(ext)
    return extensions


# Only python >= 3.5 supports the recursive glob, hence
# defining the function that works on all reasonable pythons
# http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
def recursive_glob(rootdir='.', patterns=None):
    import fnmatch
    if not patterns:
        patterns = ['*']
    return [pjoin(looproot, filename)
            for looproot, _, filenames in os.walk(rootdir)
            for filename in filenames for p in patterns
            if fnmatch.fnmatch(filename, p)]


def install_required():
    install_args = ['build', 'build_ext', 'build_clib',
                    'install', 'bdist']

    for arg in install_args:
        if arg in sys.argv:
            return True

    return False


# Taken from numpy setup.py
def setup_packages():
    '''
    Custom setup for Corrfunc package.

    Optional: Set compiler via 'CC=/path/to/compiler' or
              'CC /path/to/compiler' or 'CC = /path/to/compiler'
              All the CC options are removed from sys.argv after
              being parsed.
    '''

    # protect the user in case they run python setup.py not from root directory
    src_path = dirname(abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # create a list of the python extensions
    python_dirs = ["theory/python_bindings",
                   "mocks/python_bindings"]
    extensions = generate_extensions(python_dirs)

    # check requirement for extensions and set the compiler if specified
    # in command-line
    common_dict = requirements_check()

    # Some command options require headers/libs to be generated
    # so that the following dirs_patters supplies them.
    if install_required():
        from distutils.sysconfig import get_config_var
        if get_config_var('SHLIB_EXT') != '".so"' and version_info[0] == 2:
            msg = "The extensions all get the `.so` automatically. "\
                  "However, python expects the extension to be `{0}`"\
                  .format(get_config_var('SHLIB_EXT'))
            raise ValueError(msg)

        # global variable compiler is set if passed in
        # command-line
        extra_string = ''
        if compiler != '':
            extra_string = 'CC={0}'.format(compiler)

        command = "make libs {0}".format(extra_string)
        run_command(command)

    else:
        # not installing. Check if creating source distribution
        # in that case run distclean to delete auto-generated C
        # files
        if 'sdist' in sys.argv:
            command = "make distclean"
            run_command(command)


    # find all the data-files required.
    # Now the lib + associated header files have been generated
    # and put in lib/ and include/
    # This step must run after ``make install``
    dirs_patterns = {'theory/tests/data': ['*.ff', '*.txt',
                                              '*.txt.gz', '*.dat'],
                     'mocks/tests/data': ['*.ff', '*.txt',
                                             '*.txt.gz', '*.dat'],
                     'theory/tests': ['Mr19*', 'bins*', 'cmass*'],
                     'mocks/tests': ['Mr19*', 'bins*', 'angular_bins*'],
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
    long_description = read_text_file('README.rst')
    min_np_major = int(common_dict['MIN_NUMPY_MAJOR'][0])
    min_np_minor = int(common_dict['MIN_NUMPY_MINOR'][0])

    # All book-keeping is done.
    # base_url = "https://github.com/manodeep/Corrfunc"
    classifiers = ['Development Status :: 4 - Beta',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: MIT License',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: C',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10']
    metadata = dict(
            name=projectname,
            version=version,
            author='Manodeep Sinha',
            author_email='manodeep@gmail.com',
            maintainer='Manodeep Sinha',
            maintainer_email='manodeep@gmail.com',
            url=base_url,
            download_url='{0}/archive/{1}-{2}.tar.gz'.format(
                    base_url, projectname, version),
            description='Blazing fast correlation functions on the CPU',
            long_description=long_description,
            classifiers=classifiers,
            license='MIT',
            # Solaris might work, Windows will almost certainly not work
            platforms=["Linux", "Mac OSX", "Unix"],
            keywords=['correlation functions', 'simulations',
                      'surveys', 'galaxies'],
            provides=[projectname],
            packages=find_packages(),
            ext_package=projectname,
            ext_modules=extensions,
            package_data={'': data_files},
            include_package_data=True,
            setup_requires=['setuptools',
                            'numpy>={0}.{1}'.format(min_np_major,
                                                    min_np_minor)],
            install_requires=['numpy>={0}.{1}'.format(min_np_major,
                                                      min_np_minor),
                              'future',
                              'wurlitzer'],
            python_requires='>=2.7,!=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*, !=3.6.*,, <4',
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
    del builtins.__CORRFUNC_SETUP__
