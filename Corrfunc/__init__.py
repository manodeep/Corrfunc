#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Corrfunc is a set of high-performance routines for
computing clustering statistics on a distribution of
points.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os
import sys
import warnings

__version__ = "2.5.0"
__author__ = "Manodeep Sinha <manodeep@gmail.com>"


if sys.version_info[0] < 3:
    warnings.warn('Python 2 support is deprecated and will be dropped in Corrfunc v2.6', FutureWarning)

try:
        __CORRFUNC_SETUP__
except NameError:
    __CORRFUNC_SETUP__ = False

if not __CORRFUNC_SETUP__:
    from . import io
    from . import utils
    from . import theory
    from . import mocks

    utils.check_runtime_env()


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


def write_text_file(filename, contents, encoding="utf-8"):
    """
    Writes a file under python3 with encoding (default UTF-8).
    Also works under python2, without encoding.
    Uses the EAFP (https://docs.python.org/2/glossary.html#term-eafp)
    principle.
    """
    try:
        with open(filename, 'w', encoding) as f:
            f.write(contents)

    except TypeError:
        with open(filename, 'w') as f:
            f.write(contents)


def which(program, mode=os.F_OK | os.X_OK, path=None):
    """
    Mimics the Unix utility which.
    For python3.3+, shutil.which provides all of the required functionality.
    An implementation is provided in case shutil.which does
    not exist.

    :param program: (required) string
           Name of program (can be fully-qualified path as well)
    :param mode: (optional) integer flag bits
           Permissions to check for in the executable
           Default: os.F_OK (file exists) | os.X_OK (executable file)
    :param path: (optional) string
           A custom path list to check against. Implementation taken from
           shutil.py.

    Returns:
           A fully qualified path to program as resolved by path or
           user environment.
           Returns None when program can not be resolved.
    """
    try:
        from shutil import which as shwhich
        return shwhich(program, mode, path)
    except ImportError:
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, _ = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            if path is None:
                path = os.environ.get("PATH", os.defpath)
            if not path:
                return None

            path = path.split(os.pathsep)
            for pathdir in path:
                pathdir = pathdir.strip('"')
                exe_file = os.path.join(pathdir, program)
                if is_exe(exe_file):
                    return exe_file

        return None
