#!/usr/bin/env python3
"""Collection of file manipulation utilities."""

import os
import shutil
import subprocess

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

FILE_EXT_DELIM = "."
FILE_DELIM = "\t"
FILE_NEWLINE = "\n"
FILE_SPACE = " "
GZ_EXTENSION = "gz"


def safe_remove(paths, force_remove=False):
    """Safely remove a file or directory.

    :param tuple paths: filename(s) to remove
    :param bool force_remove: force removal of non-empty directories
    """

    for path in paths:
        if path != "/":
            if os.path.isdir(path):
                try:
                    os.rmdir(path)
                except OSError:
                    if force_remove:
                        shutil.rmtree(path)
                    else:
                        raise NotImplementedError("Directory to be removed is not empty. For safety, did not remove. "
                                                  "If you would like to proceed, set force_remove=True.")
            elif os.path.exists(path):
                os.unlink(path)


def flush_files(fls):
    """Flushes buffer to disk.

    :param tuple fls: file objects to flush
    """

    for fl in fls:
        fl.flush()
        os.fsync(fl.fileno())
        fl.seek(0)


def add_extension(filename, ext):
    """Adds an extension.

    :param str filename: file path
    :param str ext: extension to add
    """

    ext_res = ".".join((filename, ext,))
    return ext_res


def remove_extension(filename):
    """Removes the extension of a filename.

    :param str filename: file path
    :return str: filename without extension
    """

    split = os.path.splitext(filename)
    ext_res = split[0]
    return ext_res


def get_extension(filename):
    """Gets the extension of a filename.

    :param str filename: file path
    :return str: extension
    """

    ext = str(filename).split(FILE_EXT_DELIM)[-1:][0]
    return ext


def replace_extension(filename, ext, ignore_exts=(".gz", ".bz", ".bz2",)):
    """Replaces extension of a filename.

    :param str filename: file path
    :param str ext: extension to add
    :param tuple ignore_exts: extensions to strip before replacing the extension
    :return str: filepath with new extension
    """

    split = os.path.splitext(filename)

    if split[1] in set(ignore_exts):
        split = os.path.splitext(split[0])

    ext_res = add_extension(split[0], ext)
    return ext_res


def gzip_file(filename, force=False):
    """Gzips a file.

    :param str filename: file path
    :param bool force: force overwrite? Default False.
    :return str: path of the gzipped file
    """

    if os.path.exists(filename):

        gzip_call = ["gzip"]
        if force:
            gzip_call.append("--force")
        gzip_call.append(filename)

        subprocess.call(tuple(gzip_call))
        return add_extension(filename, GZ_EXTENSION)
    else:
        raise RuntimeError("Filename %s not found." % filename)


def gunzip_file(filename):
    """Gzips a file.

    :param str filename: file path
    :return str: path of the unzipped file
    """

    if os.path.exists(filename):
        subprocess.call(("gunzip", "--keep", filename))
        return remove_extension(filename)
    else:
        raise RuntimeError("Filename %s not found." % filename)


def is_gzipped(filename):
    """Determines if a file is gzipped based on the extension.

    :param str filename: file path
    :return bool: whether the file is gzipped or not
    """

    split = os.path.splitext(filename)

    if split[1] == ".gz":
        return True

    return False
