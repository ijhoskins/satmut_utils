#!/usr/bin/env/python
"""Collection of file manipulation utilities."""

import glob
import gzip
import os
import shutil
import subprocess

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

    ext_res = ".".join([filename, ext])
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

    ext_res = "{}.{}".format(split[0], ext)
    return ext_res


def create_str_from_iter(iterable):
    """Creates a tab-delimited string from an iterable.

    :param iter iterable: iterable with each element to be a field
    :return str: iterable joined by FILE_DELIM and ending with FILE_NEWLINE
    """

    res = FILE_DELIM.join(list(map(str, iterable))) + FILE_NEWLINE
    return res


def gzip_file(filename):
    """Gzips a file.

    :param str filename: file path
    :return str: path of the gzipped file
    """

    if os.path.exists(filename):
        subprocess.call(("gzip", filename))
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


def get_file_list(pattern="*fastq*", indir=".", is_paired=True, recursive=False):
    """Gets a list of sample input files in indir that match the pattern.

    :param str pattern: glob-style pattern for matching. Default find zipped or unzipped FASTQs
    :param str indir: input directory to search for files
    :param bool is_paired: look for paired files (R1, R2) must exist in filenames. Default True.
    :param bool recursive: find matching files recursively? Default False.
    :return list: list of tuples; if paired end data, corresponding files will be paired; otherwise tuple will be len 1
    """

    if not is_paired:
        search_pattern = os.path.join(indir, pattern)
        matching_files = [tuple(e) for e in glob.glob(search_pattern, recursive=recursive)]
        return matching_files

    # Get paired files
    r1_search_pattern = os.path.join(indir, "*{}{}".format("[rR]1", pattern))
    r2_search_pattern = os.path.join(indir, "*{}{}".format("[rR]2", pattern))
    r1_matches = sorted(glob.glob(r1_search_pattern, recursive=recursive))
    r2_matches = sorted(glob.glob(r2_search_pattern, recursive=recursive))

    if len(r1_matches) != len(r2_matches):
        raise RuntimeError(
            "Unable to confidently pair files as R1 and R2 matches differed. "
            "R1 matches were %s and R2 matches were %s" % (",".join(r1_matches), ",".join(r2_matches)))

    paired_files = list(zip(r1_matches, r2_matches))
    return paired_files


def filter_file_by_keys(filter_file, key_file, filter_key_pos, outfile, delim=",", is_gzipped=False):
    """Filters a file by matching keys from a key file.

    :param str filter_file: file containing data to filter
    :param str key_file: file with keys, one per line
    :param int filter_key_pos: field number for the key in the filter file
    :param str outfile: output file to write to
    :param str delim: delimiter for the filter file. Default comma-delimited
    :param bool is_gzipped: is the filter_file gzipped? Default False.
    """

    open_func = open if not is_gzipped else gzip.open

    with open_func(filter_file, "r") as filt_fh, \
            open(key_file, "r") as key_fh, \
            open(outfile, "w") as out_fh:

        keys = {k.strip() for k in key_fh}

        for line in filt_fh:
            line_str = line.decode("UTF-8")
            line_split = line_str.split(delim)
            filter_key = line_split[filter_key_pos]
            if filter_key in keys:
                out_fh.write(line_str)
