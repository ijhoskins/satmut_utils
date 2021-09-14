#!/usr/bin/env/python
"""Python script template."""

import sys
import logging
import argparse

# Great APIs for direct calls to underlying C libraries for samtools and bedtools (no subprocess calls necessary)
# import pysam
# import pybedtools

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

# Log messages from this script will be prefixed by the script filename
__logger = logging.getLogger(__name__)

# Consider putting the following in a logging config file
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
__formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(__formatter)
__logger.addHandler(__fhandler)

FIELD_DELIM = "\t"
VCF_HEADER_CHAR = "#"
SAM_HEADER_CHAR = "@"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="")

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-i", "--infile", type=str, help='Input file')
    parser.add_argument("-f", "--flag", action="store_true", help='Flag on')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


# Keep the workflow out of main() so workflow can be imported as a whole from other modules
def workflow(infile, flag):
    """Runs the workflow.

    :param str infile: input filename
    :param bool flag: flag for option
    :return:
    """

    # Always use a context manager (with statement) so that open files are closed and automatic cleanup occurs whenever
    # an error occurs within the context
    with open(infile, "r") as in_fh:

        for i, line in enumerate(in_fh):
            line_strip = line.strip("\n")

            if len(line_strip) == 0:
                __logger.warning("Skipped empty line %i" % (i + 1))
                continue

            # Useful for skipping VCF and SAM headers
            if line_strip.startswith(VCF_HEADER_CHAR) or line_strip.startswith(SAM_HEADER_CHAR):
                __logger.warning("Skipped header line %i: %s" % (i + 1, line_strip))
                continue

            tokens = line_strip.split(FIELD_DELIM)  # file fields in a list


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])
    workflow(infile=parsed_args["infile"], flag=parsed_args["flag"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
