#!/usr/bin/env python3
"""Collection of string manipulation utilities."""

import random
import string

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

DEFAULT_RANDOMER_LETTERS = string.ascii_uppercase + string.digits
DEFAULT_RANDOMER_DIGITS = string.digits


def make_random_str(str_len=10, prefix="", letters=DEFAULT_RANDOMER_DIGITS):
    """Make a random string from a selection of the letters.

    :param int str_len: length of randomer
    :param str prefix: optional prefix to include for string
    :param str letters: string of characters to use to generate randomer
    :return str: randomer string
    """

    randomer = prefix + "".join(random.choice(letters) for _ in range(str_len))
    return randomer


def make_unique_ids(n_ids, str_len=10, prefix="", letters=DEFAULT_RANDOMER_DIGITS):
    """Make a list of randomer strings.

    :param int n_ids: number of IDs to create
    :param int str_len: length of randomer
    :param str prefix: optional prefix to include for string
    :param str letters: string of characters to use to generate randomer
    :return list: randomer strings
    """

    uniq_ids = set()
    while len(uniq_ids) < n_ids:
        uniq_ids.add(make_random_str(str_len, prefix, letters))
    uniq_ids = list(uniq_ids)

    return uniq_ids


def is_number(s):
    """Checks if a string can be cast to a number.

    :param str s: string to check
    :return bool: whether or not the string represents a number

    """
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


def none_or_str(value):
    """Converts None string to proper type.

    :param str value: str
    :return str | None: appropriate type
    """

    if value == 'None':
        return None
    return value


def none_or_int(value):
    """Converts None string to proper type.

    :param str value: str
    :return int | None: appropriate type
    """

    if value == 'None':
        return None
    return int(value)
