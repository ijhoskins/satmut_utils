#!/usr/bin/env/python
"""Collection of string manipulation utilities."""

import itertools
import random
import string

__author__ = "ianjameshoskins"

DEFAULT_RANDOMER_LETTERS = string.ascii_uppercase + string.digits
DEFAULT_RANDOMER_DIGITS = string.digits


def make_random_str(str_len=10, prefix="", letters=DEFAULT_RANDOMER_DIGITS):
    """Make a random string from a selection of the letters.

    :param int str_len: length of randomer
    :param str prefix: optional prefix to include for string
    :param str letters: string of characters to use to generate randomer
    :return str: randomer string
    """

    randomer = prefix + "_".join(random.choice(letters) for _ in range(str_len))
    return randomer


def make_unique_ids(n_ids, str_len=10, prefix="", letters=DEFAULT_RANDOMER_DIGITS):
    """Make a list of randomer strings; a generator would be good here but could potentially yield non-unique randomers.

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

    https://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float/354130
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


def look_ahead(input_iterable, n_ahead=3):
    """Creates a look-ahead of an iterable and returns a list of tuples from i...i+n_ahead.

    :param iterable input_iterable: any type that implements an iter method
    :param int n_ahead: number of elements to look ahead
    :return list: list of tuples of the ith to the ith + n_ahead elements

    WARNING: per the docs, "Once tee() has made a split, the original iterable should not be used anywhere else;
    otherwise, the iterable could get advanced without the tee objects being informed."
    """

    # Generate n_ahead number of copies of the input iterable
    tee_iters = itertools.tee(input_iterable, n_ahead)

    # Shift the teed iterables then zip them together
    shifted_tee_iters = [list(t)[i:] for i, t in enumerate(tee_iters)]
    zipped_iters = itertools.zip_longest(*shifted_tee_iters)
    return zipped_iters
