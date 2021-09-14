#!/usr/bin/env/python
""" Tests for generating in silico NGS reads
"""

import unittest
from read_generators import *


class TestReadGenerator(unittest.TestCase):
    """ Tests for ReadGenerator
    """

    def setUpClass(cls):
        pass

    def tearDownClass(cls):
        pass

    def test_make_random_ends(self):
        pass

    def test_make_reads_from_fragment(self):
        pass

    def test_fastqs_from_features(self):
        pass



class TestDnaReadGenerator(unittest.TestCase):
    """ Tests for DnaReadGenerator
    """

    def setUpClass(cls):

        with read_enumerators.ReadEnumerator(r1=) as rg:
            cls.rg = rg


    def tearDownClass(cls):
        pass

    def test_weight_features(self):
        pass

    def test_make_random_ends(self):
        pass

    def test_make_reads_from_fragment(self):
        pass

    def test_fastqs_from_features(self):
        pass


class TestRnaReadGenerator(unittest.TestCase):
    """ Tests for RnaReadGenerator
    """

    def setUpClass(cls):
        with read_enumerators.ReadEnumerator(r1=) as rg:
            cls.rg = rg

    def tearDownClass(cls):
        pass

    def test_make_random_ends(self):
        pass

    def test_make_reads_from_fragment(self):
        pass

    def test_fastqs_from_features(self):
        pass