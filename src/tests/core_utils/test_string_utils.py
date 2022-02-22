#!/usr/bin/env python3
"""Tests for core_utils.string_utils."""

import unittest

import core_utils.string_utils as su


class TestStringUtils(unittest.TestCase):
    """Tests for string utilities."""

    def test_make_random_str(self):
        """Test that we can make a randomer string."""

        str_len = 9
        prefix = "woohoo"
        obs = su.make_random_str(str_len, prefix)
        self.assertTrue(obs.startswith(prefix) and len(obs) == len(prefix) + str_len)

    def test_make_unique_ids(self):
        """Test that we can make a set of unique IDs."""

        n_ids = 5
        obs = su.make_unique_ids(n_ids)
        self.assertEqual(len(obs), n_ids)

    def test_is_number_positive(self):
        """Test that we can identify cast-able strings."""

        test_numbers = ["1.0", "2", ".3"]
        obs = [su.is_number(n) for n in test_numbers]
        self.assertTrue(all(obs))

    def test_is_number_negative(self):
        """Test that we can identify non-cast-able strings."""

        non_number = "A"
        obs = su.is_number(non_number)
        self.assertFalse(obs)
