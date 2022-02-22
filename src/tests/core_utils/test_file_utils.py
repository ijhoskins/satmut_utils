#!/usr/bin/env python3
"""Tests for core_utils.file_utils."""

import tempfile
import unittest

import core_utils.file_utils as fu
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR


class TestFileUtils(unittest.TestCase):
    """Tests for file utilities."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestFileUtils."""

        cls.temp_dir = tempfile.mkdtemp()
        cls.test_file = tempfile.NamedTemporaryFile(suffix=".test.file", delete=False, dir=cls.temp_dir).name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestFileUtils."""

        fu.safe_remove((cls.temp_dir,), force_remove=True)

    def test_safe_remove_dir_not_empty(self):
        """Test that we prohibit deletion of a non-empty directory by default."""

        with self.assertRaises(NotImplementedError):
            fu.safe_remove((self.temp_dir,))

    def test_safe_remove_file(self):
        """Test that we can remove a file safely."""

        test_file = tempfile.NamedTemporaryFile(suffix=".test.file", delete=False).name
        fu.safe_remove((test_file,))
        self.assertFalse(os.path.exists(test_file))

    def test_safe_remove_dir_force(self):
        """Test that we can force removal of a directory that's not empty."""

        # Create a new tempdir so we keep the initialized one intact for test_safe_remove_dir_not_empty
        temp_dir = tempfile.mkdtemp(suffix=".test.dir")
        _ = tempfile.NamedTemporaryFile(suffix=".test.file", delete=False, dir=temp_dir).name

        fu.safe_remove((temp_dir,), force_remove=True)
        self.assertFalse(os.path.exists(temp_dir))

    def test_flush_files(self):
        """Test that we can flush/sync file buffer to disc."""

        with open(self.test_file, "w+") as test_file:
            test_file.write("woohoo")
            fu.flush_files((test_file,))
            obs = test_file.read()

        self.assertTrue(obs == "woohoo")

    def test_remove_extension(self):
        """Test that we can properly remove a file extension."""

        obs = fu.remove_extension("fakefile.txt")
        self.assertEqual(obs, "fakefile")

    def test_get_extension(self):
        """Test that we can properly get a file extension."""

        obs = fu.get_extension("fakefile.txt")
        self.assertEqual(obs, "txt")

    def test_replace_extension(self):
        """Test that we can properly replace a file extension."""

        obs = fu.replace_extension("fakefile.txt", "bam")
        self.assertEqual(obs, "fakefile.bam")
