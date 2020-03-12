import exceRNApipeline.tasks.task_preprocess_stats as task
from unittest import TestCase
import sys
import re
import os
import shutil


_FILE_PATH = os.path.dirname(__file__)

class TestTaskPreProcessStats(TestCase):
    """"Test task_preprocess_stats.py"""
    def setUp(self):
        super().setUp()
        os.mkdir(os.path.join(_FILE_PATH, "_test_temp"))
        self.sys_argv = sys.argv
    
    def tearDown(self):
        super().tearDown()
        shutil.rmtree(os.path.join(_FILE_PATH, "_test_temp"))
        sys.argv = self.sys_argv
    
    def test_preprocess_stats(self):
        """Test fore preprocess stats"""
        sys.argv = [
            "",
            "-l", f"{_FILE_PATH}/data/sample1.htsStats.log",
                  f"{_FILE_PATH}/data/sample2.htsStats.log",
                  f"{_FILE_PATH}/data/sample3.htsStats.log",
                  f"{_FILE_PATH}/data/sample4.htsStats.log",
            "-n", "sample1", "sample2", "sample3", "sample4",
            "-o", f"{_FILE_PATH}/_test_temp/sample.tsv",
            "-p", f"{_FILE_PATH}/_test_temp/sample.png"
        ]
        task.main()
