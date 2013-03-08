#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "06.03. 2013"

import sys
sys.path.append("..")
import random
import unittest
from src.Sample import Sample

class TestSampleFunctions(unittest.TestCase):
    def setUp(self):
      pass

    def test_coverage(self):
      sample = Sample("../inputs/tests/sample_coverage.bam", "../inputs/references/GRCh37/human_g1k_v37.fasta.gz")
      self.assertEqual(sample.getCoverage('1', 39), 5)
      self.assertEqual(sample.getCoverage('1', 111), 5)
      self.assertEqual(sample.getCoverage('1', 224), 2)
      self.assertEqual(sample.getCoverage('2', 0), 1)
      self.assertEqual(sample.getCoverage('2', 117), 0)
      self.assertEqual(sample.getCoverage('2', 200), 0)
      self.assertEqual(sample.getCoverage('3', 78), 1)
      self.assertEqual(sample.getCoverage('8', 12), 0)
      self.assertEqual(sample.getCoverage('8', 121), 1)
      sample.close()

    def test_insertSize(self):
      sample = Sample("../inputs/tests/sample_insertSize.bam", "../inputs/references/GRCh37/human_g1k_v37.fasta.gz")
      self.assertEqual(sample.getSumInsertSize(), 1214626805)
      self.assertEqual(sample.getAvgInsertSize(), 242925361)
      sample.close()

    def tearDown(self):
      pass

if __name__ == '__main__':
    unittest.main(verbosity=2)
