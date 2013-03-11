#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "06.03. 2013"

import sys
sys.path.append("..")
import pysam
import random
import unittest
from src.Sample import Sample

class TestSampleFunctions(unittest.TestCase):
    def setUp(self):
      self.__reference = pysam.Fastafile("../inputs/references/GRCh37/human_g1k_v37.fasta")

    def test_coverage(self):
      sample = Sample("../inputs/tests/sample_coverage.bam", self.__reference)
      sample.preprocessing()
      self.assertEqual(sample.getCoverage(0, 39), 5)
      self.assertEqual(sample.getCoverage(0, 111), 5)
      self.assertEqual(sample.getCoverage(0, 224), 2)
      self.assertEqual(sample.getCoverage(1, 0), 1)
      self.assertEqual(sample.getCoverage(1, 117), 0)
      self.assertEqual(sample.getCoverage(1, 200), 0)
      self.assertEqual(sample.getCoverage(2, 78), 1)
      self.assertEqual(sample.getCoverage(7, 12), 0)
      self.assertEqual(sample.getCoverage(7, 121), 1)
      sample.close()

    def test_insertSize(self):
      sample = Sample("../inputs/tests/sample_insertSize.bam", self.__reference)
      sample.preprocessing()
      self.assertEqual(sample.getAvgInsertSize(), 172.2)
      sample.close()

    def test_isFirst(self):
      filename = "../inputs/tests/sample_isFunctions.bam"
      sample = Sample(filename, self.__reference)
      reads = pysam.Samfile(filename)
      results = []

      for read in reads.fetch():
        results.append(sample.isFirst(read, reads.mate(read)))

      self.assertEqual(results, [True, True, False, False, True, False])
      sample.close()

    def test_isInvert(self):
      filename = "../inputs/tests/sample_isFunctions.bam"
      sample = Sample(filename, self.__reference)
      reads = pysam.Samfile(filename)
      results = []

      for read in reads.fetch():
        first = sample.isFirst(read, reads.mate(read))
        results.append(sample.isInverted(read, first))

      self.assertEqual(results, [True, False, False, True, True, True])
      sample.close()

    def test_isRearranged(self):
      filename = "../inputs/tests/sample_isFunctions.bam"
      sample = Sample(filename, self.__reference)
      reads = pysam.Samfile(filename)
      results = []

      for read in reads.fetch():
        results.append(sample.isRearranged(read, reads.mate(read)))

      self.assertEqual(results, [False, False, False, False, True, True])
      sample.close()

    def tearDown(self):
      pass

if __name__ == '__main__':
    unittest.main(verbosity=2)
