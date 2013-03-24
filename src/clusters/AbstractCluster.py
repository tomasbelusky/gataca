#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "13.03. 2013"

import re
from abc import abstractmethod
from src.Variation import Variation

class AbstractCluster:
  """
  Represents cluster of reads
  """
  noInfoChars = re.compile(r'[\[\]\s]') # characters which wouldn't be printed in VCF output

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    self._rname = rname
    self._rindex = sample.getRefIndex(rname)
    self._sample = sample
    self._variations = []
    self._start = 0
    self._end = 0
    self._alleles = []

  def getStart(self):
    """
    Return index of leftmost base from variations in cluster
    """
    return self._start

  def getEnd(self):
    """
    Return index of rightmost base from variations in cluster
    """
    return self._end

  def add(self, variation):
    """
    Add variation into cluster
    """
    self._variations.append(variation)
    self._process()

  def remove(self, variation):
    """
    Remove variation from cluster
    """
    if variation in self.__variations:
      self._variations.remove(variation)
      self._process()

  @abstractmethod
  def _process(self):
    """
    Join variations together for printing them in VCF output
    """
    pass

  @abstractmethod
  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    pass
