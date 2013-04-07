#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from abc import abstractmethod
import re
import types

class AbstractCluster:
  """
  Represents cluster of base variations
  """
  __noInfoChars = re.compile(r'[\[\]\s\']') # characters which wouldn't be printed in VCF output

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    self._rname = rname
    self._rindex = sample.getRefIndex(rname)
    self._sample = sample
    self._variations = []
    self._start = 0
    self._actualStart = 0
    self._end = 0

  def infoString(self, info):
    """
    Create string from info fields
    """
    if not info:
      return ""

    result = ""

    for key, value in info.items(): # print common info
      if type(value) == types.BooleanType:
        if value:
          result += "%s;" % (key.upper())
      else:
        result += "%s=%s;" % (key.upper(), self.__noInfoChars.sub('', str(value)))

    return result[:-1]

  def getActualStart(self):
    """
    Return actual index of leftmost base from variations in cluster
    """
    return self._actualStart

  def getStart(self):
    """
    Return possible index of leftmost base from variations in cluster
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
    return

  @abstractmethod
  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    return False
