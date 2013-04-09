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
  NDIGITS = 3 # number of digits after floating point to display
  __noInfoChars = re.compile(r'[\[\]\s\']') # characters which wouldn't be printed in VCF output
  __keysToPrint = ['imprecise', 'svtype', 'cpos', 'end', 'max',
                   'cend', 'svlen', 'cilen', 'conf', 'trachrom',
                   'trapos', 'tracpos', 'traend', 'tracend'] # keys that represents info in VCF

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    self._rname = rname
    self._rindex = sample.getRefIndex(rname)
    self._sample = sample
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

    for key in self.__keysToPrint:
      if key in info:
        if type(info[key]) == types.BooleanType:
          if info[key]:
            result += "%s;" % (key.upper())
        else:
          result += "%s=%s;" % (key.upper(), self.__noInfoChars.sub('', str(info[key])))

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

  @abstractmethod
  def add(self, variation):
    """
    Try to add variation into cluster and return if it fits into cluster
    """
    return False
