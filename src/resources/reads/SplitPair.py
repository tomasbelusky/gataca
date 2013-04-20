#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "11.04. 2013"

class SplitPair:
  """
  Represents pair where one read has split parts
  """

  def __init__(self, readFirst, read, splitread):
    """
    Initialize variables
    """
    self.__readFirst = readFirst
    self.__read = read
    self.__splitread = splitread

  @property
  def read(self):
    """
    Return not split read
    """
    return self.__read

  @property
  def splitread(self):
    """
    Return split read
    """
    return self.__splitread

  def isReadFirst(self):
    """
    Test if not split read is before split read
    """
    return self.__readFirst

  def partsEnclosePair(self):
    """
    Test if split parts enclose not split read
    """
    return (self.__readFirst and self.__splitread.remapped.pos < self.__read.pos) or \
           (not self.__readFirst and self.__read.end < self.__splitread.remapped.end)

  def hasOverlap(self):
    """
    Test if remapped split part overlap with not split read
    """
    return (self.__read.pos <= self.__splitread.remapped.pos and \
            ((self.__read.end <= self.__splitread.remapped.end and \
              self.__splitread.remapped.pos <= self.__read.end) or \
             (self.__splitread.remapped.end <= self.__read.end))) or \
            (self.__splitread.remapped.pos <= self.__read.pos and \
            ((self.__splitread.remapped.end <= self.__read.end and \
              self.__read.pos <= self.__splitread.remapped.end) or \
             (self.__read.end <= self.__splitread.remapped.end)))
