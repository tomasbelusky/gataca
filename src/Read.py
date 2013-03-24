#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "23.03. 2013"

class Read:
  """
  Represents paired reads and single read when mate is None
  """
  _rearrangedRef = None # reference to function which tests if reads are rearranged
  _invertedRef = None #  reference to function which tests if read is inverted

  def __init__(self, read, mate, reflengths):
    """
    Initialize reads
    """
    self.__read = read
    self.__mate = mate
    self.__reflengths = reflengths

  @staticmethod
  def isFirst(read, mate):
    """
    Test if read is before his mate
    """
    return (read.pos <= mate.pos and read.tid == mate.tid) or read.tid < mate.tid

  def read(self):
    """
    Return read
    """
    return self.__read

  def mate(self):
    """
    Return mate
    """
    return self.__mate

  def isSingle(self):
    """
    Test if read is single
    """
    return not self.__mate

  def length(self):
    """
    Return template length information from read
    """
    return self.__read.tlen

  def actualLength(self):
    """
    Return counted template length if is 0 in read
    """
    if self.__read.tlen or self.isSingle():
      return self.__read.tlen
    else:
      if not mate.rlen:
        mate.rlen = len(mate.seq)

      if self.__read.tid == self.__mate.tid: # same chromosome
        return mate.pos + mate.rlen - read.pos
      else: # another chromosome
        lengthBetween = sum(self.__reflengths[read.tid:mate.tid])
        return lengthBetween + mate.pos + mate.rlen - read.pos

  def _isInvertedFR(self, read, first):
    """
    Test if read is inverted in forward/reverse policy
    """
    return (first and read.is_reverse) or (not first and not read.is_reverse)

  def _isInvertedRF(self, read, first):
    """
    Test if read is inverted in reverse/forward policy
    """
    return (first and not read.is_reverse) or (not first and read.is_reverse)

  def isReadInverted(self):
    """
    Test if read is inverted
    """
    return self._invertedRef(self.__read, True)

  def isMateInverted(self):
    """
    Test if mate is inverted
    """
    return self._invertedRef(self.__mate, False)

  def _isRearrangedFR(self, read, mate):
    """
    Test if reads are rearranged in forward/reverse policy
    """
    return read.is_reverse and not mate.is_reverse

  def _isRearrangedRF(self, read, mate):
    """
    Test if reads are rearranged in reverse/forward policy
    """
    return not read.is_reverse and mate.is_reverse

  def isRearranged(self):
    """
    Test if reads are rearranged
    """
    return self._rearrangedRef(self.__read, self.__mate)
