#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "23.03. 2013"

from interface.interface import *
from Read import Read
from SplitPair import SplitPair
from SplitRead import SplitRead

class Paired:
  """
  Represents paired reads and single read when mate is None
  """
  _readStrand = True # strand of read
  _mateStrand = False # strand of mate
  rtype = enum(# type of read
               NORMAL=0,
               SINGLE=1,
               READ_SPLIT=2,
               MATE_SPLIT=3,
               FILTERED=4
              )

  def __init__(self, read, mate, reflengths, refnames, splitparts):
    """
    Initialize reads and find out type of reads
    """
    if mate and Paired.isFirst(mate, read):
      self.__read = Read(mate, True, Paired._readStrand, refnames)
      self.__mate = Read(read, False, Paired._mateStrand, refnames)
    else:
      self.__read = Read(read, True, Paired._readStrand, refnames)
      self.__mate = Read(mate, False, Paired._mateStrand, refnames)

    self.__qname = read.qname
    self.__reflengths = reflengths
    self.__splitpair = None

    if not self.__read.hasMinQuality(): # read doesn't have minimal mapping quality
      if self.__mate.isUnmapped() or not self.__mate.hasMinQuality(): # both don't have minimal mapping quality
        self.__type = Paired.rtype.FILTERED
      else: # make mate single
        self.__read = self.__mate
        self.__mate = None
        self.__type = Paired.rtype.SINGLE
    elif self.__mate.isUnmapped(): # mate unmapped
      self.__type = Paired.rtype.SINGLE
    elif not self.__mate.hasMinQuality(): # mate doesn't have minimal mapping quality
      self.__type = Paired.rtype.SINGLE
      self.__mate = None
    elif self.__read.isDuplicate(): # read is duplicate
      if self.__mate.isDuplicate(): # both are duplicated -> filter out
        self.__type = Paired.rtype.FILTERED
      else: # only read is duplicte -> make mate single
        self.__type = Paired.rtype.SINGLE
        self.__read = self.__mate
        self.__mate = None
    elif self.__mate.isDuplicate(): # mate is duplicate -> make read single
      self.__type = Paired.rtype.SINGLE
      self.__mate = None
    elif self.__read.isSplit(): # read is split
      if self.__mate.isInverted() or self.__mate.hasGaps() or self.isInterchromosomal(): # filter both
        self.__type = Paired.rtype.FILTERED
      else: # split read
        self.__type = Paired.rtype.READ_SPLIT
        self.__splitpair = SplitPair(False, self.__mate, SplitRead(self.__read, Paired._readStrand, splitparts))
    elif self.__mate.isSplit(): # mate is split
      if self.__read.isInverted() or self.__read.hasGaps() or self.isInterchromosomal(): # filter both
        self.__type = Paired.rtype.FILTERED
      else: # split mate
        self.__type = Paired.rtype.MATE_SPLIT
        self.__splitpair = SplitPair(True, self.__read, SplitRead(self.__mate, Paired._mateStrand, splitparts))
    elif not self.__read.hasGaps() and not self.__mate.hasGaps(): # paired without any gaps
      self.__type = Paired.rtype.NORMAL
    else:
      self.__type = Paired.rtype.FILTERED

    if self.__splitpair and (not self.__splitpair.splitread.hasMinQuality() or not self.__splitpair.splitread.hasMinLengths() or self.actualSize() <= 0):
      self.__type = Paired.rtype.FILTERED

  @staticmethod
  def isFirst(read, mate):
    """
    Test if read is before his mate
    """
    return (read.pos <= mate.pos and read.tid == mate.tid) or read.tid < mate.tid

  @property
  def qname(self):
    """
    Return query name
    """
    return self.__qname

  @property
  def read(self):
    """
    Return read
    """
    return self.__read

  @property
  def mate(self):
    """
    Return mate
    """
    return self.__mate

  @property
  def splitpair(self):
    """
    Return split read
    """
    return self.__splitpair

  def size(self):
    """
    Return insert size information from read
    """
    if not self.isNormal() or self.hasOverlap() or self.isRearranged() or self.__read.isInverted() or self.__mate.isInverted() or self.isInterchromosomal():
      return 0

    return self.__read.sam.tlen - (self.__read.end - self.__read.pos) - (self.__mate.end - self.__mate.pos)

  def actualSize(self):
    """
    Return counted insert size if template length is zero
    """
    size = self.size()

    if size or self.isSingle(): # size from read or read is single
      return size
    elif self.hasOverlap(): # overlapping pair
      if self.isRearranged(): # and also rearranged
        return self.__read.pos - self.__mate.end

      return self.__mate.pos - self.__read.end
    else: # normal or rearranged pair
      lengthBetween = 0

      if self.__read.tid != self.__mate.tid: # another chromosome
        lengthBetween = sum(self.__reflengths[self.__read.tid:self.__mate.tid])

      result = lengthBetween + self.__mate.pos - self.__read.end - 1

      if self.isRearranged():
        return -(result + self.__read.len + self.__mate.len)

      return result

  def hasOverlap(self):
    """
    Test if reads overlap
    """
    return self.__read.pos <= self.__mate.pos and \
            ((self.__read.end <= self.__mate.end and \
              self.__mate.pos <= self.__read.end) or \
             (self.__mate.end <= self.__read.end))

  def isNormal(self):
    """
    Test if reads are normal
    """
    return self.__type == Paired.rtype.NORMAL

  def isSingle(self):
    """
    Test if read is single
    """
    return self.__type == Paired.rtype.SINGLE

  def isReadSplit(self):
    """
    Test if read is split
    """
    return self.__type == Paired.rtype.READ_SPLIT

  def isMateSplit(self):
    """
    Test if mate is split
    """
    return self.__type == Paired.rtype.MATE_SPLIT

  def isFiltered(self):
    """
    Test if read is filtered
    """
    return self.__type == Paired.rtype.FILTERED

  def isRearranged(self):
    """
    Test if reads are rearranged
    """
    return self.__read.isInverted() and self.__mate.isInverted()

  def isInterchromosomal(self):
    """
    Test if reads are on different chromosomes
    """
    return self.__read.tid != self.__mate.tid
