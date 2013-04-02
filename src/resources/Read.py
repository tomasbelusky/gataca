#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "23.03. 2013"

from interface.interface import *

class Read:
  """
  Represents paired reads and single read when mate is None
  """
  __minimalQuality = 30 # minimal quality to not be filtered out
  _rearrangedRef = None # reference to function which tests if reads are rearranged
  _readInvertedRef = None # reference to function which tests if read is inverted
  _mateInvertedRef = None # reference to function which tests if mate is inverted
  rtype = enum(# type of read
               NORMAL=0,
               SINGLE=1,
               READ_SPLIT=2,
               MATE_SPLIT=3,
               UNUSABLE=4,
               FILTERED=5
               )
  ctype = enum(# cigar operations
               ALIGNMENT=0,
               INSERTION=1,
               DELETION=2,
               SKIPPED=3,
               SOFTCLIP=4,
               HARDCLIP=5,
               PADDING=6,
               MATCH=7,
               MISMATCH=8)

  def __init__(self, read, mate, reflengths, refnames):
    """
    Initialize reads and find out type of reads
    """
    if mate and Read.isFirst(mate, read):
      self.__read = mate
      self.__mate = read
    else:
      self.__read = read
      self.__mate = mate

    self.__reflengths = reflengths
    self.__refnames = refnames
    self.__readEnd = Read.calculateReadEnd(self.__read)
    self.__mateEnd = Read.calculateReadEnd(self.__mate)
    readHasGaps = set(dict(self.__read.cigar)).intersection((Read.ctype.SKIPPED, Read.ctype.SOFTCLIP, Read.ctype.HARDCLIP, Read.ctype.PADDING))

    if self.__mate:
      mateHasGaps = set(dict(self.__mate.cigar)).intersection((Read.ctype.SKIPPED, Read.ctype.SOFTCLIP, Read.ctype.HARDCLIP, Read.ctype.PADDING))

    if not self.hasMinQuality(self.__read): # read doesn't have minimal mapping quality
      if not self.__mate or not self.hasMinQuality(self.__mate): # both don't have minimal mapping quality
        self.__type = Read.rtype.FILTERED
      else: # make mate single
        self.__read = self.__mate
        self.__mate = None
        self.__type = Read.rtype.SINGLE
    elif not self.__mate: # no mate or mate unmapped
      self.__type = Read.rtype.SINGLE
    elif not self.hasMinQuality(self.__mate): # mate doesn't have minimal mapping quality
      self.__type = Read.rtype.SINGLE
      self.__mate = None
    elif self.__read.is_duplicate: # read is duplicate
      if self.__mate.is_duplicate: # both are duplicated -> filter out
        self.__type = Read.rtype.FILTERED
      else: # only read is duplicte -> make mate single
        self.__type = Read.rtype.SINGLE
        self.__read = self.__mate
        self.__mate = None
    elif self.__mate.is_duplicate: # mate is duplicate -> make read single
      self.__type = Read.rtype.SINGLE
      self.__mate = None
    elif Read.ctype.SOFTCLIP in dict(self.__read.cigar): # read is split
      if self.isMateInverted() or mateHasGaps: # mate has gaps -> unusable
        self.__type = Read.rtype.UNUSABLE
      else: # use split read
        self.__type = Read.rtype.READ_SPLIT
    elif Read.ctype.SOFTCLIP in dict(self.__mate.cigar): # mate is split
      if self.isReadInverted() or readHasGaps: # read has gaps -> unusable
        self.__type = Read.rtype.UNUSABLE
      else: # use split read
        self.__type = Read.rtype.MATE_SPLIT
    elif not readHasGaps and not mateHasGaps: # paired without any gaps
      self.__type = Read.rtype.NORMAL
    else:
      self.__type = Read.rtype.FILTERED

  @staticmethod
  def calculateReadEnd(read):
    """
    Calculate end of read if doesn't exist
    """
    if read:
      if read.aend:
        return int(read.aend)
      else:
        return int(read.pos + (read.rlen if read.rlen else len(read.seq)))

    return 0

  @staticmethod
  def isFirst(read, mate):
    """
    Test if read is before his mate
    """
    return (read.pos <= mate.pos and read.tid == mate.tid) or read.tid < mate.tid

  @staticmethod
  def hasMinQuality(read):
    """
    Check if read has minimal quality
    """
    return read.mapq == 0 or Read.__minimalQuality <= read.mapq

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

  def readEnd(self):
    """
    Return read's end index
    """
    return self.__readEnd

  def mateEnd(self):
    """
    Return mate's end index
    """
    return self.__mateEnd

  def size(self):
    """
    Return insert size information from read
    """
    if not self.isNormal() or self.hasOverlap() or self.isRearranged() or self.isReadInverted() or self.isMateInverted() or self.isInterchromosomal():
      return 0

    return self.__read.tlen - (self.__readEnd - self.__read.pos) - (self.__mateEnd - self.__mate.pos)

  def actualSize(self):
    """
    Return counted insert size if template length is zero
    """
    size = self.size()

    if size or self.isSingle(): # size from read or read is single
      return size
    elif self.hasOverlap(): # overlapping pair
      if self.isRearranged(): # and also rearranged
        return self.__read.pos - self.__mateEnd

      return self.__mate.pos - self.__readEnd
    else: # normal or rearranged pair
      lengthBetween = 0

      if self.__read.tid != self.__mate.tid: # another chromosome
        lengthBetween = sum(self.__reflengths[self.__read.tid:self.__mate.tid])

      result = lengthBetween + self.__mate.pos - self.__readEnd

      if self.isRearranged():
        return -(result + self.__readEnd - self.__read.pos + self.__mateEnd - self.__mate.pos)

      return result

  def getReadReference(self):
    """
    Return reference of read
    """
    return self.__refnames[self.__read.tid]

  def getMateReference(self):
    """
    Return reference of read
    """
    return self.__refnames[self.__mate.tid]

  def hasOverlap(self):
    """
    Test if reads overlap
    """
    return self.__read.pos <= self.__mate.pos and self.__readEnd <= self.__mateEnd and self.__mate.pos <= self.__readEnd

  def isNormal(self):
    """
    Test if reads are normal
    """
    return self.__type == Read.rtype.NORMAL

  def isSingle(self):
    """
    Test if read is single
    """
    return self.__type == Read.rtype.SINGLE

  def isReadSplit(self):
    """
    Test if read is split
    """
    return self.__type == Read.rtype.READ_SPLIT

  def isMateSplit(self):
    """
    Test if mate is split
    """
    return self.__type == Read.rtype.MATE_SPLIT

  def isUsable(self):
    """
    Test if reads are usable for detecting variations
    """
    return self.__type != Read.rtype.UNUSABLE

  def isFiltered(self):
    """
    Test if read is filtered
    """
    return self.__type == Read.rtype.FILTERED

  def _isReadInvertedFR(self):
    """
    Test if read is inverted in forward/reverse policy
    """
    return self.__read.is_reverse

  def _isMateInvertedFR(self):
    """
    Test if read is inverted in forward/reverse policy
    """
    return not self.__mate.is_reverse

  def _isReadInvertedRF(self):
    """
    Test if read is inverted in reverse/forward policy
    """
    return not self.__read.is_reverse

  def _isMateInvertedRF(self):
    """
    Test if read is inverted in reverse/forward policy
    """
    return self.__mate.is_reverse

  def isReadInverted(self):
    """
    Test if read is inverted
    """
    return self._readInvertedRef()

  def isMateInverted(self):
    """
    Test if mate is inverted
    """
    return self._mateInvertedRef()

  def _isRearrangedFR(self):
    """
    Test if reads are rearranged in forward/reverse policy
    """
    return self.__read.is_reverse and not self.__mate.is_reverse

  def _isRearrangedRF(self):
    """
    Test if reads are rearranged in reverse/forward policy
    """
    return not self.__read.is_reverse and self.__mate.is_reverse

  def isRearranged(self):
    """
    Test if reads are rearranged
    """
    return self._rearrangedRef()

  def isInterchromosomal(self):
    """
    Test if reads are on different chromosomes
    """
    return self.__read.tid != self.__mate.tid
