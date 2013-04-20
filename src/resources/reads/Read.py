#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "11.04. 2013"

from src.interface.Settings import Settings
from src.interface.interface import *
from SplitPart import SplitPart
from Cigar import Cigar

class Read:
  """
  Represents sequenced read
  """
  COUNT_SEPARATOR = "_"
  ptype = enum(# policy type
               FR=0, # forward/reverse
               RF=1) # reverse/forward

  def __init__(self, read, first, strand, refnames):
    """
    Initialize variables
    """
    self.__read = read
    self.__first = first
    self.__strand = strand
    self.__refnames = refnames

    if not self.isUnmapped():
      self.__end = self.calculateEnd(self.__read)
      self.__len = self.__end - self.__read.pos
      self.__reference = self.__refnames[self.__read.tid]

  @staticmethod
  def calculateEnd(read):
    """
    Calculate end of read if doesn't exist
    """
    if read:
      if read.aend:
        return int(read.aend)
      else:
        return int(read.pos + (read.rlen if read.rlen else len(read.seq)))

    return 0

  @property
  def sam(self):
    """
    Return SAM read object
    """
    return self.__read

  @property
  def tid(self):
    """
    Return id of reference
    """
    return self.__read.tid

  @property
  def pos(self):
    """
    Return start position
    """
    return self.__read.pos

  @property
  def end(self):
    """1
    Return end position
    """
    return self.__end

  @property
  def len(self):
    """
    Return length
    """
    return self.__len

  @property
  def reference(self):
    """
    Return reference
    """
    return self.__reference

  def isUnmapped(self):
    """
    Test if read is unmapped
    """
    return self.__read is None

  def isDuplicate(self):
    """
    Test if read is duplicate
    """
    return self.__read.is_duplicate

  def isInverted(self):
    """
    Test if read is inverted
    """
    return (not self.__read.is_reverse) != self.__strand

  def isSplit(self):
    """
    Test if read is split
    """
    return Cigar.op.SOFTCLIP in dict(self.__read.cigar)

  def hasGaps(self):
    """
    Test if read has gaps
    """
    gaps = (Cigar.op.SKIPPED, Cigar.op.SOFTCLIP, Cigar.op.HARDCLIP, Cigar.op.PADDING)
    return set(dict(self.__read.cigar)).intersection(gaps)

  def hasMinQuality(self):
    """
    Check if read has minimal quality
    """
    return self.__read.mapq == 0 or Settings.MIN_QUALITY <= self.__read.mapq

  def getMappedParts(self):
    """
    Return array of mapped parts
    """
    pos = self.__read.pos
    start = 0
    end = 0
    parts = []
    actualCigar = []

    for operator, length in self.__read.cigar:
      if operator == Cigar.op.SOFTCLIP:
        if actualCigar: # append part
          parts.append(SplitPart(pos,
                                 self.__read.seq[start:end],
                                 actualCigar,
                                 self.__read.mapq,
                                 self.__read.is_reverse,
                                 False,
                                 False))
          pos += length
          actualCigar = []

        start += length
      elif operator in Cigar.sums:
        end = start + length
        actualCigar.append((operator, length))

    if actualCigar:
      parts.append(SplitPart(pos,
                             self.__read.seq[start:end],
                             actualCigar,
                             self.__read.mapq,
                             self.__read.is_reverse,
                             False,
                             False))

    return parts

  def getFastqSplits(self):
    """
    Return read's split parts in FASTQ format
    """
    count = 0
    start = 0
    end = 0
    value = ''

    for operator, length in self.__read.cigar:
      if operator != Cigar.op.SOFTCLIP:
        start += length
      else:
        end = start + length
        value += '@%s%s%d/%d\n%s\n+\n%s\n' % (self.__read.qname,
                                              Read.COUNT_SEPARATOR,
                                              count,
                                              1 if self.__first else 2,
                                              self.__read.seq[start:end],
                                              self.__read.qual[start:end])

      count += 1

    return value

  def getFastq(self):
    """
    Return read in FASTQ format
    """
    return '@%s/%d\n+\n%s\n%s\n' % (self.__read.qname, 1 if self.__first else 2, self.__read.seq, self.__read.qual)

  def getSam(self):
    """
    Return read in SAM format
    """
    value = '%s\t%d\t%s\t%d\t%d\t' % (self.__read.qname, self.__read.flag, self.__reference, self.__read.pos + 1, self.__read.mapq)

    if self.__read.cigar is None:
      value += '*'
    else:
      for operator, length in self.__read.cigar:
        value += '%d%s' % (length, Cigar.abbr[operator])

    value += '\t%s\t%d\t%d\t%s\t%s' % (self.__refnames[self.__read.rnext], self.__read.pnext + 1, self.__read.tlen, self.__read.seq, self.__read.qual)

    for tagname, tagvalue in self.__read.tags:
      if isinstance(tagvalue, int):
        t = 'i'
      elif isinstance(tagvalue, float):
        t = 'f'
      elif isinstance(tagvalue, str):
        if len(tagvalue) == 1:
          t = 'A'
        else:
          t = 'Z'
      else:
        t = 'B'

      value += '\t%s:%s:%s' % (tagname, t, tagvalue)

    return '%s\n' % value
