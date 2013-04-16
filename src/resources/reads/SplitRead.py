#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "12.04. 2013"

from interface.interface import *
from Cigar import Cigar

class SplitRead:
  """
  Represents read with two split parts
  """
  MIN_PART_LENGTH = 10

  def __init__(self, original, strand, parts):
    """
    Initialize variables
    """
    self.__original = original
    self.__strand = strand
    self.__parts = parts
    self.__usable = len(self.__parts) == 2
    sortedParts = sorted(self.__parts, key=lambda p: p.pos)
    self.__rearranged = self.__parts != sortedParts
    self.__primary = None
    self.__remapped = None
    self.__left = None
    self.__right = None

    if self.__usable:
      self.__left = sortedParts[0]
      self.__left.setStrand(self.__strand)
      self.__right = sortedParts[1]
      self.__right.setStrand(self.__strand)
      self.__remapped = self.__right if self.__right.isRemapped() else self.__left
      self.__primary = self.__left if self.__right.isRemapped() else self.__right

  @property
  def original(self):
    """
    Return original read
    """
    return self.__original

  @property
  def left(self):
    """
    Return left part of read
    """
    return self.__left

  @property
  def right(self):
    """
    Return right part of read
    """
    return self.__right

  @property
  def primary(self):
    """
    Return primary part
    """
    return self.__primary

  @property
  def remapped(self):
    """
    Return remapped part
    """
    return self.__remapped

  def hasMinQuality(self):
    """
    Test if whole read has minimal quality
    """
    return not len([part for part in self.__parts if not part.hasMinQuality()])

  def hasMinLengths(self):
    """
    Test if parts have minimal lengths
    """
    fullLength = 0

    for operator, length in self.__original.sam.cigar:
      if operator == Cigar.op.SOFTCLIP:
        if fullLength:
          if fullLength < SplitRead.MIN_PART_LENGTH:
            return False

        fullLength = 0

        if length < SplitRead.MIN_PART_LENGTH:
          return False
      elif operator in Cigar.sums:
        fullLength += length

    return fullLength == 0 or SplitRead.MIN_PART_LENGTH <= fullLength

  def hasInsertion(self):
    """
    Test if read span insertion
    """
    dLeft = dict(self.__left.cigar)
    dRight = dict(self.__right.cigar)
    return Cigar.op.SOFTCLIP in dLeft or Cigar.op.HARDCLIP in dLeft or \
           Cigar.op.SOFTCLIP in dRight or Cigar.op.HARDCLIP in dRight

  def hasGap(self):
    """
    Test if thete is gap between parts
    """
    return (self.__left.end + 1) < self.__right.pos

  def hasOverlap(self):
    """
    Test if parts overlap
    """
    return self.__left.pos <= self.__right.pos and \
            ((self.__left.end <= self.__right.end and \
              self.__right.pos <= self.__left.end) or \
             (self.__right.end <= self.__left.end))

  def isUsable(self):
    """
    Test if parts are usable
    """
    return self.__usable

  def isRearranged(self):
    """
    Rest if parts are rearranged
    """
    return self.__rearranged
