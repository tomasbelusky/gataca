#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "12.04. 2013"

from Cigar import Cigar
from src.interface.Settings import Settings

class SplitPart:
  """
  Represents part of split read
  """

  def __init__(self, pos, seq, cigar, mapq, reverse, unmapped, remapped):
    """
    Initialize variables
    """
    self.__pos = pos
    self.__seq = seq
    self.__cigar = cigar
    self.__mapq = mapq
    self.__reverse = reverse
    self.__unmapped = unmapped
    self.__remapped = remapped
    self.__strand = True
    self.__len = len(self.__seq)
    self.__end = self.__pos + self.__len

    if self.__cigar and self.__cigar[0][0] == Cigar.op.SOFTCLIP:
      self.__pos -= self.__cigar[0][1]

  def hasMinQuality(self):
    """
    Test if part has minimal mapping quality
    """
    return self.__mapq == 0 or Settings.MIN_QUALITY <= self.__mapq

  @property
  def pos(self):
    """
    Return start position
    """
    return self.__pos

  @property
  def end(self):
    """
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
  def seq(self):
    """
    Return sequence
    """
    return self.__seq

  @property
  def cigar(self):
    """
    Return cigar
    """
    return self.__cigar

  def isInverted(self):
    """
    Test if part is inverted
    """
    return (not self.__reverse) != self.__strand

  def isUnmapped(self):
    """
    test if part is unmapped
    """
    return self.__unmapped

  def isRemapped(self):
    """
    Test if part is remapped
    """
    return self.__remapped

  def setStrand(self, strand):
    """
    Set strand where part should be mapped
    """
    self.__strand = strand
