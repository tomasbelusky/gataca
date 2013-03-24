#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from interface import *

class Variation:
  """
  Represents finded variation in genome
  """
  vtype = enum( # type of variation
    SNP=0,
    DEL=1,
    INS=2,
    INV=3,
    DUP=4,
    DUT=5,
    TRA=6)
  svtype = ( # type of variation for VCF output
    "DEL",
    "INS",
    "INV",
    "BND")
  mtype = enum( # type of method which find variation
    CIGAR_MD=0,
    READ_PAIR=1,
    SPLIT_READ=2,
    JOINED=3)

  def __init__(self, vtype, reference, start, seq, refseq, mtype, *args, **kwargs):
    """
    Initialize variables
    """
    self.__type = vtype
    self.__reference = reference
    self.__method = mtype
    self.__start = start
    self.__seq = seq
    self.__refseq = refseq
    self.__clusters = []
    self.__reads = args
    self.__info = kwargs
    self.__info['depth'] = 1

    if vtype in (Variation.vtype.DEL, Variation.vtype.INV, Variation.vtype.INS):
      self.__info['svtype'] = Variation.svtype[vtype - 1]
    elif vtype in (Variation.vtype.DUP, Variation.vtype.DUT, Variation.vtype.TRA):
      self.__info['svtype'] = "BND"

  @staticmethod
  def joinInfo(info1, info2):
    """
    Join and return informations from two variations
    TODO: joining different informations
    """
    info = {}

    for key, value, in info1.items():
      if key in info2: # same key in info2
        if value == info2[key]: # same value
          info[key] = value
        else: # join different info
         pass
      else: # only in first info
        info[key] = value

    for key, value, in info2.items(): # aad info which are only in info2
      if key not in info1: # unique info
        info[key] = value

    return info

  def contain(self, var):
    """
    Test if variation contains another variation
    """
    return self.getStart() <= var.getStart() and var.getEnd() <= self.getEnd()

  def overlap(self, var):
    """
    Test if variation overlap with another variation on right side
    """
    return self.getStart() < var.getStart() and self.getEnd() < var.getEnd() and var.getStart() <= self.getEnd()

  def allele(self, var):
    """
    Test if another variation represents another allele of actual variation
    """
    return self.getStart() == var.getStart() and not self.getInfo('imprecise') and not var.getInfo('imprecise')

  def getType(self):
    """
    Return type of variation
    """
    return self.__type

  def getMethod(self):
    """
    Return method which find variation
    """
    return self.__method

  def getReference(self):
    """
    Return reference name
    """
    return self.__reference

  def getStart(self):
    """
    Return lefmost index of variation
    """
    return self.__start

  def getEnd(self):
    """
    Return right most index of variation if exist, else return leftmost index
    """
    if 'end' in self.__info:
      return self.__info['end']
    else:
      return self.__start

  def getSequence(self):
    """
    Return variation sequence
    """
    return self.__seq

  def getReferenceSequence(self):
    """
    Return reference sequence
    """
    return self.__refseq

  def getInfo(self, key=None):
    """
    Return all informations or only one information if key is not None
    """
    if key:
      return self.__info.get(key, None)
    else:
      return self.__info

  def incDepth(self):
    """
    Increment depth
    """
    self.__info['depth'] += 1

  def decDepth(self):
    """
    Decrement depth
    """
    self.__info['depth'] -= 1
