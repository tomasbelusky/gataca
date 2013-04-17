#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import copy

from src.interface.interface import *

class Variation:
  """
  Represents finded variation in genome
  """
  vtype = enum(# type of variation
               SNP=0,
               DEL=1,
               INS=2,
               INV=3,
               DUP=4,
               DUT=5,
               TRA=6)
  svtype = (# type of variation for VCF output
            "DEL",
            "INS",
            "INV",
            "DUP",
            "DUP:TANDEM",
            "INS:TRA")
  mtype = enum(# type of method which find variation
               CIGAR_MD=0,
               READ_PAIR=1,
               SPLIT_READ=2,
               JOINED=3)

  def __init__(self, vtype, reference, start, seq, refseq, mtype, info={}):
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
    self.__info = info

    if 'depth' not in self.__info:
      self.__info['depth'] = 1

    if self.__type != Variation.vtype.SNP: # set SVTYPE
      self.__info['svtype'] = Variation.svtype[self.__type - 1]

      if not self.__seq: # create symbolic sequence for alternate
        if self.__type == Variation.vtype.INV:
          self.__seq = "<%s>" % self.__info['svtype']
        else:
          self.__seq = "%s<%s>" % (self.__refseq, self.__info['svtype'])

  def contain(self, var):
    """
    Test if variation contains another variation
    """
    return self.getStart() <= var.getStart() and var.getEnd() <= self.getEnd()

  def overlap(self, var):
    """
    Test if variation overlap with another variation on right side
    """
    return (self.getStart() <= var.getStart() and var.getStart() <= self.getEnd()) or \
      (var.getStart() <= self.getStart() and self.getStart() <= var.getEnd())

  def allele(self, var):
    """
    Test if another variation represents another allele of actual variation
    """
    return self.getStart() == var.getStart() and not self.isImprecise() and not var.isImprecise()

  def getType(self):
    """
    Return type of variation
    """
    return self.__type

  def getSvtype(self):
    """
    Return type of SV (usefull for debugging)
    """
    return self.__info.get('svtype', 'SNP')

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
    Return rightmost index of variation if exist, else return leftmost index
    """
    return self.__info.get('end', self.__start)

  def getMaxStart(self):
    """
    Return max lefmost index of variation
    """
    return self.__start + self.__info.get('cpos', 0)

  def getMaxEnd(self):
    """
    Return max rightmost index of variation if exist, else return leftmost index
    """
    return self.__info.get('end', self.__start) + self.__info.get('cend', 0)

  def isImprecise(self):
    """
    Test if variation is precise
    """
    return self.__info.get('imprecise', False)

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
      return copy.deepcopy(self.__info)

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
