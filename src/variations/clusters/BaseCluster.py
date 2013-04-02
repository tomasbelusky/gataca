#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import copy
import types

from variations.Variation import Variation

class BaseCluster:
  """
  Represents cluster of base variations
  """
  noInfoChars = re.compile(r'[\[\]\s\']') # characters which wouldn't be printed in VCF output

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    self._rname = rname
    self._rindex = sample.getRefIndex(rname)
    self._sample = sample
    self._variations = []
    self._start = 0
    self._end = 0
    self._alleles = []

  @staticmethod
  def joinSame(first, second):
    """
    Join variations with same start and type
    INFO: Not completed method and maybe will be deleted
    """
    var = []

    if first.getSequence().startswith(second.getSequence()):
      var.append(first)
    elif second.getSequence().startswith(first.getSequence()):
      var.append(second)
    else:
      var.append(first)
      var.append(second)

    return var

  @staticmethod
  def joinContain(first, second):
    """
    Join variations where first contains second
    INFO: Not completed method and maybe will be deleted
    """
    var = []

    # remove SNP or INV in DEL or INS
    if first.getType() in (Variation.vtype.DEL, Variation.vtype.INS) and second.getType() in (Variation.vtype.SNP, Variation.vtype.INV):
      pass

    return var

  @staticmethod
  def joinOverlap(first, second):
    """
    Join variations that overlap
    INFO: Not completed method and maybe will be deleted
    """
    var = []
    return var

  @staticmethod
  def join(first, second):
    """
    Join variations
    INFO: Not completed method and maybe will be deleted
    """
    var = []

    if first.getType() == second.getType() and first.getStart() == second.getStart(): # same type and start
      var = BaseCluster.joinSame(first, second)
    elif first.overlap(second): # overlap (first on the left)
      var = BaseCluster.joinOverlap(first, second)
    elif second.overlap(first): # overlap (second on the left)
      var = BaseCluster.joinOverlap(second, first)
    elif first.contain(second): # first contain second
      var = BaseCluster.joinContain(first, second)
    elif second.contain(first): # second contain first
      var = BaseCluster.joinContain(second, first)

    return var

  @staticmethod
  def getAlleles(variations):
    """
    Return leftmost index of alleles and list of alleles
    """
    if not len(variations):
      return 0, []

    first = variations.pop(0)
    alleles = [first]

    for j in range(len(variations)):
      second = variations.pop(0)
      append = True

      for i in range(len(alleles)): # find same variations and increment depth
        if alleles[i].getSequence() == second.getSequence() and alleles[i].getType() == second.getType():
          append = False
          alleles[i].incDepth()
          break

      if append: # append new variation
        alleles.append(second)

    return first.getStart(), alleles

  @staticmethod
  def infoString(info):
    """
    Create string from info fields
    """
    if not info:
      return ""

    result = ""

    for key, value in info.items(): # print common info
      if type(value) == types.BooleanType:
        result += "%s;" % (key.upper())
      else:
        result += "%s=%s;" % (key.upper(), BaseCluster.noInfoChars.sub('', str(value)))

    return result[:-1]

  def getStart(self):
    """
    Return index of leftmost base from variations in cluster
    """
    return self._start

  def getEnd(self):
    """
    Return index of rightmost base from variations in cluster
    """
    return self._end

  def add(self, variation):
    """
    Add variation into cluster
    """
    self._variations.append(variation)
    self._process()

  def remove(self, variation):
    """
    Remove variation from cluster
    """
    if variation in self.__variations:
      self._variations.remove(variation)
      self._process()

  def _process(self):
    """
    Join variations together for printing them in VCF output
    """
    (self._start, self._alleles) = self.getAlleles(copy.deepcopy(self._variations))
    self._end = self._start

    for allele in self._alleles:
      self._end = max(self._end, allele.getEnd())

    """
    INFO: use with join* methods
    for start in alleles: # lists of alleles excluding first list
      for i in range(len(result)): # alleles from first list
        for allele in alleles[start]: # alleles from second list
          #self._alleles.extend(self.join(allele, result.pop(0)))
          self._alleles.extend([allele, result.pop(0)])

    for allele in self._alleles:
      self._end = max(self._end, allele.getEnd())
    """

  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    for var in self._variations:
      if var.allele(variation):# or len(self.join(var, variation)) == 1:
        return True

    return False

  def __str__(self):
    """
    Print cluster in VCF format
    """
    if not len(self._alleles):
      return ""

    refseqs = {}
    fulldepth = self._sample.getExactCoverages(self._rindex, self._start, self._end)

    for var in self._alleles: # get reference sequence and common info
      refseq = var.getReferenceSequence()

      if refseq in refseqs: # reference sequence exists
        for key, value in refseqs[refseq].items():
          if key == 'depth': # add depth of allele
            refseqs[refseq]['depth'].append(var.getInfo(key))
          elif var.getInfo(key) != value: # remove not common info
            del refseqs[refseq][key]
      else: # new reference sequence
        refseqs[refseq] = var.getInfo()
        refseqs[refseq]['depth'] = [refseqs[refseq]['depth']]

    result = []

    for refseq, info in refseqs.items(): # print row for every reference sequence
      info['fulldepth'] = fulldepth
      sequences = []

      for var in self._alleles: # get alleles with same reference sequence
        if var.getReferenceSequence() == refseq:
          sequences.append(var.getSequence())

      result.append("%s\t%s\t.\t%s\t%s\t.\t.\t%s" % (self._rname,
                                                       self._start,
                                                       refseq,
                                                       ','.join(sequences),
                                                       self.infoString(info)))

    return '\n'.join(result)
