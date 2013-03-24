#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

_author_ = "Tomáš Beluský"
_date_ = "05.03. 2013"

import sys
import copy
import types
from collections import OrderedDict
from AbstractCluster import AbstractCluster
from src.Variation import Variation

class BaseCluster(AbstractCluster):
  """
  Represents cluster of base variations
  """

  def _init_(self, rname, sample):
    """
    Initialize variables
    """
    AbstractCluster._init_(self, rname, sample)

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
    fulldepth = self._sample.getExactCoverage(self._rindex, self._start, self._end)

    for var in self._alleles: # get reference sequence and common info
      refseq = var.getReferenceSequence()

      if refseq in refseqs: # same ref seq
        for key, value in refseqs[refseq].items(): # remove not common info
          if key == 'depth':
            refseqs[refseq]['depth'] = '%s,%s' % (refseqs[refseq]['depth'], var.getInfo(key))
          elif var.getInfo(key) != value:
            del refseqs[refseq][key]
      else: # new ref seq
        refseqs[refseq] = var.getInfo()

    result = []

    for refseq, info in refseqs.items(): # print row for every reference sequence
      info['fulldepth'] = fulldepth
      sequences = []

      for var in self._alleles: # get alleles with same reference sequence
        if var.getReferenceSequence() == refseq:
          sequences.append(var.getSequence())

      r = "%s\t%s\t.\t%s\t%s\t.\t.\t" % (self._rname, self.getStart(), refseq, ','.join(sequences))

      for key, value in info.items(): # print common info
        if type(value) == types.BooleanType:
          r += "%s;" % (key.upper())
        else:
          r += "%s=%s;" % (key.upper(), AbstractCluster.noInfoChars.sub('', str(value)))

      result.append(r[:-1])

    return '\n'.join(result)
