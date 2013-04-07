#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import copy

from AbstractCluster import AbstractCluster

class BaseCluster(AbstractCluster):
  """
  Represents cluster of base variations
  """

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    AbstractCluster.__init__(self, rname, sample)
    self._alleles = []

  def _process(self):
    """
    Create alleles, count start and end position
    """
    if not len(self._variations):
      return

    variations = copy.deepcopy(self._variations)
    first = variations.pop(0)
    self._alleles = [first]

    for j in range(len(variations)):
      second = variations.pop(0)
      append = True

      for i in range(len(self._alleles)): # find same variations and increment depth
        if self._alleles[i].getSequence() == second.getSequence() and self._alleles[i].getType() == second.getType():
          append = False
          self._alleles[i].incDepth()
          break

      if append: # append new variation
        self._alleles.append(second)

    self._start = first.getStart()
    self._actualStart = self._start
    self._end = self._start

    for allele in self._alleles:
      self._end = max(self._end, allele.getEnd())

  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    for var in self._variations:
      if var.allele(variation):
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
