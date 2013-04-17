#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from AbstractCluster import AbstractCluster
from src.variations.Variation import Variation

class SnpCluster(AbstractCluster):
  """
  Represents cluster of snp variations
  """

  def __init__(self, rname, sample, variation):
    """
    Initialize variables
    """
    AbstractCluster.__init__(self, rname, sample)
    self.__alleles = [variation]
    self._reference = variation.getReference()
    self._start = variation.getStart()
    self._actualStart = self._start
    self._end = variation.getEnd()

  def add(self, variation):
    """
    Try to add variation into cluster and return if it fits into cluster
    """
    if variation.getType() == Variation.vtype.SNP and variation.getStart() == self._start:
      append = True

      for i in range(len(self.__alleles)): # find same variations and increment depth
        if self.__alleles[i].getSequence() == variation.getSequence() and self.__alleles[i].getType() == variation.getType():
          append = False
          self.__alleles[i].incDepth()
          break

      if append: # append new variation
        self.__alleles.append(variation)
        self._end = max(self._end, variation.getEnd())

      return True

    return False

  def __str__(self):
    """
    Print cluster in VCF format
    """
    if not len(self.__alleles):
      return ""

    refseqs = {}
    fulldepth = self._sample.getExactCoverages(self._rindex, self._start, self._end)

    for var in self.__alleles: # get reference sequence and common info
      refseq = var.getReferenceSequence()

      if refseq in refseqs: # reference sequence exists
        for key, value in refseqs[refseq].items():
          if key == 'depth': # add depth of allele
            refseqs[refseq]['conf'].append(self.countConfidence(var.getInfo('depth'), fulldepth))
          elif var.getInfo(key) != value: # remove not common info
            del refseqs[refseq][key]
      else: # new reference sequence
        refseqs[refseq] = var.getInfo()
        refseqs[refseq]['conf'] = [self.countConfidence(refseqs[refseq]['depth'], fulldepth)]

    result = []

    for refseq, info in refseqs.items(): # print row for every reference sequence
      sequences = []

      for var in self.__alleles: # get alleles with same reference sequence
        if var.getReferenceSequence() == refseq:
          sequences.append(var.getSequence())

      result.append("%s\t%s\t.\t%s\t%s\t.\t.\t%s" % (self._rname,
                                                     self._start,
                                                     refseq,
                                                     ','.join(sequences),
                                                     self.infoString(info)))

    return '\n'.join(result)
