#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "04.04. 2013"

import copy
import operator

from variations.clusters.AbstractCluster import AbstractCluster
from variations.factories.JoinFactory import JoinFactory
from variations.Variation import Variation

class StructuralCluster(AbstractCluster):
  """
  Represents cluster of structural variations
  """

  def __init__(self, rname, sample, variation):
    """
    Initialize variables
    """
    AbstractCluster.__init__(self, rname, sample)
    self.__consensus = variation
    self.__processConsensus()
    self._createJoinTable(sample)

  def __processConsensus(self):
    """
    Get informations from consensus sequence
    """
    if self.__consensus:
      self.__type = self.__consensus.getType()
      self._reference = self.__consensus.getReference()
      self._start = self.__consensus.getMaxStart()
      self._end = self.__consensus.getMaxEnd()
      self._actualStart = self.__consensus.getStart()
      self._actualEnd = self.__consensus.getEnd()

  def _createJoinTable(self, sample):
    """
    Create table of joining variations
    """
    self._joinFactory = JoinFactory(sample)
    self._joinFactoryRef = {}
    self._joinFactoryRef[Variation.vtype.INS] = {Variation.vtype.INS : self._joinFactory.insertion,
                                                 Variation.vtype.TRA : self._joinFactory.translocationInsertion}
    self._joinFactoryRef[Variation.vtype.DEL] = {Variation.vtype.DEL : self._joinFactory.deletion}
    self._joinFactoryRef[Variation.vtype.INV] = {Variation.vtype.INV : self._joinFactory.inversion}
    self._joinFactoryRef[Variation.vtype.DUP] = {Variation.vtype.DUP : self._joinFactory.duplication,
                                                 Variation.vtype.DUT : self._joinFactory.tandemDuplication}
    self._joinFactoryRef[Variation.vtype.DUT] = {Variation.vtype.DUT : self._joinFactory.tandemDuplication,
                                                 Variation.vtype.DUP : self._joinFactory.tandemDuplication}
    self._joinFactoryRef[Variation.vtype.TRA] = {Variation.vtype.TRA : self._joinFactory.translocation,
                                                 Variation.vtype.INS : self._joinFactory.translocationInsertion}

  def _join(self, first, second):
    """
    Join two variations into one
    """
    if second.getStart() < first.getStart():
      firstCopy = copy.deepcopy(second)
      secondCopy = copy.deepcopy(first)
    else:
      firstCopy = copy.deepcopy(first)
      secondCopy = copy.deepcopy(second)

    return self._joinFactoryRef[first.getType()][second.getType()](firstCopy, secondCopy)

  def add(self, variation):
    """
    Try to add variation into cluster and return if it fits into cluster
    """
    if not self._joinFactoryRef[self.__type].get(variation.getType(), None):
      return False

    newConsensus = self._join(self.__consensus, variation)

    if newConsensus is not None:
      self.__consensus = newConsensus
      self.__processConsensus()
      return True

    return False

  def __str__(self):
    """
    Print cluster in VCF format
    """
    if not self.__consensus:
      return ""

    info = self.__consensus.getInfo()
    info['intervals'] = sorted(info['intervals'], key=operator.itemgetter(0))
    depth = len(info['intervals'])
    fulldepth = 0
    intervals = {}

    while info['intervals']: # get overlaped intervals
      int1 = info['intervals'].pop(0)
      count = 1

      for int2 in info['intervals']:
        if (int1[0] <= int2[0] and int2[0] <= int1[1]) or (int2[0] <= int1[0] and int1[0] <= int2[1]):
          breakends = zip(int1, int2)
          int1 = [min(breakends[0]), max(breakends[1])]
          count += 1

      intervals[tuple(int1)] = count

    """
    INFO: don't count in test phase
    for start, end in intervals: # get fulldepth in intervals
      fulldepth += self._sample.getExactCoverages(self._rindex, start, end)
    """

    info['conf'] = self.countConfidence(depth, fulldepth)
    return "%s\t%s\t.\t%s\t%s\t.\t.\t%s" % (self._rname,
                                            self._actualStart,
                                            self.__consensus.getReferenceSequence(),
                                            self.__consensus.getSequence(),
                                            self.infoString(info))
