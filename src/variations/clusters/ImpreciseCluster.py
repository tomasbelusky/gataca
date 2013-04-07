#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "04.04. 2013"

import copy

from variations.clusters.AbstractCluster import AbstractCluster
from variations.factories.JoinFactory import JoinFactory
from variations.Variation import Variation

class ImpreciseCluster(AbstractCluster):
  """
  Represents cluster of imprecise variations
  """

  def __init__(self, rname, sample):
    """
    Initialize variables
    """
    AbstractCluster.__init__(self, rname, sample)
    self._consensus = None
    self.__type = None
    self._joinFactory = JoinFactory(sample)
    self.__joinFactoryRef = {Variation.vtype.INS : {Variation.vtype.INS : self._joinFactory.insertion,
                                                    Variation.vtype.TRA : self._joinFactory.insertionTranslocation},
                             Variation.vtype.DEL : {Variation.vtype.DEL: self._joinFactory.deletion},
                             Variation.vtype.INV : {Variation.vtype.INV: self._joinFactory.inversion},
                             Variation.vtype.DUP : {Variation.vtype.DUP: self._joinFactory.duplication},
                             Variation.vtype.TRA : {Variation.vtype.TRA: self._joinFactory.translocation,
                                                    Variation.vtype.INS : self._joinFactory.insertionTranslocation}}

  def join(self, first, second):
    """
    Join two variations into one
    """
    if second.getStart() < first.getStart():
      firstCopy = copy.deepcopy(second)
      secondCopy = copy.deepcopy(first)
    else:
      firstCopy = copy.deepcopy(first)
      secondCopy = copy.deepcopy(second)

    return self.__joinFactoryRef[first.getType()][second.getType()](firstCopy, secondCopy)

  def _process(self):
    """
    Join variations together, count start and end position
    """
    if not len(self._variations):
      return

    variations = sorted(copy.deepcopy(self._variations), key=lambda v: v.getStart())

    while len(variations) > 1: # join each variation together
      first = variations.pop(0)
      second = variations.pop(0)
      variations.insert(0, self.join(first, second))

    self._consensus = variations.pop()
    self.__type = self._consensus.getType()
    self._start = self._consensus.getMaxStart()
    self._end = self._consensus.getMaxEnd()
    self._actualStart = self._consensus.getStart()

  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    if not self.__joinFactoryRef[self.__type].get(variation.getType(), None):
      return False

    for var in self._variations:
      if self.join(var, variation) is not None:
        return True

    return False

  def __str__(self):
    """
    Print cluster in VCF format
    """
    if not self._consensus:
      return ""

    info = self._consensus.getInfo()
    info['fulldepth'] = self._sample.getExactCoverages(self._rindex, self._start, self._end)
    return "%s\t%s\t.\t%s\t%s\t.\t.\t%s" % (self._rname,
                                            self._actualStart,
                                            self._consensus.getReferenceSequence(),
                                            self._consensus.getSequence(),
                                            self.infoString(info))
