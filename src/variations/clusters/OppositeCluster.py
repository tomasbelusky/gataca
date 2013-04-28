#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "31.03. 2013"

import sys

from src.variations.Variation import Variation
from StructuralCluster import StructuralCluster

class OppositeCluster(StructuralCluster):
  """
  Represents cluster with variations where only one variation is true
  """

  def __init__(self, sample, *args):
    """
    Initialize variables
    """
    self.__variations = list(args)
    self._sample = sample
    self.__unused = set()
    self.__helpers = set()
    self.__fromClusters = set()
    self.__winner = None
    self.__others = [[] for v in self.__variations]
    self.__reference = self.__variations[0].getReference()
    self._createJoinTable(self._sample)

  def __coverageProcess(self):
    """
    Find out winner variation with coverage
    """
    delWinnerIndex = -1
    delWinnerCoverage = 0
    dupWinnerIndex = -1
    dupWinnerCoverage = sys.maxint

    for index, variation in enumerate(self.__variations):
      if variation.getType() in (Variation.vtype.DEL, Variation.vtype.DUP):
        coverage = self._sample.getInexactCoverage(variation.getReference(), variation.getStart(), variation.getEnd())

        if variation.getType() == Variation.vtype.DEL: # deletion
          if coverage < self._sample.getMinCoverage() or coverage < delWinnerCoverage:
            delWinnerIndex = index
        elif variation.getType() == Variation.vtype.DUP: # duplication
          if self._sample.getMaxCoverage() < coverage or dupWinnerCoverage < coverage:
            dupWinnerIndex = index

    if delWinnerIndex != -1: # deletion
      if dupWinnerIndex == -1: # duplication can't have winner
        self.__winner = self.__variations[delWinnerIndex]
        self.__helpers = set(self.__others.pop(delWinnerIndex))
    elif dupWinnerIndex != -1: # duplication
      self.__winner = self.__variations[dupWinnerIndex]
      self.__helpers = set(self.__others.pop(dupWinnerIndex))
    else: # other variation than deletion and duplication
      self.__findOutWinner(True)

  def __findOutWinner(self, skipCNV=False):
    """
    Find out winner in variations
    """
    bestIndex = -1
    bestCount = 0
    count = 0

    for index, variations in enumerate(self.__others): # find variations with most helpers
      if skipCNV and self.__variations[index].getType in (Variation.vtype.DEL, Variation.vtype.DUP):
        continue

      actualCount = len(variations)

      if actualCount > bestCount: # more helpers
        bestIndex = index
        bestCount = len(self.__others[index])
        count = 1
      elif actualCount == bestCount: # same count of helpers
        count += 1

    if count == 1: # winner exists
      self.__winner = self.__variations[bestIndex]
      self.__helpers = set(self.__others.pop(bestIndex))

  def process(self):
    """
    Do some processing in cluster
    """
    self.__findOutWinner()

    if self.__winner is None: # find out winner with coverage
      self.__coverageProcess()

    for variations in self.__others: # store unused variations
      for var in variations:
        if var not in self.__helpers and var not in self.__fromClusters:
          self.__unused.add(var)

  def getVariations(self):
    """
    Return variations
    """
    return self.__variations

  def helpDecide(self, index, variation, fromCluster):
    """
    Another variationb try to help decide with variation is true
    """
    newVariation = None

    if self._joinFactoryRef[self.__variations[index].getType()].get(variation.getType(), None):
      newVariation = self._join(self.__variations[index], variation)

      if newVariation: # joined
        if fromCluster:
          self.__fromClusters.add(variation)

        self.__variations[index] = newVariation
        self.__others[index].append(variation)

    return newVariation

  def getWinner(self):
    """
    Return winner variation
    """
    return self.__winner

  def getUnused(self):
    """
    Return unused variations
    """
    return self.__unused

  def getHelpers(self):
    """
    Return variation which helped to win
    """
    return self.__helpers
