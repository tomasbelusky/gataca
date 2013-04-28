#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "31.03. 2013"

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

  def __coverageProcess(self, bestCount):
    """
    Find out winner variation with coverage
    """
    winnerIndex = -1

    for index, variations in enumerate(self.__others):
      if len(variations) == bestCount:
        actualWinnerIndex = -1
        variation = self.__variations[index]
        coverage = self._sample.getInexactCoverage(variation.getReference(), variation.getStart(), variation.getEnd())

        if variation.getType() == Variation.vtype.DEL: # deletion
          if coverage < self._sample.getMinCoverage():
            actualWinnerIndex = index
        elif variation.getType() == Variation.vtype.DUP: # duplication
          if self._sample.getMaxCoverage() < coverage:
            actualWinnerIndex = index
        elif self._sample.getMinCoverage() <= coverage and coverage <= self._sample.getMaxCoverage(): # other variation
          actualWinnerIndex = index

        if actualWinnerIndex != -1: # founded winner
          if winnerIndex == -1:
            winnerIndex = actualWinnerIndex
          else: # can't have more winners
            return

    if winnerIndex != -1: # winner exists
      self.__winner = self.__variations[winnerIndex]
      self.__helpers = set(self.__others.pop(winnerIndex))

  def __findOutWinner(self):
    """
    Find out winner in variations
    """
    bestIndex = -1
    bestCount = 0
    count = 0

    for index, variations in enumerate(self.__others): # find variations with most helpers
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

    return bestCount

  def process(self):
    """
    Do some processing in cluster
    """
    bestCount = self.__findOutWinner()

    if self.__winner is None: # find out winner with coverage
      self.__coverageProcess(bestCount)

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
