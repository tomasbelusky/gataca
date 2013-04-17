#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "31.03. 2013"

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
    self._createJoinTable(self._sample)
    self.__reference = self.__variations[0].getReference()

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

  def process(self):
    """
    Find out winner variations with it's helpers and also store unused variations
      TODO check coverage
           return INS/DEL if available and can't be decided which variation is true
    """
    bestIndex = -1
    bestCount = 0
    count = 0

    for index, variations in enumerate(self.__others): # find variations with most helepers
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

    for variations in self.__others: # store unused variations
      for var in variations:
        if var not in self.__helpers and var not in self.__fromClusters:
          self.__unused.add(var)

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
