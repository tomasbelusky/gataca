#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "31.03. 2013"

from ImpreciseCluster import ImpreciseCluster

class OppositeCluster(ImpreciseCluster):
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

  def helpDecide(self, variation):
    """
    Help cluster to decide which variation is true
    """
    value = False

    for index, var in enumerate(self.__variations): # try to join with all variations
      if self._joinFactoryRef[var.getType()].get(variation.getType(), None):
        newVariation = self._join(var, variation)

        if newVariation: # joined
          self.__variations[index] = newVariation
          self.__others[index].append(variation)
          value = True

    return value

  def extend(self, cluster):
    """
    Extend actual cluster by another cluster
    """
    value = False
    tmp = set()

    for variation in cluster.__variations: # try to join with all variations
      tmp.add(variation)

      if self.helpDecide(variation):
        value = True

    if value:
      self.__fromClusters |= tmp

    return value

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
