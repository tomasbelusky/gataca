#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "31.03. 2013"

class OppositeCluster:
  """
  Represents cluster with variations where only one variation is true
  """

  def __init__(self, *args):
    self.__variations = args

  def process(self, variation):
    return False

  def extend(self, cluster):
    return False

  def getWinner(self):
    return None

  def getUnused(self):
    return []
