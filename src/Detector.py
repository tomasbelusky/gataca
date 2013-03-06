#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from Sample import Sample
from Cluster import Cluster

class Detector:
  def __init__(self, sample, reference, policy="fr", region=None):
    self.__sample = Sample(sample, reference)
    self.__policy = policy
    self.__region = region

  def start(self):
    return 0
