#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from Variation import Variation

class Cluster:
  def __init__(self, type, start, end):
    self.__start = start
    self.__end = end
