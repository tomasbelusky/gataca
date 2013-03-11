#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from Variation import Variation

class Cluster:
  def __init__(self, vtype, variation):
    self.__type = vtype
    self.__variations = [variation]
    self.__start = variation.getStart()
    self.__end = variation.getEnd()
    self.__info = {}

  def add(self, variation):
    self.__variations.append(variation)

  def remove(self, variation):
    if variation in self.__variations:
      self.__variations.remove(variation)

  def compare(variation):
    pass

  def __str__(self):
    var = self.__variations[0]
    result = "%s\t%s\t.\t%s\t%s\t.\t.\t" % (var.getReference(), var.getStart(), var.getReferenceSequence(), var.getSequence())
    first = True

    for key, value in self.__info:
      if not first: # separate values
        output.write(";")

      result += "%s=%s" % (key, value)

    return result
