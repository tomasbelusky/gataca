#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from interface import *

class Variation:
  varType = enum(SNP=0,
                 DELETION=1,
                 INSERTION=2,
                 INVERSION=3,
                 DUPLICATION_TANDEM=4,
                 DUPLICATION_INTERSPERSED=5,
                 TRANSLOCATION=6)
  methodType = enum(READ_PAIR=0,
                    SPLIT_READ=1)

  def __init__(self, type, start, end, sequence, method, *args):
    self.__type = type
    self.__method = method
    self.__start = start
    self.__end = end
    self.__clusters = []
    self.__reads = None

    for arg in args: # save reads
      self.__read.append(arg)
