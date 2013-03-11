#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

from interface import *

class Variation:
  vtype = enum(SNP=0,
               DEL=1,
               INS=2,
               INV=3,
               DUP=4,
               DUT=5,
               TRA=6)
  mtype = enum(CIGAR_MD=0,
               READ_PAIR=1,
               SPLIT_READ=2)

  def __init__(self, vtype, reference, start, end, seq, mtype, *args, **kwargs):
    self.__type = vtype
    self.__reference = reference
    self.__method = mtype
    self.__start = start
    self.__end = end
    self.__seq = seq
    self.__refseq = kwargs.get("refseq", "")
    self.__len = 0
    self.__clusters = []
    self.__reads = []
    self.__cipos = kwargs.get("cipos", [0, 0])
    self.__ciend = kwargs.get("ciend", [0, 0])
    self.__cprefix = kwargs.get("cprefix", 0)
    self.__csuffix = kwargs.get("csuffix", 0)

    for arg in args: # save reads
      self.__reads.append(arg)

  def getType(self):
    return self.__type

  def getMethod(self):
    return self.__method

  def getReference(self):
    return self.__reference

  def getStart(self):
    return self.__start

  def getEnd(self):
    return self.__end

  def getSequence(self):
    return self.__seq

  def getReferenceSequence(self):
    return self.__refseq

  def getCipos(self):
    return self.__cipos

  def getCiend(self):
    return self.__ciend

  def getCprefix(self):
    return self.__cprefix

  def getCsuffix(self):
    return self.__csuffix

