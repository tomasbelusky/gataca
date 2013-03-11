#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import pysam
import math
from interface import *

class Sample:
  window = 100
  core = 0.1
  minCoreCount = 10
  reGCcontent = re.compile('[G|C]', re.I)
  ptype = enum(FR=0,
               RF=1)

  def __init__(self, filename, refgenome, policy=ptype.FR):
    self.__reads = pysam.Samfile(filename)
    self.__refgenome = refgenome
    self.__policy = policy
    self.__minInsertSize = 0
    self.__maxInsertSize = 0
    self.__avgInsertSize = 0
    self.__stdInsertSize = 0
    self.__insertSizes = []
    self.__rearrangedRef = None
    self.__invertedRef = None
    self.__coverage = dict((ref, {}) for ref in range(self.__reads.nreferences))

    # set references on methods which are based on policy
    if self.__policy == Sample.ptype.FR:
      self.__rearrangedRef = self.__isRearrangedFR
      self.__invertedRef = self.__isInvertedFR
    elif self.__policy == Sample.ptype.RF:
      self.__rearrangedRef = self.__isRearrangedRF
      self.__invertedRef = self.__isInvertedRF

  def getWindow(self, position):
    return int(math.floor(position / (Sample.window + 0.0))) * Sample.window

  def getMinInsertSize(self):
    return self.__minInsertSize

  def getMaxInsertSize(self):
    return self.__maxInsertSize

  def getAvgInsertSize(self):
    return self.__avgInsertSize

  def getCoverage(self, reference, position):
    pos = self.getWindow(position)

    if pos in self.__coverage[reference]:
      return self.__coverage[reference][pos][0]
    else:
      return 0

  def getMate(self, read):
    try:
      return self.__reads.mate(read)
    except ValueError:
      return None

  def getRefSequences(self):
    if 'SQ' in self.__reads.header:
      for record in self.__reads.header['SQ']:
        yield record

  def getReferences(self):
    return self.__reads.references

  def getRefName(self, rindex):
    return self.__reads.references[rindex]

  def isFirst(self, read, mate):
    return (read.pos <= mate.pos and read.tid == mate.tid) or read.tid < mate.tid

  def __isInvertedFR(self, first, read):
    return (first and read.is_reverse) or (not first and not read.is_reverse)

  def __isInvertedRF(self, first, read):
    return (first and not read.is_reverse) or (not first and read.is_reverse)

  def isInverted(self, read, first):
    return self.__invertedRef(first, read)

  def __isRearrangedFR(self, read, mate, first):
    return (first and read.is_reverse and not mate.is_reverse) or (not first and not read.is_reverse and mate.is_reverse)

  def __isRearrangedRF(self, read, mate, first):
    return (first and not read.is_reverse and mate.is_reverse) or (not first and read.is_reverse and not mate.is_reverse)

  def isRearranged(self, read, mate):
    return self.__rearrangedRef(read, mate, self.isFirst(read, mate))

  def preprocessing(self):
    for read, mate in self.fetchPairs():
      if read and mate: # count coverage and insert size
        self.countCoverage(read)
        self.countCoverage(mate)
        self.countInsertSize(read, mate)
      else: # count only coverage for read
        self.countCoverage(read)

    allCount = len(self.__insertSizes)
    allHalf = int(allCount / 2)
    coreCount = max(Sample.minCoreCount, int(allCount * Sample.core))

    if allCount < coreCount: # give all values into core
      coreInsertSizes = self.__insertSizes
    else: # create core from middle of values
      coreHalf = int(math.ceil(coreCount / 2.0))
      coreInsertSizes = sorted(self.__insertSizes)[allHalf-coreHalf:allHalf+coreHalf]

    count = float(len(coreInsertSizes))

    if count: # there are paired reads
      self.__sumInsertSize = sum(coreInsertSizes)
      self.__avgInsertSize = self.__sumInsertSize / count
      tmpStdSum = sum([pow(x-self.__avgInsertSize, 2) for x in coreInsertSizes])
      self.__stdInsertSize = math.sqrt(tmpStdSum / count)
      std3 = self.__stdInsertSize * 3
      self.__minInsertSize = self.__avgInsertSize - std3
      self.__maxInsertSize = self.__avgInsertSize + std3

  def countCoverage(self, read):
    position = self.getWindow(read.pos)

    if position in self.__coverage[read.tid]:
      self.__coverage[read.tid][position][0] += 1
    else:
      content = self.countGCcontent(read.tid, position, position + Sample.window)
      self.__coverage[read.tid][position] = [1, content]

  def countInsertSize(self, read, mate):
    if not read.tlen or self.isRearranged(read, mate) or self.isInverted(read, True) or self.isInverted(mate, False):
      return # don't count rearranged and inverted reads

    if read.rlen == 0:
      read.rlen = len(read.query)
    if mate.rlen == 0:
      mate.rlen = len(mate.query)

    self.__insertSizes.append(read.tlen)# - (read.rlen + mate.rlen))
    """
    else: # count insert size
      if read.tid == mate.tid: # same chromosome
        self.__insertSizes[read.qname] = mate.pos - read.pos - read.rlen
      else: # another chromosome
        lengthBetween = sum(self.__reads.lengths[read.tid:mate.tid])
        self.__insertSizes[read.qname] = lengthBetween + mate.pos - read.pos - read.rlen
    """

  def countGCcontent(self, reference, start, end):
    sequence = self.fetchReference(reference, start, end)
    gcCount = len(Sample.reGCcontent.findall(sequence))
    atCount = end - start - gcCount
    return gcCount / (gcCount + atCount + 0.0) * 100

  def fetchReference(self, rindex, start, end):
    reference = self.__reads.references[rindex]
    return self.__refgenome.fetch(reference=reference, start=start, end=end)

  def fetchPairs(self, reference=None, start=None, end=None):
    for read in self.__reads.fetch(reference=reference, start=start, end=end):
      if not read.is_unmapped: # first read must be mapped
        mate = None

        if (reference is None or read.rnext == reference) and (start is None or read.pnext >= start) and (end is None or read.pnext < end):
          mate = self.getMate(read)

        if not mate or self.isFirst(read, mate):
          yield read, mate

  def close(self):
    self.__refgenome.close()
    self.__reads.close()

  def __exit__(self, type, value, traceback):
    self.close()
