#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import pysam
import math

class Sample:
  window = 100
  reGCcontent = re.compile('[G|C]', re.I)

  def __init__(self, filename, reference):
    self.__reads = pysam.Samfile(filename)
    self.__reference = pysam.Fastafile(reference)
    self.__minInsertSize = 0
    self.__maxInsertSize = 0
    self.__avgInsertSize = 0
    self.__stdInsertSize = 0
    self.__sumInsertSize = 0
    self.__insertSizes = {}
    self.__coverage = {}
    self.__preprocessing()

  def getWindow(self, position):
    return int(math.floor(position / (Sample.window + 0.0))) * Sample.window

  def getMinInsertSize(self):
    return self.__minInsertSize

  def getMaxInsertSize(self):
    return self.__maxInsertSize

  def getStdInsertSize(self):
    return self.__stdInsertSize

  def getAvgInsertSize(self):
    return self.__avgInsertSize

  def getSumInsertSize(self):
    return self.__sumInsertSize

  def getCoverage(self, reference, position):
    pos = self.getWindow(position)

    if pos in self.__coverage[reference]:
      return self.__coverage[reference][pos][0]
    else:
      return 0

  def __preprocessing(self):
    for index, chrom in enumerate(self.__reads.references):
      self.__coverage[chrom] = {}

    for read in self.__reads.fetch():
      if not read.is_unmapped:
        self.countCoverage(read)
  
        if read.is_paired and not read.mate_is_unmapped:
          self.countInsertSize(read)

    if self.__insertSizes: # count statistics on insert size
      count = len(self.__insertSizes) + 0.0
      self.__avgInsertSize = self.__sumInsertSize / count
      tmpStdSum = sum([pow(x-self.__avgInsertSize, 2) for x in self.__insertSizes.values()])
      self.__stdInsertSize = math.sqrt(tmpStdSum / count)
      std3 = self.__stdInsertSize * 3
      self.__minInsertSize = self.__avgInsertSize - std3
      self.__maxInsertSize = self.__avgInsertSize + std3

  def countCoverage(self, read):
    position = self.getWindow(read.pos)
    rindex = self.__reads.references[read.tid]

    if position in self.__coverage[rindex]:
      self.__coverage[rindex][position][0] += 1
    else:
      content = self.countGCcontent(rindex, position, position + Sample.window)
      self.__coverage[rindex][position] = [1, content]

  def countInsertSize(self, read):
    mate = self.__reads.mate(read)

    if (read.pos > mate.pos and read.tid == mate.tid) or read.tid > mate.tid: # counted
      return

    if read.rlen == 0:
      read.rlen = len(read.query)

    if read.tlen: # insert length is present
      if mate.rlen == 0:
        mate.rlen = len(mate.query)

      self.__insertSizes[read.qname] = abs(read.tlen) - (read.rlen + mate.rlen)
    else: # count insert size
      if read.tid == mate.tid: # same chromosome
        self.__insertSizes[read.qname] = mate.pos - read.pos - read.rlen
      else: # another chromosome
        lengthBetween = sum(self.__reads.lengths[read.tid:mate.tid])
        self.__insertSizes[read.qname] = lengthBetween + mate.pos - read.pos - read.rlen

    self.__sumInsertSize += self.__insertSizes[read.qname]

  def countGCcontent(self, reference, start, end):
    sequence = self.__reference.fetch(reference=reference, start=start, end=end)
    gcCount = len(Sample.reGCcontent.findall(sequence))
    atCount = end - start - gcCount
    return gcCount / (gcCount + atCount + 0.0) * 100

  def close(self):
    self.__reference.close()
    self.__reads.close()

  def __exit__(self, type, value, traceback):
    self.close()
