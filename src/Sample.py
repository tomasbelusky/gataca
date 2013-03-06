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
    self.__readsCount = 0
    self.__coverage = {}
    self.__preprocessing()

  def getWindow(self, position):
    return int(math.floor(position / (Sample.window + 0.0))) * Sample.window

  def getCoverage(self, reference, position):
    return self.__coverage[reference][self.getWindow(position)][0]

  def __preprocessing(self):
    countedReads = set()

    for index, chrom in enumerate(self.__reads.references):
      self.__coverage[chrom] = {}

    for read in self.__reads.fetch():
      if not read.is_unmapped:
        self.countCoverage(read)

        if read.is_paired and not read.mate_is_unmapped and read not in countedReads:
          self.countInsertSize(read)
          countedReads.add(self.__reads.mate(read))

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

    if read.tlen: # insert length is present
      if read.rlen == 0:
        read.rlen = len(read.query)
      if mate.rlen == 0:
        mate.rlen = len(mate.query)

      self.__readsCount += 1
      self.__sumInsertSize += abs(read.tlen) - (read.rlen + mate.rlen)
    else: # count insert size
      pass

  def countGCcontent(self, reference, start, end):
    sequence = self.__reference.fetch(reference=reference, start=start, end=end)
    gcCount = len(Sample.reGCcontent.findall(sequence))
    atCount = end - start - gcCount
    return gcCount / (gcCount + atCount + 0.0) * 100

  def close():
    self.__reference.close()
    self.__reads.close()

  def __exit__(self, type, value, traceback):
    self.close()
