#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import pysam
import math
from interface import *
from Read import Read

class Sample:
  """
  Represents sample and also hold file with reference genome
  """
  window = 100 # length of window for getting coverage
  core = 0.1 # core from <0,1> interval for getting insert length
  minCoreCount = 10 # minimal count of iterms in core
  reGCcontent = re.compile('[G|C]', re.I) # re for getting length of GC content
  ptype = enum( # policy type
    FR=0, # forward/reverse
    RF=1) # reverse/forward

  def __init__(self, filename, refgenome, policy=ptype.FR):
    """
    Initialize variables
    """
    self.__reads = pysam.Samfile(filename)
    self.__refgenome = refgenome
    self.__policy = policy
    self.__minInsertLength = 0
    self.__maxInsertLength = 0
    self.__avgInsertLength = 0
    self.__stdInsertLength = 0
    self.__insertLengths = []
    self.__coverage = dict((ref, {}) for ref in range(self.__reads.nreferences))
    self.__gcContent = {}

    # set references on methods which are based on policy
    if self.__policy == Sample.ptype.FR:
      Read._rearrangedRef = Read._isRearrangedFR
      Read._invertedRef = Read._isInvertedFR
    elif self.__policy == Sample.ptype.RF:
      Read._rearrangedRef = Read._isRearrangedRF
      Read._invertedRef = Read._isInvertedRF

  def getWindow(self, position):
    """
    Get window where postiion belongs
    """
    return int(math.floor(position / (Sample.window + 0.0))) * Sample.window

  def getMinInsertLength(self):
    """
    Return minimal insert length
    """
    return self.__minInsertLength

  def getMaxInsertLength(self):
    """
    Return maximal insert length
    """
    return self.__maxInsertLength

  def getAvgInsertLength(self):
    """
    Return average insert length
    """
    return self.__avgInsertLength

    return self.__avgInsertLength

  def getInexactCoverage(self, reference, position):
    """
    Return inexact repaired coverage from GC content in reference and position
    """
    pos = self.getWindow(position)
    return self.__coverage[reference].get(pos, 0)

  def getInexactCoverage(self, reference, start, end):
    """
    Return inexact repaired coverage from GC content in reference and interval
      Faster then exact method
    """
    startPos = self.getWindow(start)
    endPos = self.getWindow(end) + 1
    step = (endPos - startPos) / Sample.window
    coverage = 0
    count = 0

    for i in range(startPos, endPos, Sample.window):
      coverage += self.__coverage[reference].get(i, 0)
      count += 1

    return int(float(coverage) / count if count > 0 else 0)

  def getExactCoverage(self, reference, start, end):
    """
    Return exact coverage in reference and interval
      Slower then inexact method
    """
    count = 0

    for read in self.__reads.fetch(reference=self.__reads.references[reference], start=start, end=end+1):
      if not read.is_unmapped:
        count += 1

    return count

  def getMate(self, read):
    """
    Return mate of read: None if mate is unmapped or read is single
    """
    try:
      return self.__reads.mate(read)
    except ValueError:
      return None

  def getRefSequences(self):
    """
    Return reference sequences from BAM header
    """
    if 'SQ' in self.__reads.header:
      for record in self.__reads.header['SQ']:
        yield record

  def getReferences(self):
    """
    Return list of reference names
    """
    return self.__reads.references

  def getRefName(self, rindex):
    """
    Return name of reference on rindex
    """
    return self.__reads.references[rindex]

  def getRefIndex(self, ref):
    """
    Get index of reference name ref
    """
    return self.__reads.references.index(ref)

  def preprocessing(self, reference=None, start=None, end=None):
    """
    Count and repair coverage from GC content and count insert length's statistics
    """
    for paired in self.fetchPairs(reference, start, end):
      if paired.isSingle(): # count only coverage for read
        self.__addCoverage(paired.read())
      else: # count coverage and insert length
        self.__addCoverage(paired.read())
        self.__addCoverage(paired.mate())
        self.__addInsertLength(paired)

    self.__repairGCcontent()
    allCount = len(self.__insertLengths)
    allHalf = int(allCount / 2)
    coreCount = max(Sample.minCoreCount, int(allCount * Sample.core))

    if allCount < coreCount: # give all values into core
      coreInsertLengths = self.__insertLengths
    else: # create core from middle of values
      coreHalf = int(math.ceil(coreCount / 2.0))
      coreInsertLengths = sorted(self.__insertLengths)[allHalf-coreHalf:allHalf+coreHalf]

    count = float(len(coreInsertLengths))

    if count: # there are paired reads
      self.__sumInsertLength = sum(coreInsertLengths)
      self.__avgInsertLength = self.__sumInsertLength / count
      tmpStdSum = sum([pow(x-self.__avgInsertLength, 2) for x in coreInsertLengths])
      self.__stdInsertLength = math.sqrt(tmpStdSum / count)
      std3 = self.__stdInsertLength * 3
      self.__minInsertLength = self.__avgInsertLength - std3
      self.__maxInsertLength = self.__avgInsertLength + std3

  def __addCoverage(self, read):
    """
    Add coverage into window or create new one
    """
    position = self.getWindow(read.pos)

    if position in self.__coverage[read.tid]: # add
      self.__coverage[read.tid][position] += 1
    else: # new
      self.__coverage[read.tid][position] = 1
      self.__addGCcontent(read.tid, position, position + Sample.window)

  def __addInsertLength(self, paired):
    """
    Add insert length into list
    """
    if paired.length() and not paired.isRearranged() and not paired.isReadInverted() and not paired.isMateInverted():
      self.__insertLengths.append(paired.length())

  def __addGCcontent(self, reference, start, end):
    """
    Add GC content into dict or create new one
    """
    sequence = self.fetchReference(reference, start, end)
    gcCount = len(Sample.reGCcontent.findall(sequence))
    atCount = end - start - gcCount
    gcContent = int(round(gcCount / (gcCount + atCount + 0.0) * 100))
    self.__gcContent[gcContent] = self.__gcContent.get(gcContent, []) + [(reference, start)]

  def __repairGCcontent(self):
    """
    Repair coverage form GC content
    """
    allCoverages = []

    for values in self.__coverage.values(): # get all coverage values
      allCoverages.extend(values.values())

    if len(allCoverages): # coverage values exist
      median = findMedian(sorted(allCoverages))

      for content, windows in self.__gcContent.items(): # repair GC content for all windows
        values = []

        for ref, pos in windows: # get median of all windows with same GC content
          values.append(self.__coverage[ref][pos])

        coeficient = median / float(findMedian(sorted(values)))

        for ref, pos in windows: # repair each window with same GC content
          self.__coverage[ref][pos] = int(round(self.__coverage[ref][pos] * coeficient))

  def fetchReference(self, rindex, start, end):
    """
    Fetch sequence of reference genome
    """
    reference = self.__reads.references[rindex]
    return self.__refgenome.fetch(reference=reference, start=start, end=end)

  def fetchPairs(self, reference=None, start=None, end=None):
    """
    Fetch paired reads
      Could also fetch single read if mate is unmapped or read is single
    """
    for read in self.__reads.fetch(reference=reference, start=start, end=end):
      if not read.is_unmapped: # first read must be mapped
        mate = None

        if (reference is None or read.rnext == reference) and (start is None or read.pnext >= start) and (end is None or read.pnext < end):
          mate = self.getMate(read)

        if not mate or Read.isFirst(read, mate):
          yield Read(read, mate, self.__reads.lengths)

  def close(self):
    """
    Close and free resources
    """
    self.__refgenome.close()
    self.__reads.close()
