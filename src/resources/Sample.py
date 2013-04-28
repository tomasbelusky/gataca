#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import sys
import math
import pysam
import re
import os
import glob

from src.interface.Settings import Settings
from reads.Paired import Paired
from reads.Read import Read
from reads.SplitPart import SplitPart
from src.tools.Bwa import Bwa
from src.interface.interface import *

class Sample:
  """
  Represents sample and also hold file with reference genome
  """
  TMP_PATH = "tmp"
  FILENAME_TEMPLATE = "%s/%d_" % (TMP_PATH, os.getpid())
  reGCcontent = re.compile('[G|C]', re.I) # re for getting length of GC content

  def __init__(self, filename, refgenome):
    """
    Initialize variables
    """
    self.__reads = pysam.Samfile(filename)
    self.__refgenome = refgenome
    self.__bwa = Bwa()
    self.__splitParts = {}

    self.__minInsertSize = Settings.MIN_INSERT
    self.__maxInsertSize = Settings.MAX_INSERT
    self.__countInsertSize = not (self.__minInsertSize and self.__maxInsertSize)

    self.__minCoverage = Settings.MIN_COVERAGE
    self.__maxCoverage = Settings.MAX_COVERAGE
    self.__coverage = dict((ref, {}) for ref in range(self.__reads.nreferences))
    self.__countCoverage = not (self.__minCoverage and self.__maxCoverage)
    self.__gcContent = {}

    # set references on methods which are based on policy
    if Settings.POLICY == Read.ptype.FR:
      Paired._readStrand = True
      Paired._mateStrand = False
    elif Settings.POLICY == Read.ptype.RF:
      Paired._readStrand = False
      Paired._mateStrand = True

  def __countInterval(self, values, core, minCount):
    """
    Count interval of allowed values
    """
    allCount = len(values)
    allHalf = int(allCount / 2)
    coreCount = max(minCount, int(allCount * core))

    if allCount < coreCount: # give all values into core
      coreValues = values
    else: # create core from middle of values
      coreHalf = int(math.ceil(coreCount / 2.0))
      coreValues = sorted(values)[allHalf-coreHalf:allHalf+coreHalf]

    count = float(len(coreValues))

    if not count:
      return 0, 0

    average = sum(coreValues) / count
    tmpStdSum = sum([pow(x-average, 2) for x in coreValues])
    std3 = math.sqrt(tmpStdSum / count) * 3
    return int(math.ceil(average - std3)), int(math.floor(average + std3))

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

    if not len(allCoverages): # coverage values exist
      return

    coverages = []
    median = findMedian(sorted(allCoverages))

    for windows in self.__gcContent.values(): # repair GC content for all windows
      values = []

      for ref, pos in windows: # get median of all windows with same GC content
        values.append(self.__coverage[ref][pos])

      coeficient = median / float(findMedian(sorted(values)))

      for ref, pos in windows: # repair each window with same GC content
        c = int(round(self.__coverage[ref][pos] * coeficient))
        self.__coverage[ref][pos] = c
        coverages.append(c)

    self.__gcContent.clear()

    if self.__countCoverage:
      (self.__minCoverage, self.__maxCoverage) = self.__countInterval(coverages, Settings.COVERAGE_CORE, Settings.MIN_COVERAGE_COUNT)

  def __addCoverage(self, read):
    """
    Add coverage into window or create new one
    """
    position = self.getWindow(read.pos)

    if position in self.__coverage[read.tid]: # add
      self.__coverage[read.tid][position] += 1
    else: # new
      self.__coverage[read.tid][position] = 1
      self.__addGCcontent(read.tid, position, position + Settings.WINDOW_SIZE)

  def __estimateCoverage(self):
    """
    Count depth of coverages
    """
    for paired in self.fetchPairs(Settings.REFERENCE):
      if paired.isSingle():
        self.__addCoverage(paired.read)
      else:
        self.__addCoverage(paired.read.sam)
        self.__addCoverage(paired.mate.sam)

  def __estimateInterval(self):
    """
    Estimate interval of insert size of properly aligned reads
    """
    insertSizes = []
    count = 0

    for paired in self.fetchPairs(Settings.REFERENCE):
      size = paired.size()

      if size and paired.isNormal(): # must be normal paired with size > 0
        insertSizes.append(size)
        count += 1

        if count == Settings.INSERT_READS: # check limit
          break

    (self.__minInsertSize, self.__maxInsertSize) = self.__countInterval(insertSizes, Settings.INSERT_CORE, Settings.MIN_INSERT_COUNT)

  def __remapping(self):
    """
    Remapping of soft-clipped sequences
    """
    usedFiles = set()
    readDescriptors = {}

    for ref in self.__reads.references: # create tmo FASTQ files
      readname = "%s%s" % (Sample.FILENAME_TEMPLATE, ref)
      readFile = open("%s.fastq" % readname, 'w')
      readDescriptors[ref] = readFile

    for paired in self.fetchPairs(Settings.REFERENCE, Settings.START, Settings.END): # fetch split reads for remapping
      if paired.isReadSplit() or paired.isMateSplit():
        readRef = paired.read if paired.isReadSplit() else paired.mate
        readDescriptors[readRef.reference].write(readRef.getFastqSplits())
        usedFiles.add(readRef.reference)
        self.__splitParts[readRef.sam.qname] = self.__splitParts.get(readRef.sam.qname, []) + readRef.getMappedParts()

    for readFile in readDescriptors.values(): # close FASTQ files
      readFile.close()

    for ref in usedFiles: # do remapping
      actualRef = "%s%s" % (Sample.FILENAME_TEMPLATE, ref)

      with open("%s.fasta" % actualRef, 'w') as refFile: # create file with actual reference
        rindex = self.getRefIndex(ref)
        refFile.write(">%s\n" % ref)
        refFile.write(self.fetchReference(rindex, 0, self.__reads.lengths[rindex]))

      self.__bwa.index(actualRef)
      self.__bwa.align(actualRef)

      with pysam.Samfile("%s.sam" % actualRef) as splitsam:
        for read in splitsam.fetch(): # save remapped parts
          name = read.qname.split(Read.COUNT_SEPARATOR)
          self.__splitParts[name[0]].insert(int(name[1]), SplitPart(read.pos,
                                                                    read.seq,
                                                                    read.cigar,
                                                                    read.mapq,
                                                                    read.is_reverse,
                                                                    read.is_unmapped,
                                                                    True))

      for f in glob.glob("%s*" % actualRef): # remove temporary files
        os.remove(f)

    for f in glob.glob("%s*" % Sample.FILENAME_TEMPLATE): # remove unused temporary FASTQ files
      os.remove(f)

  def preprocessing(self):
    """
    Count and repair coverage from GC content and count insert length's statistics
    """
    self.__remapping()

    if self.__countInsertSize:
      self.__estimateInterval()

    self.__estimateCoverage()
    self.__repairGCcontent()

  def fetchReference(self, rindex, start, end):
    """
    Fetch sequence of reference genome
    """
    return self.__refgenome.fetch(reference=self.getRefName(rindex), start=start, end=end)

  def fetchPairs(self, reference=None, start=None, end=None):
    """
    Fetch paired reads
      Could also fetch single read if mate is unmapped or read is singleton
    """
    for read in self.__reads.fetch(reference=reference, start=start, end=end):
      if not read.is_unmapped: # first read must be mapped
        mate = self.getMate(read)

        # prevent fetching two times same pair
        if not mate or Paired.isFirst(read, mate) or self.readOutOfRegion(reference, start, end, mate):
          paired = Paired(read, mate, self.__reads.lengths, self.__reads.references, self.__splitParts.get(read.qname, []))

          if not paired.isFiltered():
            yield paired

  def fetchTuplePairs(self, reference=None, start=None, end=None):
    """
    Fetch paired reads as Read instances rather than as Paired instance
      Fetch also low quality reads, first read unmapped
    """
    for read in self.__reads.fetch(reference=reference, start=start, end=end):
      if not read.mate_is_unmapped:
        mate = self.getMate(read)

        # prevent fetching two times same pair
        if not mate or read.is_unmapped or Paired.isFirst(read, mate) or self.readOutOfRegion(reference, start, end, mate):
          yield Read(read, True, Paired._readStrand, self.__reads.references), Read(mate, False, Paired._mateStrand, self.__reads.references)

  def readOutOfRegion(self, reference, start, end, read):
    """
    Check if read is out of specified region
    """
    if reference is None:
      return False

    if start is None:
      start = 0

    if end is None:
      end = sys.maxint

    return reference != self.getRefName(read.tid) or \
           (reference == self.getRefName(read.tid) and \
            (end < read.pos or Read.calculateEnd(read) < start))

  def clusterOutOfRegion(self, reference, start, end, cluster):
    """
    Check if cluster is out of specified region
    """
    return (reference is not None and reference != cluster.getReference()) or \
           (end is not None and end < cluster.getActualStart()) or \
           (start is not None and cluster.getEnd() < start)

  def getWindow(self, position):
    """
    Get window where postiion belongs
    """
    return int(math.floor(position / (Settings.WINDOW_SIZE + 0.0))) * Settings.WINDOW_SIZE

  def getMinInsertSize(self):
    """
    Return minimal insert length
    """
    return self.__minInsertSize

  def getMaxInsertSize(self):
    """
    Return maximal insert length
    """
    return self.__maxInsertSize

  def getMinCoverage(self):
    """
    Return minimal coverage value
    """
    return self.__minCoverage

  def getMaxCoverage(self):
    """
    Return maximal coverage value
    """
    return self.__maxCoverage

  def getInexactCoverage(self, reference, start, end):
    """
    Return inexact repaired coverage from GC content
    """
    tid = self.getRefIndex(reference)
    startPos = self.getWindow(start)
    endPos = self.getWindow(end) + 1
    coverage = 0
    count = 0

    for i in range(startPos, endPos, Settings.WINDOW_SIZE):
      coverage += self.__coverage[tid].get(i, 0)
      count += 1

    return int(float(coverage) / count if count > 0 else 0)

  def getExactCoverage(self, reference, start, end):
    """
    Return exact coverage
    """
    count = 0

    for read in self.__reads.fetch(reference=self.__reads.references[reference], start=start, end=end+1):
      r = Read(read, True, True, self.__reads.references)

      if not r.isUnmapped() and not r.isDuplicate() and r.hasMinQuality():
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

  def getHeader(self):
    """
    Return header of BAM file
    """
    return self.__reads.header

  def getLengths(self):
    """
    Return list of reference lengths
    """
    return self.__reads.lengths

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

  def close(self):
    """
    Close and free resources
    """
    self.__refgenome.close()
    self.__reads.close()
