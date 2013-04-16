#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import math
import pysam
import re
import os
import glob

from reads.Paired import Paired
from reads.Read import Read
from reads.SplitPart import SplitPart
from tools.Bwa import Bwa
from interface.interface import *

class Sample:
  """
  Represents sample and also hold file with reference genome
  """
  filenameTemplate = "../tmp/%d_" % os.getpid()
  window = 100 # length of window for getting coverage
  core = 0.1 # core from <0,1> interval for getting insert size
  minCoreCount = 10 # minimal count of iterms in core
  reGCcontent = re.compile('[G|C]', re.I) # re for getting length of GC content
  ptype = enum(# policy type
               FR=0, # forward/reverse
               RF=1) # reverse/forward

  def __init__(self, filename, refgenome, policy=ptype.FR, minInsertSize=0, maxInsertSize=0, countCoverage=True):
    """
    Initialize variables
    """
    self.__reads = pysam.Samfile(filename)
    self.__refgenome = refgenome
    self.__policy = policy
    self.__minInsertSize = minInsertSize
    self.__maxInsertSize = maxInsertSize
    self.__insertSizes = []
    self.__coverage = dict((ref, {}) for ref in range(self.__reads.nreferences))
    self.__gcContent = {}
    self.__bwa = Bwa()
    self.__refFiles = []
    self.__usedFiles = set()
    self.__readFiles = []
    self.__readDescriptors = {}
    self.__splitParts = {}

    if countCoverage:
      self.__coverageRef = self.__countCoverage
    else:
      self.__coverageRef = self.__passCounting

    if self.__minInsertSize and self.__maxInsertSize:
      self.__insertSizeRef = self.__passCounting
    else:
      self.__insertSizeRef = self.__countInsertSize

    # set references on methods which are based on policy
    if self.__policy == Sample.ptype.FR:
      Paired._readStrand = True
      Paired._mateStrand = False
    elif self.__policy == Sample.ptype.RF:
      Paired._readStrand = False
      Paired._mateStrand = True

  def __countCoverage(self, read):
    """
    Add depth of coverage
    """
    self.__addCoverage(read)

  def __countInsertSize(self, paired):
    """
    Add insert size
    """
    self.__addInsertSize(paired)

  def __passCounting(self, *args):
    """
    Don't count coverage/insert size
    """
    pass

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

  def __addInsertSize(self, paired):
    """
    Add insert length into list
    """
    size = paired.size()

    if size and paired.isNormal():
      self.__insertSizes.append(size)

  def __estimateInterval(self):
    """
    Estimate interval of insert size of properly aligned reads
    """
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
      sumIs = sum(coreInsertSizes)
      avgIs = sumIs / count
      tmpStdSum = sum([pow(x-avgIs, 2) for x in coreInsertSizes])
      std3 = math.sqrt(tmpStdSum / count) * 3
      self.__minInsertSize = int(math.ceil(avgIs - std3))
      self.__maxInsertSize = int(math.floor(avgIs + std3))

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

      for windows in self.__gcContent.values(): # repair GC content for all windows
        values = []

        for ref, pos in windows: # get median of all windows with same GC content
          values.append(self.__coverage[ref][pos])

        coeficient = median / float(findMedian(sorted(values)))

        for ref, pos in windows: # repair each window with same GC content
          self.__coverage[ref][pos] = int(round(self.__coverage[ref][pos] * coeficient))

  def __openSplitFastqFiles(self):
    """
    Creating temporary files
    """
    for ref in self.__reads.references:
      readname = "%s%s" % (Sample.filenameTemplate, ref)
      readFile = open("%s.fastq" % readname, 'w')
      self.__readDescriptors[ref] = readFile
      self.__readFiles.append(readname)

  def __remapping(self):
    """
    Remapping of soft-clipped sequences
    """
    for readFile in self.__readDescriptors.values(): # close read files
      readFile.close()

    for ref in self.__usedFiles: # do remapping
      actualRef = "%s%s" % (Sample.filenameTemplate, ref)

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

    for f in glob.glob("%s*" % Sample.filenameTemplate): # remove unused temporary FASTA files
      os.remove(f)

  def preprocessing(self, reference=None, start=None, end=None):
    """
    Count and repair coverage from GC content and count insert length's statistics
    """
    self.__openSplitFastqFiles()

    for paired in self.fetchPairs(reference, start, end):
      if paired.isSingle(): # count only coverage for read
        self.__coverageRef(paired.read)
      else: # count coverage and insert length, save split reads into appropriate files
        if paired.isReadSplit() or paired.isMateSplit():
          readRef = paired.read if paired.isReadSplit() else paired.mate
          self.__readDescriptors[readRef.reference].write(readRef.getFastqSplits())
          self.__usedFiles.add(readRef.reference)
          self.__splitParts[paired.qname] = self.__splitParts.get(paired.qname, []) + readRef.getMappedParts()
        else:
          self.__insertSizeRef(paired)

        self.__coverageRef(paired.read.sam)
        self.__coverageRef(paired.mate.sam)

    self.__remapping()

    # repair content
    if self.__coverageRef != self.__passCounting:
      self.__repairGCcontent()

    # estimate insert length's min and max
    if self.__insertSizeRef != self.__passCounting:
      self.__estimateInterval()

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
          r = Paired(read, mate, self.__reads.lengths, self.__reads.references, self.__splitParts.get(read.qname, []))

          if not r.isFiltered():
            yield r

  def readOutOfRegion(self, reference, start, end, read):
    """
    Check if read is out of specified region
    """
    return (reference is not None and reference != self.getRefName(read.tid)) or \
           (end is not None and end < read.pos) or \
           (start is not None and Read.calculateEnd(read) < start)

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
    return int(math.floor(position / (Sample.window + 0.0))) * Sample.window

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

  def getAvgInsertSize(self):
    """
    Return average insert length
    """
    return self.__avgInsertSize

    return self.__avgInsertSize

  def getInexactCoverage(self, reference, position):
    """
    Return inexact repaired coverage from GC content in reference and position
    """
    pos = self.getWindow(position)
    return self.__coverage[reference].get(pos, 0)

  def getInexactCoverages(self, reference, start, end):
    """
    Return inexact repaired coverage from GC content in reference and interval
      Faster then exact method
    """
    startPos = self.getWindow(start)
    endPos = self.getWindow(end) + 1
    coverage = 0
    count = 0

    for i in range(startPos, endPos, Sample.window):
      coverage += self.__coverage[reference].get(i, 0)
      count += 1

    return int(float(coverage) / count if count > 0 else 0)

  def getExactCoverage(self, reference, start):
    """
    Return exact coverage in reference and position
      Slower then inexact method
    """
    count = 0

    for read in self.__reads.fetch(reference=self.__reads.references[reference], start=start, end=start+1):
      r = Read(read, True, True, self.__reads.references)

      if not r.isUnmapped() and not r.isDuplicate() and r.hasMinQuality():
        count += 1

    return count

  def getExactCoverages(self, reference, start, end):
    """
    Return exact coverage in reference and interval
      Slower then inexact method
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
