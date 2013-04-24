#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import copy
from datetime import date
from bx.intervals.intersection import Intersecter, Interval

from Variation import Variation
from clusters.SnpCluster import SnpCluster
from clusters.StructuralCluster import StructuralCluster
from factories.PairFactory import PairFactory
from factories.SplitFactory import SplitFactory
from src.interface.Settings import Settings
from src.interface.interface import *
from src.resources.VcfCreator import VcfCreator
from src.variations.factories.CigarFactory import CigarFactory

class Detector:
  """
  Detector of variations
  """

  def __init__(self, sample, refgenome, output):
    """
    Initialize variables
    """
    self.__sample = sample
    self.__cigarFactory = CigarFactory(self.__sample)
    self.__pairFactory = PairFactory(self.__sample)
    self.__splitFactory = SplitFactory(self.__sample)

    self.__vcfCreator = VcfCreator(refgenome.filename, output)
    self.__vcfCreator.addHeader('fileDate', date.today().strftime('%Y%m%d'))
    self.__vcfCreator.addHeader('source', 'gataca')
    self.__vcfCreator.addInfo('SVTYPE', 1, 'String', 'Type of structural variant')
    self.__vcfCreator.addInfo('END', 1, 'Integer', 'End position of the variant')
    self.__vcfCreator.addInfo('IMPRECISE', 0, 'Flag', 'Imprecise structural variation')
    self.__vcfCreator.addInfo('CILEN', 2, 'Integer', 'Confidence interval of length of variation')
    self.__vcfCreator.addInfo('CPOS', 1, 'Integer', 'Confidence of POS')
    self.__vcfCreator.addInfo('CEND', 1, 'Integer', 'Confidence of END')
    self.__vcfCreator.addInfo('TRACHROM', 1, 'String', 'Chromosome of translocated sequence')
    self.__vcfCreator.addInfo('TRAPOS', 1, 'Integer', 'Start position of translocated sequence')
    self.__vcfCreator.addInfo('TRAEND', 1, 'Integer', 'End position of translocated sequence')
    self.__vcfCreator.addInfo('TRACPOS', 1, 'Integer', 'Confidence of start position of translocated sequence')
    self.__vcfCreator.addInfo('TRACEND', 1, 'Integer', 'Confidence of end position of translocated sequence')
    self.__vcfCreator.addInfo('MAX', 1, 'Integer', 'Maximal end position of deleted sequence')
    self.__vcfCreator.addInfo('CONF', 'A', 'Integer', 'Confidence of each variation')
    self.__vcfCreator.addAlt('DEL', 'Deletion')
    self.__vcfCreator.addAlt('INS', 'Insertion')
    self.__vcfCreator.addAlt('INV', 'Inversion')
    self.__vcfCreator.addAlt('DUP', 'Duplication')
    self.__vcfCreator.addAlt('DUP:TANDEM', 'Tandem duplication')
    self.__vcfCreator.addAlt('INS:TRA', 'Translocation')

    for contig in self.__sample.getRefSequences(): # add contigs
      self.__vcfCreator.addContig(contig)

    self.__vcfCreator.writeHeader()

  def __emptyVariables(self):
    self.__variations = dict((x, []) for x in self.__sample.getReferences())

    self.__clusters = dict((x, {}) for x in self.__sample.getReferences())
    self.__finalClusters = dict((x, []) for x in self.__sample.getReferences())
    self.__clusterTree = dict((x, Intersecter()) for x in self.__sample.getReferences())

    self.__opposites = copy.deepcopy(self.__clusters)
    self.__finalOpposites = []
    self.__oppositeTree = copy.deepcopy(self.__clusterTree)

  def __changeInterval(self, item, oldStart, oldEnd, newStart, newEnd, intervals, tree):
    """
    Change interval where item belongs
    """
    if newStart != oldStart or newEnd != oldEnd:
      intervals[(oldStart, oldEnd)].remove(item)

      if not intervals.get((newStart, newEnd), None):
        intervals[(newStart, newEnd)] = [item]
        tree.add_interval(Interval(newStart, newEnd))
      else:
        intervals[(newStart, newEnd)].append(item)

  def __removeInterval(self, intervals, start, end):
    """
    Remove interval if it's empty
    """
    if not len(intervals.get((start, end), [None])):
      del intervals[(start, end)]

  def __addVariationIntoOpposite(self, var, fromCluster):
    """
    Add variation into oppostine cluster
    """
    ref = var.getReference()
    added = False

    for interval in self.__oppositeTree[ref].find(var.getMaxStart()-1, var.getMaxEnd()+1):
      items = self.__opposites[ref].get((interval.start, interval.end), [])

      for i in reversed(range(len(items))): # go through all cluster's variations
        newVariation = items[i][0].helpDecide(items[i][1], var, fromCluster)

        if newVariation: # joined
          added = True
          self.__changeInterval(items[i], interval.start, interval.end,
                                newVariation.getMaxStart(), newVariation.getMaxEnd(),
                                self.__opposites[ref], self.__oppositeTree[ref])

      self.__removeInterval(self.__opposites[ref], interval.start, interval.end)

    return added

  def __addOppositeCluster(self, cluster):
    """
    Add oposite cluster into list or extend existing one
    """
    append = True

    for var in cluster.getVariations(): # try to add variation into opposite
      append &= not self.__addVariationIntoOpposite(var, True)

    if append: # append cluster
      self.__finalOpposites.append(cluster)

      for index, var in enumerate(cluster.getVariations()):
        ref = var.getReference()
        interval = (var.getMaxStart(), var.getMaxEnd())
        self.__opposites[ref][interval] = self.__opposites[ref].get(interval, []) + [(cluster, index)]
        self.__oppositeTree[ref].add_interval(Interval(interval[0], interval[1]))

  def __processOppositeClusters(self):
    """
    Find winning variations and append them with unused variations into all variations
    """
    for ref in self.__variations: # help clusters to decide about winning variations
      for i in range(len(self.__variations[ref])):
        variation = self.__variations[ref].pop(0)

        if not self.__addVariationIntoOpposite(variation, False):
          self.__variations[ref].append(variation)

    unused = set()
    helpers = set()

    while len(self.__finalOpposites):
      cluster = self.__finalOpposites.pop()
      cluster.process()
      winner = cluster.getWinner()

      if winner: # append winner
        self.__variations[winner.getReference()].append(winner)
        helpers |= cluster.getHelpers()

      unused |= cluster.getUnused()

    for var in unused.difference(helpers):
      self.__variations[var.getReference()].append(var)

  def __makeClusters(self):
    """
    Make clusters from finded variations to join alleles together
    """
    for ref in sorted(self.__variations, key=lambda r: self.__sample.getRefIndex(r)): # all references
      for variation in sorted(self.__variations[ref], key=lambda v: v.getStart()): # all variations
        ref = variation.getReference()
        start = variation.getMaxStart()
        end = variation.getMaxEnd()
        added = False

        for interval in self.__clusterTree[ref].find(start-1, end+1): # compare similarity with overlapped clusters
          clusters = self.__clusters[ref].get((interval.start, interval.end), [])

          for i in reversed(range(len(clusters))):
            if clusters[i].add(variation): # can add into cluster
              added = True
              self.__changeInterval(clusters[i], interval.start, interval.end,
                                    clusters[i].getStart(), clusters[i].getEnd(),
                                    self.__clusters[ref], self.__clusterTree[ref])

          self.__removeInterval(self.__clusters[ref], interval.start, interval.end)

        if not added: # create new cluster
          if variation.getType() == Variation.vtype.SNP:
            cluster = SnpCluster(ref, self.__sample, variation)
          else:
            cluster = StructuralCluster(ref, self.__sample, variation)

          self.__finalClusters[ref].append(cluster)
          self.__clusters[ref][(start, end)] = self.__clusters[ref].get((start, end), []) + [cluster]
          self.__clusterTree[ref].add_interval(Interval(start, end))

  def __write(self):
    """
    Write clusters with founded variations in VCF format
    """
    for ref in self.__sample.getReferences(): # add all clusters
      for cluster in sorted(self.__finalClusters[ref], key=lambda c: c.getActualStart()):
        if not self.__sample.clusterOutOfRegion(Settings.REFERENCE, Settings.START, Settings.END, cluster):
          self.__vcfCreator.writeRecord(cluster.toString())

  def __getSnpIndels(self, read):
    """
    Get SNPs and indels from CIGAR string and MD tag
    """
    for variation in self.__cigarFactory.snpIndels(read):
      self.__variations[variation.getReference()].append(variation)

  def start(self):
    """
    Start finding variations
    """
    self.__sample.preprocessing()
    self.__emptyVariables()

    # fetch paired reads (can be also singleton or unmapped mate)
    for paired in self.__sample.fetchPairs(reference=Settings.REFERENCE, start=Settings.START, end=Settings.END):
      if paired.isSingle():
        self.__getSnpIndels(paired.read)
      else:
        self.__getSnpIndels(paired.read)
        self.__getSnpIndels(paired.mate)
        variation = None

        if paired.isNormal():
          if paired.isRearranged():
            if paired.isInterchromosomal():
              self.__addOppositeCluster(self.__pairFactory.translocationRearranged(paired))
            elif paired.hasOverlap():
              variation = self.__pairFactory.overlapRearranged(paired)
            else:
              self.__addOppositeCluster(self.__pairFactory.rearrangement(paired))
          elif paired.isInterchromosomal():
            if not paired.read.isInverted() and not paired.mate.isInverted():
              self.__addOppositeCluster(self.__pairFactory.translocation(paired))
          elif paired.read.isInverted():
            variation = self.__pairFactory.inversionRead(paired)
          elif paired.mate.isInverted():
            variation = self.__pairFactory.inversionMate(paired)
          elif paired.hasOverlap():
            variation = self.__pairFactory.overlap(paired)
          elif paired.actualSize() < self.__sample.getMinInsertSize():
            self.__addOppositeCluster(self.__pairFactory.smallInsertSize(paired))
          elif paired.actualSize() > self.__sample.getMaxInsertSize():
            self.__addOppositeCluster(self.__pairFactory.bigInsertSize(paired))
        elif (paired.isReadSplit() or paired.isMateSplit()) and paired.splitpair.splitread.isUsable():
          if paired.splitpair.splitread.left.isUnmapped():
            if not paired.splitpair.splitread.right.isInverted():
              if paired.splitpair.isReadFirst():
                variation = self.__splitFactory.leftEncloseInsertion(paired.splitpair.read, paired.splitpair.splitread)
              else:
                variation = self.__splitFactory.leftInsertion(paired.splitpair.splitread)
          elif paired.splitpair.splitread.right.isUnmapped():
            if not paired.splitpair.splitread.left.isInverted():
              if paired.splitpair.isReadFirst():
                variation = self.__splitFactory.rightInsertion(paired.splitpair.splitread)
              else:
                variation = self.__splitFactory.rightEncloseInsertion(paired.splitpair.read, paired.splitpair.splitread)
          elif paired.splitpair.splitread.left.isInverted():
            if not paired.splitpair.partsEnclosePair() and not paired.splitpair.splitread.isRearranged() and not paired.splitpair.splitread.right.isInverted():
              variation = self.__splitFactory.inversionLeft(paired.splitpair.splitread)
          elif paired.splitpair.splitread.right.isInverted():
            if not paired.splitpair.partsEnclosePair() and not paired.splitpair.splitread.isRearranged():
              variation = self.__splitFactory.inversionRight(paired.splitpair.splitread)
          elif paired.splitpair.hasOverlap():
            variation = self.__splitFactory.overlapPair(paired.splitpair, paired.splitpair.splitread)
          elif paired.splitpair.partsEnclosePair():
            if paired.splitpair.splitread.isRearranged():
              self.__addOppositeCluster(self.__splitFactory.encloseRearrangement(paired))
            else:
              self.__addOppositeCluster(self.__splitFactory.enclose(paired))
          elif paired.splitpair.splitread.isRearranged():
            if paired.splitpair.splitread.hasGap():
              self.__addOppositeCluster(self.__splitFactory.rearrangement(paired))
            else:
              variation = self.__splitFactory.overlapRearrangedParts(paired.splitpair.splitread)
          elif paired.splitpair.splitread.hasGap():
            self.__addOppositeCluster(self.__splitFactory.gap(paired))
          elif paired.splitpair.splitread.hasOverlap():
            variation = self.__splitFactory.overlapParts(paired.splitpair.splitread)
          elif paired.splitpair.splitread.hasInsertion():
            variation = self.__splitFactory.normalInsertion(paired.splitpair.splitread)

        if variation: # append variation
          self.__variations[variation.getReference()].append(variation)

    self.__processOppositeClusters()
    self.__makeClusters()
    self.__write()
    self.__vcfCreator.close()
    return 0
