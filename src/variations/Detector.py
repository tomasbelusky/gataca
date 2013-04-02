#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import copy
import re
from collections import OrderedDict
from datetime import date

from factories.PairFactory import PairFactory
from Variation import Variation
from clusters.BaseCluster import BaseCluster
from clusters.OppositeCluster import OppositeCluster
from interface.interface import *
from resources.Read import Read
from resources.Sample import Sample
from resources.VcfCreator import VcfCreator

class Detector:
  """
  Detector of variations
  """
  faStates = enum(# states of FA for finding SNPs and deletions from MD tag
                  START=0,
                  MATCH=1,
                  DELETION=2)
  bases = re.compile(r'[A-Z]', re.I) # re that folds all letters

  def __init__(self, sample, refgenome, policy=Sample.ptype.FR, reference=None, start=None, end=None):
    """
    Initialize variables
    """
    self.__sample = sample
    self.__policy = policy
    self.__reference = reference
    self.__start = start
    self.__end = end
    self.__variations = dict((x, []) for x in self.__sample.getReferences())
    self.__startClusters = OrderedDict((x, OrderedDict()) for x in self.__sample.getReferences())
    self.__endClusters = copy.deepcopy(self.__startClusters)
    self.__oppositeClusters = []
    self.__vcfCreator = VcfCreator(refgenome.filename)

  def start(self):
    """
    Start finding variations
    """
    self.__sample.preprocessing(self.__reference, self.__start, self.__end)

    # fetch paired reads (can be also singleton or unmapped mate)
    for paired in self.__sample.fetchPairs(reference=self.__reference, start=self.__start, end=self.__end):
      if paired.isSingle(): # get only SNPs and indels from single read
        #self.getSnpIndels(paired.read())
        continue

      #self.getSnpIndels(paired.mate())

      variation = None

      if paired.isNormal():
        if paired.isRearranged():
          if paired.isInterchromosomal():
            self.__addOppositeCluster(PairFactory.translocationRearranged(paired, self.__sample))
          elif paired.hasOverlap():
            pass#variation = PairFactory.overlapRearranged(paired, self.__sample)
          else:
            pass#self.__addOppositeCluster(PairFactory.rearrangement(paired, self.__sample))
        elif paired.isReadInverted():
          if paired.isInterchromosomal():
            self.__addOppositeCluster(PairFactory.translocation(paired, self.__sample))
            #variation = PairFactory.inversionReadOnly(paired.read(), self.__sample)
          else:
            pass#variation = PairFactory.inversionRead(paired, self.__sample)
        elif paired.isMateInverted():
          if paired.isInterchromosomal():
            self.__addOppositeCluster(PairFactory.translocation(paired, self.__sample))
            #variation = PairFactory.inversionReadOnly(paired.mate(), self.__sample)
          else:
            pass#variation = PairFactory.inversionMate(paired, self.__sample)
        elif paired.hasOverlap():
          pass#variation = PairFactory.overlap(paired, self.__sample)
        elif paired.actualSize() < self.__sample.getMinInsertSize():
          pass#variation = PairFactory.insertion(paired, self.__sample)
        elif paired.actualSize() > self.__sample.getMaxInsertSize():
          pass#variation = PairFactory.deletion(paired, self.__sample)
      elif paired.isReadSplit() or paired.isMateSplit():
        pass#print "split"

      if variation: # append variation
        self.__variations[variation.getReference()].append(variation)

    self.__processOppositeClusters()
    self.__makeClusters()
    return 0

  def __addOppositeCluster(self, cluster):
    append = True

    for cluster in self.__oppositeClusters: # find overlap with another cluster
      if cluster.extend(cluster):
        append = False

    if append:
      self.__oppositeClusters.append(cluster)

  def __processOppositeClusters(self):
    """
    Find winning variations and append them with unused variations into all variations
    """
    for ref in self.__variations: # help clusters to decide about winning variations
      for i in range(len(self.__variations[ref])):
        variation = self.__variations[ref].pop(0)
        append = True

        for cluster in self.__oppositeClusters: # find overlap with another cluster
          if cluster.process(variation): # help
            append = False

        if append: # not help
          self.__variations[ref].append(variation)

    for cluster in self.__oppositeClusters: # get winning variations and unused variations
      winner = cluster.getWinner()

      if winner:
        self.__variations[winner.getReference()].extend(winner)

      for unused in cluster.getUnused():
        self.__variations[unused.getReference()].append(unused)

  def __makeClusters(self):
    """
    Make clusters from finded variations to join alleles together
    """
    for ref in self.__variations: # all references
      for variation in sorted(self.__variations[ref], key=lambda v: v.getStart()): # all variations
        if Sample.variationOutOfRegion(self.__reference, self.__start, self.__end, variation):
          continue

        ref = variation.getReference()
        start = variation.getStart()
        end = variation.getEnd()
        overlaps = []
        added = False

        for i in range(start, end + 1): # search for overlaps with existing clusters
          overlaps.extend(self.__startClusters[ref].get(i, []))
          overlaps.extend(self.__endClusters[ref].get(i, []))

        processed = []

        if overlaps: # some overlaps exists
          for cluster in overlaps: # compare similarity with all clusters
            if cluster not in processed and cluster.compare(variation): # can add into cluster
              processed.append(cluster)
              added = True
              oldStart = cluster.getStart()
              oldEnd = cluster.getEnd()
              cluster.add(variation)
              newStart = cluster.getStart()
              newEnd = cluster.getEnd()

              if newStart != oldStart: # change start
                self.__startClusters[ref][oldStart].remove(cluster)
                self.__startClusters[ref][newStart] = self.__startClusters[ref].get(newStart, []) + [cluster]

              if newEnd != oldEnd: # change end
                self.__endClusters[ref][oldEnd].remove(cluster)
                self.__endClusters[ref][newEnd] = self.__endClusters[ref].get(newEnd, []) + [cluster]

        if not overlaps or not added: # create new cluster
          c = BaseCluster(ref, self.__sample)
          c.add(variation)
          self.__startClusters[ref][start] = self.__startClusters[ref].get(start, []) + [c]
          self.__endClusters[ref][end] = self.__endClusters[ref].get(end, []) + [c]

  def getSnpIndels(self, read):
    """
    Get SNPs and indels from CIGAR string and MD tag
    """
    if self.__sample.readOutOfRegion(self.__reference, self.__start, self.__end, read):
      return

    tags = dict(read.tags)
    mdtag = tags.get("MD", "")
    pos = read.pos
    insPos = pos
    queryIndex = 0
    insIndex = 0
    refname = self.__sample.getRefName(read.tid)

    for operation, length in read.cigar: # look at all operations and their lengths in CIGAR
      if operation in [Read.ctype.ALIGNMENT, Read.ctype.SOFTCLIP, Read.ctype.MATCH, Read.ctype.MISMATCH]: # shift index and position
        insIndex += length
        insPos += length

        if operation == Read.ctype.SOFTCLIP: # shift query index and skip parsing MD
          queryIndex += length
          continue
      elif operation == Read.ctype.INSERTION: # store insertion and skip parsing MD
        newIndex = insIndex + length
        tmpPos = insPos - 1
        refseq = self.__sample.fetchReference(read.tid, tmpPos, insPos)
        self.__variations[refname].append(Variation(Variation.vtype.INS, refname, tmpPos, refseq + read.seq[insIndex:newIndex], refseq, Variation.mtype.CIGAR_MD, read, end=tmpPos))
        insIndex = newIndex
        queryIndex += length
        continue
      elif operation != Read.ctype.DELETION: # skip parsing MD
        continue

      state = Detector.faStates.START
      intIndex = 0
      match = ""
      deletion = ""
      saveSNP = False
      delVariation = None

      while mdtag: # parse MD string from read tags
        sign = mdtag[0]

        if state == Detector.faStates.START: # START STATE ---------------------
          mdtag = mdtag[1:]

          if sign.isdigit(): # matches
            match += sign
            state = Detector.faStates.MATCH
          elif Detector.bases.match(sign):  # SNP
            saveSNP = True
          elif sign == '^': # deletion
            pos -= 1
            state = Detector.faStates.DELETION
        elif state == Detector.faStates.MATCH: # MATCH STATE -------------------
          if sign.isdigit(): # match
            match += sign
            mdtag = mdtag[1:]
          else: # stop matching
            offset = int(match)
            pos += offset
            queryIndex += offset
            intIndex += offset
            match = ""

            if intIndex < length: # continue parsing
              mdtag = mdtag[1:]

              if sign == '^': # deletion
                pos -= 1
                state = Detector.faStates.DELETION
              else: # SNP
                saveSNP = True
        elif state == Detector.faStates.DELETION: # DELETION STATE -------------
          if Detector.bases.match(sign): # deleted bases
            deletion += sign
            mdtag = mdtag[1:]
            newDeletion = "%s%s" % (read.seq[queryIndex], deletion)
            delVariation = Variation(Variation.vtype.DEL, refname, pos, read.seq[queryIndex], newDeletion, Variation.mtype.CIGAR_MD, read, end=pos + len(deletion))
          else: # end of deletion
            intIndex += len(deletion)
            #self.__variations[refname].append(delVariation)
            self.__addVariation(delVariation)
            deletion = ""
            state = Detector.faStates.START

        if saveSNP: # save SNP
          self.__variations[refname].append(Variation(Variation.vtype.SNP, refname, pos, read.seq[queryIndex], sign, Variation.mtype.CIGAR_MD, read))
          pos += 1
          queryIndex += 1
          intIndex += 1
          state = Detector.faStates.START
          saveSNP = False

        if intIndex >= length: # go to next cigar operation
          break

      if state == Detector.faStates.DELETION: # add deletion
        self.__variations[refname].append(delVariation)

  def write(self, filename):
    """
    Write clusters with founded variations in VCF format
    """
    self.__vcfCreator.addHeader('fileDate', date.today().strftime('%Y%m%d'))
    self.__vcfCreator.addHeader('source', 'gataca')
    self.__vcfCreator.addInfo('SVTYPE', 1, 'String', 'Type of structural variant')
    self.__vcfCreator.addInfo('SVLEN', 1, 'Integer', 'Difference in length between REF and ALT')
    self.__vcfCreator.addInfo('EVENT', 1, 'String', 'ID of event associated to breakend')
    self.__vcfCreator.addInfo('END', 1, 'Integer', 'End position of the variant')
    self.__vcfCreator.addInfo('DEPTH', 'A', 'Integer', 'Read depths supporting each variation in a row')
    self.__vcfCreator.addInfo('FULLDEPTH', 1, 'Integer', 'Full read depth')
    self.__vcfCreator.addInfo('IMPRECISE', 0, 'Flag', 'Imprecise structural variation')
    self.__vcfCreator.addInfo('CILEN', 2, 'Integer', 'Confidence interval of length')
    self.__vcfCreator.addInfo('CIPOS', 2, 'Integer', 'Confidence interval around POS')
    self.__vcfCreator.addInfo('CIEND', 2, 'Integer', 'Confidence interval around END')
    self.__vcfCreator.addAlt('DEL', 'Deletion')
    self.__vcfCreator.addAlt('INS', 'Insertion')
    self.__vcfCreator.addAlt('INV', 'Inversion')
    self.__vcfCreator.addAlt('DUP', 'Duplication')
    self.__vcfCreator.addAlt('INS:TRA', 'Translocation')

    for contig in self.__sample.getRefSequences(): # add contigs
      self.__vcfCreator.addContig(contig)

    for ref in self.__startClusters: # add all clusters
      for index in sorted(self.__startClusters[ref]):
        for cluster in self.__startClusters[ref][index]:
          self.__vcfCreator.addRecord(cluster)

    self.__vcfCreator.write(filename)
