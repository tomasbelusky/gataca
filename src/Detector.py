#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import sys
import pysam
import bisect
from collections import deque
from datetime import date
from interface import *
from VcfCreator import VcfCreator
from Variation import Variation
from Sample import Sample
from Cluster import Cluster

class Detector:
  cigarOp = enum(ALIGNMENT=0,
                 INSERTION=1,
                 DELETION=2,
                 SKIPPED=3,
                 SOFTCLIP=4,
                 HARDCLIP=5,
                 PADDING=6,
                 MATCH=7,
                 MISMATCH=8)
  faStates = enum(START=0,
                  MATCH=1,
                  DELETION=2)
  bases = re.compile(r'[ATGC]', re.I)

  def __init__(self, sample, refgenome, policy=Sample.ptype.FR, reference=None, start=None, end=None):
    self.__sample = sample
    self.__policy = policy
    self.__reference = reference
    self.__start = start
    self.__end = end
    self.__variations = []
    self.__clusters = []
    self.__vcfCreator = VcfCreator(self.__sample, refgenome)

  def start(self):
    #self.__sample.preprocessing()

    # fetch paired reads (can be also singleton or unmapped mate)
    for read, mate in self.__sample.fetchPairs(reference=self.__reference, start=self.__start, end=self.__end):
      for end in [read]: # get SNP and indels from both pairs
        if end: # read is mapped
          self.getSnpIndels(end)

      if not mate:
        continue

      """
      distance = mate.pos - read.pos

      if read.tid != mate.tid: # interchromosomal translocation?
        pass
      elif distance > self.__sample.getMaxInsertSize(): # deletion
        pass
      elif distance < self.__sample.getMinInsertSize(): # insertion
        pass
      """

    self.makeClusters()
    return 0

  def makeClusters(self):
    for variation in self.__variations:
      self.__clusters.append(Cluster(variation.getType(), variation))

    """
    ref = variation.getReference()
    start = variation.getStart()
    end = variation.getEnd() + 1
    overlaps = []

    for i in range(start, end):
      overlaps.extend(self.__clusters[ref].get(i, []))

    if not overlaps: # create new cluster
      cluster = Cluster(variation.getType(), variation)
      self.__clusters[ref][start] = [cluster]
      self.__clusters[ref][end] = [cluster]
    else: # some overlaps exists
      for cluster in overlaps:
        pass
    """

  def getSnpIndels(self, read):
    tags = dict(read.tags)
    mdtag = tags.get("MD", "")
    pos = read.pos
    insPos = pos
    queryIndex = 0
    insIndex = 0
    refname = self.__sample.getRefName(read.tid)

    for operation, length in read.cigar: # look at all operations and their lengths in CIGAR
      if operation in [Detector.cigarOp.ALIGNMENT, Detector.cigarOp.SOFTCLIP, Detector.cigarOp.MATCH, Detector.cigarOp.MISMATCH]: # shift index and position
        insIndex += length
        insPos += length

        if operation == Detector.cigarOp.SOFTCLIP: # shift query index and skip parsing MD
          queryIndex += length
          continue
      elif operation == Detector.cigarOp.INSERTION: # store insertion and skip parsing MD
        newIndex = insIndex + length
        self.__variations.append(Variation(Variation.vtype.INS, refname, insPos, insPos + length, read.seq[insIndex:newIndex], Variation.mtype.CIGAR_MD, read, refseq="."))
        insIndex = newIndex
        queryIndex += length
        continue
      elif operation != Detector.cigarOp.DELETION: # skip parsing MD
        continue

      state = Detector.faStates.START
      intIndex = 0
      match = ""
      deletion = ""

      while mdtag: # parse MD string from read tags
        sign = mdtag[0]

        if state == Detector.faStates.START: # START STATE ---------------------
          mdtag = mdtag[1:]

          if sign.isdigit(): # matches
            match += sign
            state = Detector.faStates.MATCH
          elif Detector.bases.match(sign):  # SNP -> start
            self.__variations.append(Variation(Variation.vtype.SNP, refname, pos, pos, read.seq[queryIndex], Variation.mtype.CIGAR_MD, read, refseq=sign))
            pos += 1
            queryIndex += 1
            intIndex += 1
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
                self.__variations.append(Variation(Variation.vtype.SNP, refname, pos, pos, read.seq[queryIndex], Variation.mtype.CIGAR_MD, read, refseq=sign))
                pos += 1
                queryIndex += 1
                intIndex += 1
                state = Detector.faStates.START
        elif state == Detector.faStates.DELETION: # DELETION STATE -------------
          if Detector.bases.match(sign): # deleted bases
            deletion += sign
            mdtag = mdtag[1:]
          else: # end of deletion
            offset = len(deletion)
            intIndex += len(deletion)
            deletion = "%s%s" % (read.seq[queryIndex], deletion)
            self.__variations.append(Variation(Variation.vtype.DEL, refname, pos, pos + offset, read.seq[queryIndex], Variation.mtype.CIGAR_MD, read, refseq=deletion))
            deletion = ""
            state = Detector.faStates.START

        if intIndex >= length:
          break

      if state == Detector.faStates.DELETION: # add deletion
        self.__variations.append(Variation(Variation.vtype.DEL, refname, pos, pos + len(deletion), read.seq[queryIndex], Variation.mtype.CIGAR_MD, read, refseq=deletion))

  def write(self, filename):
    self.__vcfCreator.setHeader('fileDate', date.today().strftime('%Y%m%d'))
    self.__vcfCreator.setHeader('source', 'gataca')
    self.__vcfCreator.setInfo('SVTYPE', 1, 'String', 'Type of structural variant')
    self.__vcfCreator.setInfo('SVLEN', 1, 'Integer', 'Difference in length between REF and ALT')
    self.__vcfCreator.setInfo('EVENT', 1, 'String', 'ID of event associated to breakend')
    self.__vcfCreator.setInfo('BNDLEN', 1, 'Integer', 'Length of sequence between with breakends')
    self.__vcfCreator.setInfo('END', 1, 'Integer', 'End position of the variant')
    self.__vcfCreator.setInfo('IMPRECISE', 0, 'Flag', 'Imprecise structural variation')
    self.__vcfCreator.setInfo('CIPOS', 2, 'Integer', 'Confidence interval around POS')
    self.__vcfCreator.setInfo('CIEND', 2, 'Integer', 'Confidence interval around END')
    self.__vcfCreator.setInfo('ISEQ', 1, 'String', 'Imprecise inserted sequence')
    self.__vcfCreator.setInfo('CPREFIX', 1, 'Integer', 'Confidence prefix of inserted sequence')
    self.__vcfCreator.setInfo('CSUFFIX', 1, 'Integer', 'Confidence suffix of inserted sequence')
    self.__vcfCreator.setAlt('DEL', 'Deletion')
    self.__vcfCreator.setAlt('INS', 'Insertion')
    self.__vcfCreator.setAlt('INV', 'Inversion')
    self.__vcfCreator.setAlt('BND', 'Breakend')

    for cluster in self.__clusters:
      self.__vcfCreator.setRecord(cluster)

    self.__vcfCreator.write(filename)
