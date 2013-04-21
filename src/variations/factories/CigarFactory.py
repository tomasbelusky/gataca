#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "21.04. 2013"

import re

from BaseFactory import BaseFactory
from src.interface.interface import *
from src.interface.Settings import Settings
from src.variations.Variation import Variation

class CigarFactory(BaseFactory):
  """
  Creating factory of variations in CIGAR string and MD tag
  """
  op = enum(# operations
             ALIGNMENT=0,
             INSERTION=1,
             DELETION=2,
             SKIPPED=3,
             SOFTCLIP=4,
             HARDCLIP=5,
             PADDING=6,
             MATCH=7,
             MISMATCH=8)
  sums = (op.ALIGNMENT, op.INSERTION, op.SOFTCLIP, op.MATCH, op.MISMATCH) # read seqence length
  abbr = 'MIDNSHP=X' # abbreviations
  faStates = enum(# states of FA for finding SNPs and deletions from MD tag
                  START=0,
                  MATCH=1,
                  DELETION=2)
  bases = re.compile(r'[A-Z]', re.I) # re that folds all letters

  def __init__(self, sample):
    """
    Initialize variables
    """
    BaseFactory.__init__(self, sample)

  def snpIndels(self, read):
    """
    Get SNPs and indels from CIGAR string and MD tag
    """
    if self._sample.readOutOfRegion(Settings.REFERENCE, Settings.START, Settings.END, read.sam):
      return

    tags = dict(read.sam.tags)
    mdtag = tags.get("MD", "")
    pos = read.pos
    insPos = pos
    queryIndex = 0
    insIndex = 0
    refname = self._sample.getRefName(read.tid)
    intervals = [[read.pos, read.end]]

    for operation, length in read.sam.cigar: # look at all operations and their lengths in CIGAR
      if operation in [CigarFactory.op.ALIGNMENT, CigarFactory.op.SOFTCLIP, CigarFactory.op.MATCH, CigarFactory.op.MISMATCH]: # shift index and position
        insIndex += length
        insPos += length

        if operation == CigarFactory.op.SOFTCLIP: # shift query index and skip parsing MD
          queryIndex += length
          continue
      elif operation == CigarFactory.op.INSERTION: # store insertion and skip parsing MD
        newIndex = insIndex + length
        tmpPos = insPos - 1
        refseq = self._sample.fetchReference(read.tid, tmpPos, insPos)
        yield Variation(Variation.vtype.INS, refname, tmpPos, None, refseq, Variation.mtype.CIGAR_MD, info={'svlen': length, 'intervals': intervals})
        insIndex = newIndex
        queryIndex += length
        continue
      elif operation != CigarFactory.op.DELETION: # skip parsing MD
        continue

      state = CigarFactory.faStates.START
      intIndex = 0
      match = ""
      lenDeletion = 0
      saveSNP = False
      delVariation = None

      while mdtag: # parse MD string from read tags
        sign = mdtag[0]

        if state == CigarFactory.faStates.START: # START STATE ---------------------
          mdtag = mdtag[1:]

          if sign.isdigit(): # matches
            match += sign
            state = CigarFactory.faStates.MATCH
          elif CigarFactory.bases.match(sign):  # SNP
            saveSNP = True
          elif sign == '^': # deletion
            pos -= 1
            state = CigarFactory.faStates.DELETION
        elif state == CigarFactory.faStates.MATCH: # MATCH STATE -------------------
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
                state = CigarFactory.faStates.DELETION
              else: # SNP
                saveSNP = True
        elif state == CigarFactory.faStates.DELETION: # DELETION STATE -------------
          if CigarFactory.bases.match(sign): # deleted bases
            mdtag = mdtag[1:]
            lenDeletion += 1
            delVariation = Variation(Variation.vtype.DEL, refname, pos, None, read.sam.seq[queryIndex], Variation.mtype.CIGAR_MD, info={'svlen': lenDeletion, 'intervals': intervals, 'end': pos + lenDeletion})
          else: # end of deletion
            yield delVariation
            intIndex += lenDeletion
            lenDeletion = 0
            state = CigarFactory.faStates.START

        if saveSNP: # save SNP
          yield Variation(Variation.vtype.SNP, refname, pos, read.sam.seq[queryIndex], sign, Variation.mtype.CIGAR_MD, info={'intervals': intervals})
          pos += 1
          queryIndex += 1
          intIndex += 1
          state = CigarFactory.faStates.START
          saveSNP = False

        if intIndex >= length: # go to next cigar operation
          break

      if state == CigarFactory.faStates.DELETION: # add deletion
        yield delVariation
