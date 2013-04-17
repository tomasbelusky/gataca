#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.04. 2013"

import sys

from BaseFactory import BaseFactory
from src.variations.Variation import Variation

class JoinFactory(BaseFactory):
  """
  Factory for joining same overlaping variations
  """

  def __init__(self, sample):
    """
    Initialize variables
    """
    BaseFactory.__init__(self, sample)

  def __countCpos(self, info1, info2, key='cpos'):
    """
    Count confidence start position
    """
    if key in info1:
      cpos = min(info1[key], info2[key]) if key in info2 else info1[key]
    elif key in info2:
      cpos = info2[key]
    else:
      cpos = 0

    return {} if cpos == 0 else {key : cpos}

  def __countCilen(self, info1, info2):
    """
    Count confidence length
    """
    if 'cilen' in info1:
      if 'cilen' in info2:
        cilen = [min(info1['cilen'][0], info2['cilen'][0]), max(info1['cilen'][1], info2['cilen'][1])]
      else:
        cilen = info1['cilen']
        cilen[0] = min(cilen[0], info2.get('svlen', sys.maxint))
        cilen[1] = max(info1['cilen'][1], info2.get('svlen', 0))
    elif 'cilen' in info2:
      cilen = info2['cilen']
      cilen[0] = min(cilen[0], info1.get('svlen', sys.maxint))
      cilen[1] = max(info2['cilen'][1], info1.get('svlen', 0))
    else:
      cilen = [0,0]

    return {'cilen' : cilen}

  def __countCend(self, info1, info2, key='cend'):
    """
    Count confidence end position
    """
    if key in info1:
      cend = max(info1[key], info2[key]) if key in info2 else info1[key]
    elif key in info2:
      cend = info2[key]
    else:
      cend = 0

    return {} if cend == 0 else {key : cend}

  def __shiftConfidence(self, info, plus, key):
    """
    Shift confidence position in key
    """
    if key in info:
      info[key] += plus

  def __haveOverlap(self, first, second, info1, info2):
    """
    Chech if variations overlap
    """
    return second.getMaxStart() <= first.getEnd() and second.getStart() <= first.getMaxEnd()

  def __countIntervals(self, info1, info2):
    """
    Join intervals from both variations together
    """
    return {'intervals' : info1['intervals'] + info2['intervals']}

  def __boundaryVariaton(self, first, second, vtype):
    """
    Join variation with both start and end position
    """
    info1 = first.getInfo()
    info2 = second.getInfo()

    if not self.__haveOverlap(first, second, info1, info2):
      return None

    pos = max(first.getStart(), second.getMaxStart())
    info = {}
    self.__shiftConfidence(info2, second.getStart() - pos, 'cpos')
    info['end'] = min(first.getMaxEnd(), second.getEnd())
    self.__shiftConfidence(info1, first.getEnd() - info['end'], 'cend')
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCend(info1, info2))
    info.update(self.__countCilen(info1, info2))
    info.update(self.__countIntervals(info1, info2))
    self._repairInfo(pos, info)
    refseq = self._sample.fetchReference(self._sample.getRefIndex(first.getReference()), pos, pos+1)
    return Variation(vtype, first.getReference(), pos, None, refseq, Variation.mtype.JOINED, info=info)

  def __insertedVariation(self, first, second):
    """
    Join inserted variations
    """
    info1 = first.getInfo()
    info2 = second.getInfo()

    if not self.__haveOverlap(first, second, info1, info2):
      return None, None

    self.__shiftConfidence(info1, first.getStart() - second.getStart(), 'cpos')
    pos = second.getStart()
    info = {}
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCilen(info1, info2))
    info.update(self.__countIntervals(info1, info2))
    self._repairInfo(pos, info)
    return pos, info

  def insertion(self, first, second):
    """
    Join insertions
    """
    (pos, info) = self.__insertedVariation(first, second)

    if pos is None:
      return None

    return Variation(Variation.vtype.INS, first.getReference(), pos, None, second.getReferenceSequence(), Variation.mtype.JOINED, info=info)

  def deletion(self, first, second):
    """
    Join deletions
    """
    info1 = first.getInfo()
    info2 = second.getInfo()

    if not self.__haveOverlap(first, second, info1, info2):
      return None

    self.__shiftConfidence(info1, first.getStart() - second.getStart(), 'cpos')
    pos = second.getStart()
    info = {}
    info['end'] = min(first.getMaxEnd(), second.getEnd())
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCilen(info1, info2))
    info.update(self.__countCend(info1, info2))
    info.update(self.__countIntervals(info1, info2))
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DEL, first.getReference(), pos, None, first.getReferenceSequence(), Variation.mtype.JOINED, info=info)

  def inversion(self, first, second):
    """
    Join inversions
    """
    return self.__boundaryVariaton(first, second, Variation.vtype.INV)

  def duplication(self, first, second):
    """
    Join duplications
    """
    return self.__boundaryVariaton(first, second, Variation.vtype.DUP)

  def tandemDuplication(self, first, second):
    """
    Join tandem duplications (one of them could be interspersed)
    """
    return self.__boundaryVariaton(first, second, Variation.vtype.DUT)

  def translocation(self, first, second):
    """
    Join translocations
    """
    info1 = first.getInfo()
    info2 = second.getInfo()
    traMaxPos2 = info2['trapos'] + info2.get('tracpos', 0)
    traMaxEnd1 = info1['traend'] + info1.get('tracend', 0)

    # check overlap
    if info1['trachrom'] != info2['trachrom'] or info1['traend'] < traMaxPos2 or traMaxEnd1 < info2['trapos']:
      return None

    (pos, info) = self.__insertedVariation(first, second)

    if pos is None:
      return None

    info['trachrom'] = info1['trachrom']

    # count start position of translocated sequence
    info['trapos'] = max(info1['trapos'], traMaxPos2)
    self.__shiftConfidence(info2, info2['trapos'] - info['trapos'], 'tracpos')
    info.update(self.__countCpos(info1, info2, 'tracpos'))

    # count end position of translocated sequence
    info['traend'] = min(traMaxEnd1, info2['traend'])
    self.__shiftConfidence(info1, info1['traend'] - info['traend'], 'tracend')
    info.update(self.__countCend(info1, info2, 'tracend'))

    info.update(self.__countIntervals(info1, info2))
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, first.getReference(), pos, None, second.getReferenceSequence(), Variation.mtype.JOINED, info=info)

  def translocationInsertion(self, first, second):
    """
    Join insertion witn translocation into new translocation
    """
    (pos, info) = self.__insertedVariation(first, second)

    if pos is None:
      return None

    infoTra = first.getInfo() if first.getType() == Variation.vtype.TRA else second.getInfo()
    info['trachrom'] = infoTra['trachrom']
    info['trapos'] = infoTra['trapos']
    info['traend'] = infoTra['traend']

    if 'tracpos' in infoTra:
      info['tracpos'] = infoTra['tracpos']

    if 'tracend' in infoTra:
      info['tracend'] = infoTra['tracend']

    if 'cilen' in info:
      del info['cilen']

    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, first.getReference(), pos, None, second.getReferenceSequence(), Variation.mtype.JOINED, info=info)
