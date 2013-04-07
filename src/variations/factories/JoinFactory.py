#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.04. 2013"

from BaseFactory import BaseFactory
from variations.Variation import Variation

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
    cpos = 0

    if key in info2:
      cpos = info2[key]

    if key in info1:
      if cpos == '-':
        cpos = info1[key]
      elif info1[key] != '-':
        cpos = min(cpos, info1[key])

    return {} if cpos == 0 else {key : cpos}

  def __countCilen(self, info1, info2):
    """
    Count confidence length
    """
    cilen = [0,0]

    if 'cilen' in info1 and 'cilen' in info2:
      cilen = [min(info1['cilen'][0], info2['cilen'][0]), max(info1['cilen'][1], info2['cilen'][1])]
    elif 'cilen' in info1:
      cilen = info1['cilen']
      cilen[0] = min(cilen[0], info2['svlen'])
      cilen[1] = min(cilen[1], info2['svlen'])
    elif 'cilen' in info2:
      cilen = info2['cilen']
      cilen[0] = min(cilen[0], info1['svlen'])
      cilen[1] = min(cilen[1], info1['svlen'])

    return {'cilen' : cilen}

  def __countCend(self, info1, info2, key='cend'):
    """
    Count confidence end position
    """
    cend = 0

    if key in info2:
      cend = info2[key]

    if key in info1:
      if cend == '+':
        cend = info1[key]
      elif info1[key] != '+':
        cend = max(cend, info1[key])

    return {} if cend == 0 else {key : cend}

  def __shiftConfidence(self, info, plus, key):
    """
    Shift confidence position in key
    """
    if key in info and info[key] not in ('+', '-'):
      info[key] += plus

  def __countDepth(self, info1, info2):
    """
    Count depth
    """
    return {'depth' : info1['depth'] + info2['depth']}

  def __countMax(self, info1, info2):
    """
    Count max position
    """
    if 'max' in info1 and 'max' in info2:
      return {'max' : max(info1['max'], info2['max'])}
    elif 'max' in info1:
      return {'max' : info1['max']}
    elif 'max' in info2:
      return {'max' : info2['max']}
    else:
      return {}

  def __boundaryVariaton(self, first, second, vtype):
    """
    Join variation with both start and end position
    """
    info1 = first.getInfo()
    info2 = second.getInfo()

    # check overlap
    if (first.getEnd() < second.getMaxStart() and info2.get('cpos', 0) != '-') or \
       (first.getMaxEnd() < second.getStart() and info1.get('cend', 0) != '+'):
      return None

    pos = max(first.getStart(), second.getMaxStart())
    info = {}
    self.__shiftConfidence(info2, second.getStart() - pos, 'cpos')
    info['end'] = min(first.getMaxEnd(), second.getEnd())
    self.__shiftConfidence(info1, first.getEnd() - info['end'], 'cend')
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCend(info1, info2))
    info.update(self.__countDepth(info1, info2))
    self._repairInfo(pos, info)
    refseq = self._sample.fetchReference(self._sample.getRefIndex(first.getReference()), pos, pos+1)
    return Variation(vtype, first.getReference(), pos, None, refseq, Variation.mtype.JOINED, info=info)

  def __insertedVariation(self, first, second):
    """
    Join inserted variations
    """
    info1 = first.getInfo()
    info2 = second.getInfo()

    if first.getStart() != second.getStart() and 'cpos' not in info2:
      return None, None

    if info2['cpos'] != '-' and first.getStart() < (second.getStart() + info2['cpos']):
      return None, None

    self.__shiftConfidence(info1, first.getStart() - second.getStart(), 'cpos')
    pos = second.getStart()
    info = {}
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCilen(info1, info2))
    info.update(self.__countDepth(info1, info2))
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

    # check overlap
    if info1.get('max', -1) < second.getStart() and (first.getStart() + info1.get('svlen', 0)) < second.getStart():
      return None

    self.__shiftConfidence(info1, first.getStart() - second.getStart(), 'cpos')
    pos = second.getStart()
    info = {}
    info.update(self.__countCpos(info1, info2))
    info.update(self.__countCilen(info1, info2))
    info.update(self.__countMax(info1, info2))
    info.update(self.__countDepth(info1, info2))
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

  def translocation(self, first, second):
    """
    Join translocations
    """
    info1 = first.getInfo()
    info2 = second.getInfo()
    traMaxPos2 = info2['trapos'] + (0 if info2.get('tracpos', 0) == '-' else info2.get('tracpos', 0))
    traMaxEnd1 = info1['traend'] + (0 if info1.get('tracend', 0) == '+' else info1.get('tracend', 0))

    # check overlap
    if first.getReference() != second.getReference() or info1['traend'] < traMaxPos2 or traMaxEnd1 < info2['trapos']:
      return None

    (pos, info) = self.__insertedVariation(first, second)

    if pos is None:
      return None

    # count start position of translocated sequence
    info['trapos'] = max(info1['trapos'], traMaxPos2)
    self.__shiftConfidence(info2, info2['trapos'] - info['trapos'], 'tracpos')
    info.update(self.__countCpos(info1, info2, 'tracpos'))

    # count end position of translocated sequence
    info['traend'] = min(traMaxEnd1, info2['traend'])
    self.__shiftConfidence(info1, info1['traend'] - info['traend'], 'tracend')
    info.update(self.__countCend(info1, info2, 'tracend'))

    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, first.getReference(), pos, None, second.getReferenceSequence(), Variation.mtype.JOINED, info=info)

  def insertionTranslocation(self, first, second):
    """
    Join insertion witn translocation into new translocation
    """
    (pos, info) = self.__insertedVariation(first, second)

    if pos is None:
      return None

    if second.getType() == Variation.vtype.INS: # set insertion as first
      tmp = copy.deepcopy(first)
      first = copy.deepcopy(second)
      second = tmp

    info1 = first.getInfo()
    info2 = second.getInfo()
    minTraLen = info2['traend'] - info2['trapos']
    maxTraLen = (info2['traend'] + info2.get('tracend', 0)) - (info2['trapos'] + info2.get('tracpos', 0))
    info.update(self.__countCilen(info1, {'cilen' : [minTraLen, maxTraLen]}))
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, first.getReference(), pos, None, second.getReferenceSequence(), Variation.mtype.JOINED, info=info)
