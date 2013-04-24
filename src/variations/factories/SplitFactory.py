#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "13.04. 2013"

from src.interface.interface import *
from src.variations.Variation import Variation
from src.variations.clusters.OppositeCluster import OppositeCluster
from CigarFactory import CigarFactory
from BaseFactory import BaseFactory

class SplitFactory(BaseFactory):
  """
  Creating factory of variations detected by split read method
  """

  def __init__(self, sample):
    """
    Initialize variables
    """
    BaseFactory.__init__(self, sample)

  def __variation(self, vtype, splitread, pos, info):
    """
    Create variation
    """
    refseq = self._sample.fetchReference(splitread.original.tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(vtype, splitread.original.reference, pos, None, refseq, Variation.mtype.SPLIT_READ, info=info)

  def leftEncloseInsertion(self, read, splitread):
    """
    Create insertion of left part enclosed
    """
    info = {}
    pos = splitread.right.pos - 1
    info['cilen'] = [splitread.left.len, splitread.left.len]
    size = splitread.right.pos - read.end
    info['cilen'][0] += max(0, self._sample.getMinInsertSize() - size)
    info['cilen'][1] += max(0, self._sample.getMaxInsertSize() - size)
    info['intervals'] = [[pos, splitread.right.pos], [pos, splitread.right.pos]]
    return self.__variation(Variation.vtype.INS, splitread, pos, info)
    
  def rightEncloseInsertion(self, read, splitread):
    """
    Create insertion of right part enclosed
    """
    info = {}
    pos = splitread.left.end - 1
    info['cilen'] = [splitread.right.len, splitread.right.len]
    size = read.pos - splitread.left.end
    info['cilen'][0] += max(0, self._sample.getMinInsertSize() - size)
    info['cilen'][1] += max(0, self._sample.getMaxInsertSize() - size)
    info['intervals'] = [[splitread.left.end, splitread.left.end + 1], [splitread.left.end, splitread.left.end + 1]]
    return self.__variation(Variation.vtype.INS, splitread, pos, info)

  def leftInsertion(self, splitread):
    """
    Create insertion of left part
    """
    info = {}
    pos = splitread.right.pos - 1
    info['cilen'] = [splitread.left.len, splitread.right.pos]
    info['intervals'] = [[pos, splitread.right.pos], [pos, splitread.right.pos]]
    return self.__variation(Variation.vtype.INS, splitread, pos, info)

  def rightInsertion(self, splitread):
    """
    Create insertion of right part
    """
    info = {}
    pos = splitread.left.end - 1
    info['cilen'] = [splitread.right.len, self._countMaxLength(splitread.original.tid, splitread.left.end)]
    info['intervals'] = [[splitread.left.end, splitread.left.end + 1], [splitread.left.end, splitread.left.end + 1]]
    return self.__variation(Variation.vtype.INS, splitread, pos, info)

  def normalInsertion(self, splitread):
    """
    Create insertion between split parts
    """
    info = {}
    pos = splitread.left.pos - lengthLeft - 1
    info['svlen'] = 0

    for operation, length in splitread.left.cigar.reverse:
      if operation in (CigarFactory.align.SOFTCLIP, CigarFactory.align.HARDCLIP):
        pos -= length
        info['svlen'] += length
      else:
        break

    for operation, length in splitread.right.cigar:
      if operation in (CigarFactory.align.SOFTCLIP, CigarFactory.align.HARDCLIP):
        info['svlen'] += length
      else:
        break

    info['intervals'] = [[pos, pos + info['svlen']]]
    return self.__variation(Variation.vtype.INS, splitread, pos, info)

  def deletion(self, splitread):
    """
    Create deletion
    """
    info = {}
    pos = splitread.left.end - 1
    info['end'] = splitread.right.pos
    info['svlen'] = info['end'] - pos + 1
    info['intervals'] = [[splitread.left.end, splitread.right.pos]]
    return self.__variation(Variation.vtype.DEL, splitread, pos, info)

  def inversionLeft(self, splitread):
    """
    Create inversion on left part
    """
    info = {}
    pos = splitread.remapped.pos
    info['end'] = splitread.primary.pos
    info['svlen'] = info['end'] - pos + 1
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.INV, splitread, pos, info)

  def inversionRight(self, splitread):
    """
    Create inversion on right part
    """
    info = {}
    pos = splitread.primary.end
    info['end'] = splitread.remapped.end
    info['svlen'] = info['end'] - pos + 1
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.INV, splitread, pos, info)

  def overlapPair(self, splitpair, splitread):
    """
    Create duplication from overlap of remapped part and mate pair
    """
    info = {}
    pos = max(splitpair.read.pos, splitread.remapped.pos) - 1
    info['end'] = min(splitpair.read.end, splitread.remapped.end)

    if splitpair.isReadFirst():
      if splitread.isRearranged():
        info['cend'] = splitread.right.end - info['end']
      else:
        info['cpos'] = info['end'] - splitread.right.pos
    else:
      if splitread.isRearranged():
        info['cpos'] = splitread.right.pos - pos
      else:
        info['cend'] = pos - splitread.left.end

    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info.get('cend', -info.get('cpos', 0))]
    info['intervals'] = [[pos, info['end']], [pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def overlapParts(self, splitread):
    """
    Create tandem duplication from overlapped parts
    """
    info = {}
    pos = splitread.right.pos - 1
    info['end'] = splitread.left.end
    info['svlen'] = info['end'] - pos + 1
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUT, splitread, pos, info)

  def overlapRearrangedParts(self, splitread):
    """
    Create duplication from overlapped and rearranged parts
    """
    info = {}
    pos = splitread.right.pos - 1
    info['end'] = splitread.left.end
    info['cend'] = (splitread.right.pos - splitread.left.pos) + (splitread.right.end - splitread.left.end)
    info['cpos'] = -info['cend']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def leftDuplicationEncloseRearrangement(self, splitread):
    """
    Create duplication of rearranged parts which enclose mate on left side
    """
    info = {}
    pos = splitread.left.pos - 1
    info['end'] = splitread.left.end
    info['cend'] = splitread.right.end - info['end']
    info['cpos'] = -info['cend']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def rightDuplicationEncloseRearrangement(self, splitread):
    """
    Create duplication of rearranged parts which enclose mate on right side
    """
    info = {}
    pos = splitread.right.pos - 1
    info['end'] = splitread.right.end
    info['cpos'] = splitread.left.pos - pos
    info['cend'] = -info['cpos']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def leftDuplicationEnclose(self, splitread):
    """
    Create duplication of parts which enclose mate on left side
    """
    info = {}
    pos = splitread.left.pos - 1
    info['end'] = splitread.left.end
    info['cpos'] = splitread.left.end - splitread.right.pos + splitread.left.len
    info['cend'] = -info['cpos']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def rightDuplicationEnclose(self, splitread):
    """
    Create duplication of parts which enclose mate on right side
    """
    info = {}
    pos = splitread.right.pos - 1
    info['end'] = splitread.right.end
    info['cend'] = splitread.right.pos - splitread.left.end - splitread.right.len
    info['cpos'] = -info['cend']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def leftDuplicationRearrangement(self, splitread):
    """
    Create rearranged duplication on left side
    """
    info = {}
    pos = splitread.left.pos - 1
    info['end'] = splitread.left.end
    info['cend'] = splitread.right.end - info['end']
    info['cpos'] = -info['cend']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def rightDuplicationRearrangement(self, splitread):
    """
    Create rearranged duplication on right side
    """
    info = {}
    pos = splitread.right.pos - 1
    info['end'] = splitread.right.end
    info['cpos'] = splitread.left.pos - pos
    info['cend'] = -info['cpos']
    smallestSize = info['end'] - pos
    info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, splitread, pos, info)

  def leftDuplication(self, splitread):
    """
    Create duplication on left side
    """
    info = {}
    vtype = Variation.vtype.DUP
    pos = splitread.left.pos - 1
    info['end'] = splitread.left.end
    info['cend'] = splitread.right.pos - info['end'] - splitread.left.len

    if info['cend'] < 0:
      info['end'] += info['cend']

      if not (splitread.left.len % 2) and (splitread.left.len / 2) == -info['cend']:
        vtype = Variation.vtype.DUT

      info['svlen'] = info['end'] - pos + 1
    else:
      info['cpos'] = -info['cend']
      smallestSize = info['end'] - pos
      info['cilen'] = [smallestSize, smallestSize + info['cend']]

    info['intervals'] = [[pos, info['end']]]
    return self.__variation(vtype, splitread, pos, info)

  def rightDuplication(self, splitread):
    """
    Create duplication on right side
    """
    info = {}
    vtype = Variation.vtype.DUP
    pos = splitread.right.pos - 1
    info['end'] = splitread.right.end
    info['cend'] = pos - splitread.left.end - splitread.right.len

    if info['cend'] < 0:
      pos -= info['cend'] 

      if not (splitread.right.len % 2) and (splitread.right.len / 2) == -info['cend']:
        vtype = Variation.vtype.DUT

      info['svlen'] = info['end'] - pos + 1
    else:
      smallestSize = info['end'] - pos
      info['cilen'] = [smallestSize, smallestSize + info['cend']]
      info['cpos'] = -info['cend']

    info['intervals'] = [[pos, info['end']]]
    return self.__variation(vtype, splitread, pos, info)

  def leftTranslocationEncloseRearrangement(self, read, splitread):
    """
    Create translocation of rearranged parts which enclose mate on left side
    """
    info = {}
    pos = splitread.right.end - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.left.pos - 1
    info['traend'] = splitread.left.end
    info['tracend'] = read.pos - info['traend']
    info['intervals'] = [[splitread.right.end, splitread.right.end + 1]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def rightTranslocationEncloseRearrangement(self, read, splitread):
    """
    Create translocation of rearranged parts which enclose mate on right side
    """
    info = {}
    pos = splitread.left.pos - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.right.pos - 1
    info['traend'] = splitread.right.end
    info['tracpos'] = read.end - info['trapos']
    info['intervals'] = [[splitread.left.pos - 1, splitread.left.pos]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def leftTranslocationEnclose(self, read, splitread):
    """
    Create translocation of parts which enclose mate on left side
    """
    info = {}
    pos = splitread.right.pos - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.left.pos - 1
    info['traend'] = splitread.left.end
    info['tracpos'] = (splitread.right.pos - read.end) - self._sample.getMaxInsertSize()
    info['intervals'] = [[splitread.right.pos - 1, splitread.right.pos]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def rightTranslocationEnclose(self, read, splitread):
    """
    Create translocation of parts which enclose mate on right side
    """
    info = {}
    pos = splitread.left.end - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.right.pos - 1
    info['traend'] = splitread.right.end
    info['tracend'] = self._sample.getMaxInsertSize() - (read.pos - splitread.left.end)
    info['intervals'] = [[splitread.left.end, splitread.left.end + 1]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def leftTranslocationRearrangement(self, splitread):
    """
    Create rearranged translocation on left side
    """
    info = {}
    pos = splitread.right.end - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.left.pos - 1
    info['traend'] = splitread.left.end
    info['tracend'] = splitread.right.pos - splitread.left.end
    info['intervals'] = [[splitread.right.end, splitread.right.end + 1]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def rightTranslocationRearrangement(self, splitread):
    """
    Create rearranged translocation on right side
    """
    info = {}
    pos = splitread.left.pos - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.right.pos - 1
    info['traend'] = splitread.right.end
    info['tracpos'] = splitread.left.end - splitread.right.pos
    info['intervals'] = [[splitread.left.pos - 1, splitread.left.pos]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def leftTranslocation(self, splitpair, splitread):
    """
    Create translocation on left side
    """
    info = {}
    pos = splitread.right.pos - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.left.pos - 1
    info['traend'] = splitread.left.end
    info['tracpos'] = splitpair.read.end - splitread.left.pos if splitpair.isReadFirst() else -info['trapos']
    info['intervals'] = [[splitread.right.pos - 1, splitread.right.pos]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def rightTranslocation(self, splitpair, splitread):
    """
    Create translocation on right side
    """
    info = {}
    pos = splitread.left.end - 1
    info['trachrom'] = splitread.original.reference
    info['trapos'] = splitread.right.pos - 1
    info['traend'] = splitread.right.end
    info['tracend'] = self._countMaxCend(splitread.right.tid, info['traend']) if splitpair.isReadFirst() else splitpair.read.pos - info['traend']
    info['intervals'] = [[splitread.left.end, splitread.left.end + 1]]
    return self.__variation(Variation.vtype.TRA, splitread, pos, info)

  def gap(self, paired):
    """
    Create deletion, duplication and translocation
    """
    var1 = self.deletion(paired.splitpair.splitread)
    var2 = self.leftDuplication(paired.splitpair.splitread)
    var3 = self.rightDuplication(paired.splitpair.splitread)
    var4 = self.leftTranslocation(paired.splitpair, paired.splitpair.splitread)
    var5 = self.rightTranslocation(paired.splitpair, paired.splitpair.splitread)
    return OppositeCluster(self._sample, var1, var2, var3, var4, var5)

  def rearrangement(self, paired):
    """
    Create rearranged translocations and duplications
    """
    var1 = self.leftDuplicationRearrangement(paired.splitpair.splitread)
    var2 = self.rightDuplicationRearrangement(paired.splitpair.splitread)
    var3 = self.leftTranslocationRearrangement(paired.splitpair.splitread)
    var4 = self.rightTranslocationRearrangement(paired.splitpair.splitread)
    return OppositeCluster(self._sample, var1, var2, var3, var4)

  def enclose(self, paired):
    """
    Create enclosed translocations and duplications
    """
    if paired.splitpair.isReadFirst():
      var1 = self.leftDuplicationEnclose(paired.splitpair.splitread)
      var2 = self.leftTranslocationEnclose(paired.splitpair.read, paired.splitpair.splitread)
    else:
      var1 = self.rightDuplicationEnclose(paired.splitpair.splitread)
      var2 = self.rightTranslocationEnclose(paired.splitpair.read, paired.splitpair.splitread)

    return OppositeCluster(self._sample, var1, var2)

  def encloseRearrangement(self, paired):
    """
    Create rearranged and enclosed translocations and duplications
    """
    if paired.splitpair.isReadFirst():
      var1 = self.leftDuplicationEncloseRearrangement(paired.splitpair.splitread)
      var2 = self.leftTranslocationEncloseRearrangement(paired.splitpair.read, paired.splitpair.splitread)
    else:
      var1 = self.rightDuplicationEncloseRearrangement(paired.splitpair.splitread)
      var2 = self.rightTranslocationEncloseRearrangement(paired.splitpair.read, paired.splitpair.splitread)

    return OppositeCluster(self._sample, var1, var2)
