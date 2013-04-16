#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "25.03. 2013"

from resources.reads.Read import Read
from variations.Variation import Variation
from variations.clusters.OppositeCluster import OppositeCluster
from BaseFactory import BaseFactory

class PairFactory(BaseFactory):
  """
  Creating factory of variations detected by read pair method
  """

  def __init__(self, sample):
    """
    Initialize variables
    """
    BaseFactory.__init__(self, sample)

  def __variation(self, vtype, pos, info, tid, reference):
    """
    Create variation
    """
    refseq = self._sample.fetchReference(tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(vtype, reference, pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def insertion(self, paired):
    """
    Create insertion
    """
    info = {}
    pos = paired.mate.pos - 1
    info['cilen'] = [self._sample.getMinInsertSize() - paired.actualSize(), self._sample.getMaxInsertSize() - paired.actualSize()]
    info['cpos'] = -paired.actualSize()
    info['intervals'] = [[paired.read.end, paired.mate.pos], [paired.read.end, paired.mate.pos]]
    return self.__variation(Variation.vtype.INS, pos, info, paired.read.tid, paired.read.reference)

  def deletion(self, paired):
    """
    Create deletion
    """
    info = {}
    pos = paired.read.end + self._sample.getMaxInsertSize() - 1
    info['end'] = paired.mate.pos - 1
    info['cilen'] = [paired.actualSize() - self._sample.getMaxInsertSize(), paired.actualSize() - self._sample.getMinInsertSize()]
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read.end, paired.mate.pos], [paired.read.end, paired.mate.pos]]
    return self.__variation(Variation.vtype.DEL, pos, info, paired.read.tid, paired.read.reference)

  def inversionRead(self, paired):
    """
    Create inversion of read
    """
    info = {}
    pos = paired.read.pos
    info['end'] = paired.read.end
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = paired.actualSize()
    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.INV, pos, info, paired.read.tid, paired.read.reference)

  def inversionMate(self, paired):
    """
    Create inversion of mate
    """
    info = {}
    pos = paired.mate.pos
    info['end'] = paired.mate.end
    info['cpos'] = -paired.actualSize()
    info['cend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.INV, pos, info, paired.mate.tid, paired.read.reference)

  def inversionReadOnly(self, read):
    """
    Create inversion of one read not according to it's mate
    """
    info = {}
    pos = read.pos
    info['end'] = Read.calculateEnd(read)
    info['cpos'] = -pos
    info['cend'] = self._countMaxLength(read.tid, info['end'])
    info['intervals'] = [[pos, info['end']]]
    return self.__variation(Variation.vtype.INV, pos, info, read.tid, self._sample.getRefName(read.tid))

  def overlap(self, paired):
    """
    Create duplication from overlap
    """
    info = {}
    pos = paired.mate.pos - 1
    info['end'] = min(paired.read.end, paired.mate.end)
    info['cpos'] = -self._sample.getMaxInsertSize()
    plusEnd = paired.read.end - paired.mate.end if info['end'] == paired.mate.end else 0
    info['cend'] = self._sample.getMaxInsertSize() + plusEnd
    #smallestSize = info['end'] - pos
    #info['cilen'] = [smallestSize, smallestSize + info['cend']]
    info['intervals'] = [[pos, info['end']], [pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.read.tid, paired.read.reference)

  def overlapRearranged(self, paired):
    """
    Create duplication from rearranged overlap
    """
    info = {}
    pos = paired.mate.pos - 1
    info['end'] = min(paired.read.end, paired.mate.end)
    info['cpos'] = -(paired.mate.end - paired.read.pos + paired.actualSize() + self._sample.getMaxInsertSize())
    plusEnd = paired.read.end - paired.mate.end if info['end'] == paired.mate.end else 0
    info['cend'] = -info['cpos'] + plusEnd
    info['intervals'] = [[pos, info['end']], [pos, info['end']]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.mate.tid, paired.read.reference)

  def rightSmallDuplication(self, paired):
    """
    Create duplication on right side when insert size is small
    """
    info = {}
    pos = paired.mate.pos - 1
    overlap = ((paired.actualSize() + paired.mate.len) - self._sample.getMinInsertSize()) > 0
    info['end'] = paired.read.end + self._sample.getMinInsertSize() if overlap else paired.mate.end
    info['cend'] = paired.read.end + self._sample.getMaxInsertSize() - info['end']
    info['cpos'] = -info['cend']
    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.read.tid, paired.read.reference)

  def leftSmallDuplication(self, paired):
    """
    Create duplication on right side when insert size is small
    """
    info = {}
    info['end'] = paired.read.end
    overlap = ((paired.actualSize() + paired.read.len) - self._sample.getMinInsertSize()) > 0
    pos = (paired.mate.pos - self._sample.getMinInsertSize() if overlap else paired.read.pos) - 1
    info['cpos'] = paired.mate.pos - self._sample.getMaxInsertSize() - pos
    info['cend'] = -info['cpos']
    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.mate.tid, paired.mate.reference)

  def rightBigDuplication(self, paired):
    """
    Create duplication on right side when insert size is big
    """
    info = {}
    info['end'] = paired.mate.end
    overlap = (self._sample.getMaxInsertSize() + paired.mate.len - paired.actualSize()) < 0
    pos = (paired.read.end + self._sample.getMaxInsertSize() + paired.mate.len if overlap else paired.mate.pos) - 1
    info['cpos'] = paired.read.end + self._sample.getMinInsertSize() + paired.mate.len - pos
    info['cend'] = -info['cpos']
    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.read.tid, paired.read.reference)

  def leftBigDuplication(self, paired):
    """
    Create duplication on right side when insert size is big
    """
    info = {}
    pos = paired.read.pos - 1
    overlap = (paired.actualSize() - (self._sample.getMaxInsertSize() + paired.read.len)) < 0
    info['end'] = paired.mate.pos - self._sample.getMaxInsertSize() - paired.read.len if overlap else paired.read.end
    info['cend'] = paired.mate.pos - self._sample.getMinInsertSize() - paired.read.len - info['end']
    info['cpos'] = -info['cend']
    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.mate.tid, paired.mate.reference)

  def rightDuplicationRearranged(self, paired):
    """
    Create duplication on right side with rearranged reads
    """
    info = {}
    pos = paired.mate.pos - 1
    info['end'] = paired.mate.end
    info['cpos'] = paired.actualSize() + paired.mate.len - self._sample.getMaxInsertSize()
    info['cend'] = -info['cpos']
    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.read.tid, paired.read.reference)

  def leftDuplicationRearranged(self, paired):
    """
    Create duplication on left side with rearranged reads
    """
    info = {}
    pos = paired.read.pos - 1
    info['end'] = paired.read.end
    info['cpos'] = paired.actualSize() + paired.read.len - self._sample.getMaxInsertSize()
    info['cend'] = -paired.actualSize()
    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.DUP, pos, info, paired.mate.tid, paired.mate.reference)

  def rightTranslocation(self, paired):
    """
    Create translocation on right side
    """
    info = {}
    pos = paired.read.end + self._sample.getMaxInsertSize() - 1
    info['cpos'] = -min(self._sample.getMaxInsertSize(), paired.actualSize())
    info['trachrom'] = paired.mate.reference
    info['trapos'] = paired.mate.pos
    info['traend'] = paired.mate.end

    if paired.mate.isInverted():
      info['tracpos'] = -info['trapos']
      info['tracend'] = self._sample.getMaxInsertSize()
    else:
      info['tracpos'] = -self._sample.getMaxInsertSize()
      info['tracend'] = self._countMaxLength(paired.mate.tid, info['traend'])

    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.TRA, pos, info, paired.read.tid, paired.read.reference)

  def leftTranslocation(self, paired):
    """
    Create translocation on left side
    """
    info = {}
    pos = paired.mate.pos - 1
    info['cpos'] = -min(self._sample.getMaxInsertSize(), paired.actualSize())
    info['trachrom'] = paired.read.reference
    info['trapos'] = paired.read.pos
    info['traend'] = paired.read.end

    if paired.read.isInverted():
      info['tracpos'] = -self._sample.getMaxInsertSize()
      info['tracend'] = self._countMaxLength(paired.read.tid, info['traend'])
    else:
      info['tracpos'] = -info['trapos']
      info['tracend'] = self._sample.getMaxInsertSize()

    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.TRA, pos, info, paired.mate.tid, paired.mate.reference)

  def rightTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on right side
    """
    info = {}
    pos = paired.read.pos - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = paired.mate.reference
    info['trapos'] = paired.mate.pos
    info['traend'] = paired.mate.end
    info['tracpos'] = paired.actualSize()
    info['tracend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.mate.pos, paired.mate.end]]
    return self.__variation(Variation.vtype.TRA, pos, info, paired.read.tid, paired.read.reference)

  def leftTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on left side
    """
    info = {}
    pos = paired.mate.end + self._sample.getMaxInsertSize() - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = paired.read.reference
    info['trapos'] = paired.read.pos
    info['traend'] = paired.read.end
    info['tracpos'] = -self._sample.getMaxInsertSize()
    info['tracend'] = -paired.actualSize()
    info['intervals'] = [[paired.read.pos, paired.read.end]]
    return self.__variation(Variation.vtype.TRA, pos, info, paired.mate.tid, paired.mate.reference)

  def translocation(self, paired):
    """
    Create translocations
    """
    var1 = self.leftTranslocation(paired)
    var2 = self.rightTranslocation(paired)
    return OppositeCluster(self._sample, var1, var2)

  def translocationRearranged(self, paired):
    """
    Create rearranged translocations
    """
    var1 = self.leftTranslocationRearranged(paired)
    var2 = self.rightTranslocationRearranged(paired)
    return OppositeCluster(self._sample, var1, var2)

  def bigInsertSize(self, paired):
    """
    Create deletion, translocations and duplications
    """
    var1 = self.deletion(paired)
    var2 = self.leftTranslocation(paired)
    var3 = self.rightTranslocation(paired)
    var4 = self.leftBigDuplication(paired)
    var5 = self.rightBigDuplication(paired)
    return OppositeCluster(self._sample, var1, var2, var3, var4, var5)

  def smallInsertSize(self, paired):
    """
    Create deletion, translocations and duplications
    """
    var1 = self.insertion(paired)
    var2 = self.leftTranslocation(paired)
    var3 = self.rightTranslocation(paired)
    var4 = self.leftSmallDuplication(paired)
    var5 = self.rightSmallDuplication(paired)
    return OppositeCluster(self._sample, var1, var2, var3, var4, var5)

  def rearrangement(self, paired):
    """
    Create rearranged translocations and duplications
    """
    var1 = self.leftTranslocationRearranged(paired)
    var2 = self.rightTranslocationRearranged(paired)
    var3 = self.leftDuplicationRearranged(paired)
    var4 = self.rightDuplicationRearranged(paired)
    return OppositeCluster(self._sample, var1, var2, var3, var4)
