#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "25.03. 2013"

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

  def insertion(self, paired):
    """
    Create insertion
    """
    info = {}
    pos = paired.mate().pos - 1
    info['cilen'] = [self._sample.getMinInsertSize() - paired.actualSize(), self._sample.getMaxInsertSize() - paired.actualSize()]
    info['cpos'] = -paired.actualSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INS, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def deletion(self, paired):
    """
    Create deletion
    """
    info = {}
    pos = paired.readEnd() + self._sample.getMaxInsertSize() - 1
    info['max'] = paired.mate().pos - 1
    info['cilen'] = [paired.actualSize() - self._sample.getMaxInsertSize(), paired.actualSize() - self._sample.getMinInsertSize()]
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DEL, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def inversionRead(self, paired):
    """
    Create inversion of read
    """
    info = {}
    pos = paired.read().pos
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = paired.actualSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def inversionMate(self, paired):
    """
    Create inversion of mate
    """
    info = {}
    pos = paired.mate().pos
    info['end'] = paired.mateEnd()
    info['cpos'] = -paired.actualSize()
    info['cend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def inversionReadOnly(self, read):
    """
    Create inversion of one read not according to it's mate
    """
    info = {}
    pos = read.pos
    info['end'] = read.pos + read.rlen
    info['cpos'] = "-"
    info['cend'] = "+"
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(read.tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, self._sample.getRefName(read.tid), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def overlap(self, paired):
    """
    Create duplication from overlap
    """
    info = {}
    pos = paired.mate().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def overlapRearranged(self, paired):
    """
    Create duplication from rearranged overlap
    """
    info = {}
    pos = paired.mate().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -(paired.mateEnd() - paired.read().pos + paired.actualSize() + self._sample.getMaxInsertSize())
    info['cend'] = -info['cpos']
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def rightDuplication(self, paired):
    """
    Create duplication on right side
    """
    info = {}
    pos = paired.mate().pos - 1
    info['end'] = paired.mateEnd()
    info['cpos'] = -paired.actualSize()
    info['cend'] = "+"
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def leftDuplication(self, paired):
    """
    Create duplication on right side
    """
    info = {}
    pos = paired.read().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = "-"
    info['cend'] = paired.actualSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def rightDuplicationRearranged(self, paired):
    """
    Create duplication on right side with rearranged reads
    """
    info = {}
    pos = paired.mate().pos - 1
    info['end'] = paired.mateEnd()
    info['cpos'] = paired.actualSize()
    info['cend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def leftDuplicationRearranged(self, paired):
    """
    Create duplication on left side with rearranged reads
    """
    info = {}
    pos = paired.read().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = -paired.actualSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def rightTranslocation(self, paired):
    """
    Create translocation on right side
    """
    info = {}
    pos = paired.readEnd() + self._sample.getMaxInsertSize() - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.mate().tid)
    info['trapos'] = paired.mate().pos
    info['traend'] = paired.mateEnd()

    if paired.isMateInverted():
      info['tracpos'] = "-"
      info['tracend'] = self._sample.getMaxInsertSize()
    else:
      info['tracpos'] = -self._sample.getMaxInsertSize()
      info['tracend'] = "+"

    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def leftTranslocation(self, paired):
    """
    Create translocation on left side
    """
    info = {}
    pos = paired.mate().pos - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.read().tid)
    info['trapos'] = paired.read().pos
    info['traend'] = paired.readEnd()

    if paired.isReadInverted():
      info['tracpos'] = -self._sample.getMaxInsertSize()
      info['tracend'] = "+"
    else:
      info['tracpos'] = "-"
      info['tracend'] = self._sample.getMaxInsertSize()

    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def rightTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on right side
    """
    info = {}
    pos = paired.read().pos - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.mate().tid)
    info['trapos'] = paired.mate().pos
    info['traend'] = paired.mateEnd()
    info['tracpos'] = paired.actualSize()
    info['tracend'] = self._sample.getMaxInsertSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

  def leftTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on left side
    """
    info = {}
    pos = paired.mateEnd() + self._sample.getMaxInsertSize() - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.read().tid)
    info['trapos'] = paired.read().pos
    info['traend'] = paired.readEnd()
    info['tracpos'] = -self._sample.getMaxInsertSize()
    info['tracend'] = -paired.actualSize()
    info['intervals'] = [[paired.read().pos, paired.readEnd()], [paired.mate().pos, paired.mateEnd()]]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, info=info)

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

  def duplication(self, paired):
    """
    Create duplications
    """
    var1 = self.leftDuplication(paired)
    var2 = self.rightDuplication(paired)
    return OppositeCluster(self._sample, var1, var2)

  def duplicationRearranged(self, paired):
    """
    Create rearranged duplications
    """
    var1 = self.leftDuplicationRearranged(paired)
    var2 = self.rightDuplicationRearranged(paired)
    return OppositeCluster(self._sample, var1, var2)

  def bigInsertSize(self, paired):
    """
    Create deletion, translocations and duplications
    """
    var1 = self.deletion(paired)
    var2 = self.leftTranslocation(paired)
    var3 = self.rightTranslocation(paired)
    var4 = self.leftDuplication(paired)
    var5 = self.rightDuplication(paired)
    return OppositeCluster(self._sample, var1, var2, var3, var4, var5)

  def smallInsertSize(self, paired):
    """
    Create insertion, translocations and duplications
    """
    var1 = self.insertion(paired)
    var2 = self.leftTranslocation(paired)
    var3 = self.rightTranslocation(paired)
    var4 = self.leftDuplication(paired)
    var5 = self.rightDuplication(paired)
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
