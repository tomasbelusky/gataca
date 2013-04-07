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
    info = {'imprecise': True}
    pos = paired.mate().pos - 1
    info['cilen'] = [self._sample.getMinInsertSize() - paired.actualSize(), self._sample.getMaxInsertSize() - paired.actualSize()]
    info['cpos'] = -paired.actualSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INS, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def deletion(self, paired):
    """
    Create deletion
    """
    info = {'imprecise': True}
    pos = paired.readEnd() + self._sample.getMaxInsertSize() - 1
    info['max'] = paired.mate().pos - 1
    info['cilen'] = [paired.actualSize() - self._sample.getMaxInsertSize(), paired.actualSize() - self._sample.getMinInsertSize()]
    info['cpos'] = -self._sample.getMaxInsertSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DEL, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def inversionRead(self, paired):
    """
    Create inversion of read
    """
    info = {'imprecise': True}
    pos = paired.read().pos
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = paired.actualSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), info=info)

  def inversionMate(self, paired):
    """
    Create inversion of mate
    """
    info = {'imprecise': True}
    pos = paired.mate().pos
    info['end'] = paired.mateEnd()
    info['cpos'] = -paired.actualSize()
    info['cend'] = self._sample.getMaxInsertSize()
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.mate(), info=info)

  def inversionReadOnly(self, read):
    """
    Create inversion of one read not according to it's mate
    """
    info = {'imprecise': True}
    pos = read.pos
    info['end'] = read.pos + read.rlen
    info['cpos'] = "-"
    info['cend'] = "+"
    refseq = self._sample.fetchReference(read.tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.INV, self._sample.getRefName(read.tid), pos, None, refseq, Variation.mtype.READ_PAIR, read, info=info)

  def overlap(self, paired):
    """
    Create duplication from overlap
    """
    info = {'imprecise': True}
    pos = paired.mate().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = self._sample.getMaxInsertSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    var = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var.incDepth()
    return var

    """
    INFO: possible feature
    # first possible position
    info['dupchrom'] = self._sample.getRefName(paired.mate().tid)
    info['duppos'] = paired.readEnd() + self._sample.getMaxInsertSize()
    info['cduppos'] = -self._sample.getMaxInsertSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self.repairInfo(pos, info)
    var1 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var1.incDepth()

    # second possible position
    info['dupchrom'] = self._sample.getRefName(paired.read().tid)
    info['duppos'] = paired.mate().pos - self._sample.getMaxInsertSize() + paired.actualSize()
    info['cduppos'] = self._sample.getMaxInsertSize() - self._sample.getMinInsertSize()
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self.repairInfo(pos, info)
    var2 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var2.incDepth()

    return OppositeCluster(var1, var2)
    """

  def overlapRearranged(self, paired):
    """
    Create duplication from rearranged overlap
    """
    info = {'imprecise': True}
    pos = paired.mate().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -(paired.mateEnd() - paired.read().pos + paired.actualSize() + self._sample.getMaxInsertSize())
    info['cend'] = -info['cpos']
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    var = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var.incDepth()
    return var

    """
    INFO: possible feature
    # first possible position
    info['dupchrom'] = self._sample.getRefName(paired.read().tid)
    info['duppos'] = paired.mateEnd() + self._sample.getMinInsertSize()
    info['cduppos'] = [0, self._sample.getMaxInsertSize() - self._sample.getMinInsertSize()]
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self.repairInfo(pos, info)
    var1 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var1.incDepth()

    # second possible position
    info['dupchrom'] = self._sample.getRefName(paired.mate().tid)
    info['duppos'] = paired.read().pos - self._sample.getMinInsertSize()
    info['cduppos'] = [self._sample.getMinInsertSize() - self._sample.getMaxInsertSize(), 0]
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self.repairInfo(pos, info)
    var2 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var2.incDepth()

    return OppositeCluster(var1, var2)
    """

  def rightDuplication(self, paired):
    """
    Create duplication on right side
    """
    info = {'imprecise' : True}
    pos = paired.mate().pos - 1
    info['end'] = paired.mateEnd()
    info['cpos'] = -paired.actualSize()
    info['cend'] = self._sample.getMaxInsertSize()
    """
    INFO: possible feature
    info['dupchrom'] = self._sample.getRefName(paired.read().tid)
    info['duppos'] = paired.read().pos - self._sample.getMinInsertSize() - (paired.mateEnd() - paired.mate().pos)
    info['cduppos'] = [-paired.actualSize() + self._sample.getMinInsertSize() - self._sample.getMaxInsertSize(), 0]
    """
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def leftDuplication(self, paired):
    """
    Create duplication on left side
    """
    info = {'imprecise' : True}
    pos = paired.read().pos - 1
    info['end'] = paired.readEnd()
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['cend'] = paired.actualSize()
    """
    INFO: possible feature
    info['dupchrom'] = self._sample.getRefName(paired.mate().tid)
    info['duppos'] = paired.mateEnd() + self._sample.getMinInsertSize()
    info['cduppos'] = [0, self._sample.getMaxInsertSize() - self._sample.getMinInsertSize()]
    """
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.DUP, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def rightTranslocation(self, paired):
    """
    Create translocation on right side
    """
    info = {'imprecise': True}
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

    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def leftTranslocation(self, paired):
    """
    Create translocation on left side
    """
    info = {'imprecise': True}
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

    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def rightTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on right side
    """
    info = {'imprecise': True}
    pos = paired.read().pos - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.mate().tid)
    info['trapos'] = paired.mate().pos
    info['traend'] = paired.mateEnd()
    info['tracpos'] = -paired.actualSize()
    info['tracend'] = self._sample.getMaxInsertSize()
    refseq = self._sample.fetchReference(paired.read().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def leftTranslocationRearranged(self, paired):
    """
    Create rearranged translocation on left side
    """
    info = {'imprecise': True}
    pos = paired.mateEnd() + self._sample.getMaxInsertSize() - 1
    info['cpos'] = -self._sample.getMaxInsertSize()
    info['trachrom'] = self._sample.getRefName(paired.read().tid)
    info['trapos'] = paired.read().pos
    info['traend'] = paired.readEnd()
    info['tracpos'] = -self._sample.getMaxInsertSize()
    info['tracend'] = paired.actualSize()
    refseq = self._sample.fetchReference(paired.mate().tid, pos, pos+1)
    self._repairInfo(pos, info)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  def translocation(self, paired):
    """
    Create translocations
    """
    var1 = self.leftTranslocation(paired)
    var2 = self.rightTranslocation(paired)
    return OppositeCluster(var1, var2)

  def translocationRearranged(self, paired):
    """
    Create rearranged translocations
    """
    var1 = self.leftTranslocationRearranged(paired)
    var2 = self.rightTranslocationRearranged(paired)
    return OppositeCluster(var1, var2)

  def rearrangement(self, paired):
    """
    Create rearranged translocations and duplications
    """
    var1 = self.leftTranslocationRearranged(paired)
    var2 = self.rightTranslocationRearranged(paired)
    var3 = self.leftDuplication(paired)
    var4 = self.rightDuplication(paired)
    return OppositeCluster(var1, var2, var3, var4)
