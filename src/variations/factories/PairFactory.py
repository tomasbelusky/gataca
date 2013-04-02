#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "25.03. 2013"

from variations.Variation import Variation
from variations.clusters.OppositeCluster import OppositeCluster

class PairFactory:
  """
  Creating factory of variations detected by read pair method
  """

  @staticmethod
  def insertion(paired, sample):
    info = {'imprecise': True}
    pos = paired.readEnd() - 1
    info['end'] = pos
    info['cilen'] = [sample.getMinInsertSize() - paired.actualSize(), sample.getMaxInsertSize() - paired.actualSize()]
    info['cipos'] = [0, paired.mate().pos - pos]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.INS, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def deletion(paired, sample):
    info = {'imprecise': True}
    pos = paired.readEnd() - 1
    info['end'] = paired.mate().pos
    info['cilen'] = [-(paired.actualSize() - sample.getMinInsertSize()), -(paired.actualSize() - sample.getMaxInsertSize())]
    info['cipos'] = [0, paired.mate().pos - pos + info['cilen'][1]]
    info['ciend'] = [-info['cipos'][1], 0]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.DEL, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def inversionRead(paired, sample):
    info = {'imprecise': True}
    pos = paired.read().pos - sample.getMaxInsertSize()
    info['end'] = paired.mate().pos
    info['cilen'] = [paired.readEnd() - paired.read().pos, info['end'] - pos]
    info['cipos'] = [0, info['cilen'][1] - info['cilen'][0]]
    info['ciend'] = [-info['cipos'][1], 0]
    addToLen = paired.actualSize() - sample.getMaxInsertSize()
    toRead = paired.read().pos - pos
    afterRead = -(info['end'] - paired.readEnd())

    if addToLen > 0:
      info['cilen'][0] += addToLen

    if toRead < info['cipos'][1]:
      info['cipos'][1] = toRead

    if afterRead > info['ciend'][0]:
      info['ciend'][0] = afterRead

    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), info=info)

  @staticmethod
  def inversionMate(paired, sample):
    info = {'imprecise': True}
    pos = paired.readEnd()
    info['end'] = paired.mateEnd() + sample.getMaxInsertSize()
    info['cilen'] = [paired.mateEnd() - paired.mate().pos, info['end'] - pos]
    info['cipos'] = [0, info['cilen'][1] - info['cilen'][0]]
    info['ciend'] = [-info['cipos'][1], 0]
    addToLen = paired.actualSize() - sample.getMaxInsertSize()
    toMate = paired.mate().pos - pos
    afterMate = -(info['end'] - paired.mateEnd())

    if addToLen > 0:
      info['cilen'][0] += addToLen

    if toMate < info['cipos'][1]:
      info['cipos'][1] = toMate

    if afterMate > info['ciend'][0]:
      info['ciend'][0] = afterMate

    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    return Variation(Variation.vtype.INV, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.mate(), info=info)

  @staticmethod
  def inversionReadOnly(read, sample):
    info = {'imprecise': True}
    pos = read.pos
    info['end'] = read.pos + read.alen
    info['cilen'] = [info['end'] - pos, "."]
    info['cipos'] = [".", 0]
    info['ciend'] = [0, "."]
    refseq = sample.fetchReference(read.tid, pos, pos+1)
    return Variation(Variation.vtype.INV, sample.getRefName(read.tid), pos, None, refseq, Variation.mtype.READ_PAIR, read, info=info)

  @staticmethod
  def overlap(paired, sample):
    info = {'imprecise': True}
    pos = paired.mate().pos - sample.getMaxInsertSize() - 1
    info['end'] = paired.readEnd() + sample.getMaxInsertSize()
    info['cilen'] = [paired.readEnd() - paired.mate().pos, info['end'] - pos]
    info['cipos'] = [0, paired.mate().pos - pos]
    info['ciend'] = [paired.readEnd() - info['end'], 0]

    # first possible position
    info['dupchrom'] = sample.getRefName(paired.mate().tid)
    info['duppos'] = paired.readEnd() + sample.getMinInsertSize()
    info['ciduppos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    var1 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var1.incDepth()

    # second possible position
    info['dupchrom'] = sample.getRefName(paired.read().tid)
    info['duppos'] = paired.mate().pos - sample.getMinInsertSize()
    info['ciduppos'] = [sample.getMinInsertSize() - sample.getMaxInsertSize(), 0]
    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    var2 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var2.incDepth()

    return OppositeCluster(var1, var2)

  @staticmethod
  def overlapRearranged(paired, sample):
    info = {'imprecise': True}
    pos = paired.read().pos - sample.getMaxInsertSize() - 1
    info['end'] = paired.mateEnd() + sample.getMaxInsertSize()
    info['cilen'] = [paired.readEnd() - paired.mate().pos, info['end'] - pos]
    info['cipos'] = [0, paired.mate().pos - pos]
    info['ciend'] = [paired.readEnd() - info['end'], 0]

    # first possible position
    info['dupchrom'] = sample.getRefName(paired.read().tid)
    info['duppos'] = paired.mateEnd() + sample.getMinInsertSize()
    info['ciduppos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]
    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    var1 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var1.incDepth()

    # second possible position
    info['dupchrom'] = sample.getRefName(paired.mate().tid)
    info['duppos'] = paired.read().pos - sample.getMinInsertSize()
    info['ciduppos'] = [sample.getMinInsertSize() - sample.getMaxInsertSize(), 0]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    var2 = Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)
    var2.incDepth()

    return OppositeCluster(var1, var2)

  @staticmethod
  def rightDuplication(paired, sample):
    info = {'imprecise' : True}
    pos = paired.mate().pos - 1
    info['end'] = paired.mateEnd()
    info['cilen'] = [paired.mateEnd() - paired.mate().pos, paired.mateEnd() - paired.readEnd() + sample.getMaxInsertSize()]
    info['cipos'] = [paired.readEnd() - paired.mate().pos, 0]
    info['ciend'] = [0, sample.getMaxInsertSize()]
    info['dupchrom'] = sample.getRefName(paired.read().tid)
    info['duppos'] = paired.read().pos - sample.getMinInsertSize() - (paired.mateEnd() - paired.mate().pos) - 1
    info['ciduppos'] = [paired.readEnd() - paired.mate().pos + sample.getMinInsertSize() - sample.getMaxInsertSize(), 0]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.DUP, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def leftDuplication(paired, sample):
    info = {'imprecise' : True}
    pos = paired.read().pos - 1
    info['end'] = paired.readEnd()
    info['cilen'] = [paired.readEnd() - paired.read().pos, paired.mate().pos - paired.read().pos + sample.getMaxInsertSize()]
    info['cipos'] = [-sample.getMaxInsertSize(), 0]
    info['ciend'] = [0, paired.mate().pos - paired.readEnd()]
    info['dupchrom'] = sample.getRefName(paired.read().tid)
    info['duppos'] = paired.mateEnd() + sample.getMinInsertSize()
    info['ciduppos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]
    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    return Variation(Variation.vtype.DUP, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def rearrangement(paired, sample):
    var1 = PairFactory.leftTranslocationRearranged(paired, sample)
    var2 = PairFactory.rightTranslocationRearranged(paired, sample)
    var3 = PairFactory.leftDuplication(paired, sample)
    var4 = PairFactory.rightDuplication(paired, sample)
    return OppositeCluster(var1, var2, var3, var4)

  @staticmethod
  def rightTranslocation(paired, sample):
    info = {'imprecise': True}

    if paired.isReadInverted():
      pos = paired.readEnd()
      info['cipos'] = [0, sample.getMaxInsertSize()]
    else:
      pos = paired.readEnd() + sample.getMinInsertSize() - 1
      info['cipos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]

    if paired.isMateInverted():
      info['trapos'] = paired.mate().pos
      info['tracipos'] = [".", 0]
      info['traciend'] = [0, sample.getMaxInsertSize()]
    else:
      info['trapos'] = paired.mate().pos - sample.getMaxInsertSize()
      info['tracipos'] = [0, sample.getMaxInsertSize()]
      info['traciend'] = [0, "."]

    info['end'] = pos
    info['cilen'] = [paired.mateEnd() - paired.mate().pos, "."]
    info['trachrom'] = sample.getRefName(paired.mate().tid)
    info['traend'] = paired.mateEnd()
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def leftTranslocation(paired, sample):
    info = {'imprecise': True}

    if paired.isReadInverted():
      info['tracipos'] = [-sample.getMaxInsertSize(), 0]
      info['traend'] = paired.readEnd()
      info['traciend'] = [0, "."]
    else:
      info['tracipos'] = [".", 0]
      info['traend'] = paired.readEnd() + sample.getMaxInsertSize()
      info['traciend'] = [sample.getMinInsertSize() - sample.getMaxInsertSize(), 0]

    if paired.isMateInverted():
      info['cipos'] = [0, sample.getMaxInsertSize()]
    else:
      info['cipos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]

    pos = paired.mate().pos - sample.getMaxInsertSize() - (paired.readEnd() - paired.read().pos)
    info['end'] = pos
    info['cilen'] = [paired.readEnd() - paired.read().pos, "."]
    info['trachrom'] = sample.getRefName(paired.read().tid)
    info['trapos'] = paired.read().pos
    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def translocation(paired, sample):
    return OppositeCluster(PairFactory.leftTranslocation(paired, sample), PairFactory.rightTranslocation(paired, sample))

  @staticmethod
  def rightTranslocationRearranged(paired, sample):
    info = {'imprecise': True}
    pos = paired.read().pos - sample.getMaxInsertSize() - (paired.mateEnd() - paired.mate().pos)
    info['end'] = pos
    info['cilen'] = [paired.mateEnd() - paired.mate().pos, "."]
    info['cipos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]
    info['trachrom'] = sample.getRefName(paired.mate().tid)
    info['trapos'] = paired.mate().pos
    info['traend'] = paired.mateEnd() + sample.getMaxInsertSize()
    info['tracipos'] = [".", 0]
    info['traciend'] = [sample.getMinInsertSize() - sample.getMaxInsertSize(), 0]
    refseq = sample.fetchReference(paired.read().tid, pos, pos+1)
    return Variation(Variation.vtype.TRA, paired.getReadReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def leftTranslocationRearranged(paired, sample):
    info = {'imprecise': True}
    pos = paired.mateEnd() + sample.getMinInsertSize() - 1
    info['end'] = pos
    info['cilen'] = [paired.readEnd() - paired.read().pos, "."]
    info['cipos'] = [0, sample.getMaxInsertSize() - sample.getMinInsertSize()]
    info['trachrom'] = sample.getRefName(paired.read().tid)
    info['trapos'] = paired.read().pos - sample.getMaxInsertSize()
    info['traend'] = paired.readEnd()
    info['tracipos'] = [0, sample.getMaxInsertSize()]
    info['traciend'] = [0, "."]
    refseq = sample.fetchReference(paired.mate().tid, pos, pos+1)
    return Variation(Variation.vtype.TRA, paired.getMateReference(), pos, None, refseq, Variation.mtype.READ_PAIR, paired.read(), paired.mate(), info=info)

  @staticmethod
  def translocationRearranged(paired, sample):
    return OppositeCluster(PairFactory.leftTranslocationRearranged(paired, sample), PairFactory.rightTranslocationRearranged(paired, sample))
