#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

_author_ = "Tomáš Beluský"
_date_ = "13.03. 2013"

from variations.Variation import Variation

class BreakendCluster:
  """
  Represents cluster with variations which are represents with breakends in VCF format
    INFO: not using -> delete, maybe?
  """
  __counter = {Variation.vtype.TRA : 0,
               Variation.vtype.DUP : 0}
  __svtype = {Variation.vtype.TRA : 'TRA',
              Variation.vtype.DUP : 'DUP'}

  def __init__(self, reference, sample):
    """
    Initialize variables
    """
    AbstractCluster.__init__(self, reference, sample)
    self._consensus = None
    self._copyPos = 0
    self._copyReference = 0

  def _process(self):
    """
    Join variations together for printing them in VCF output
    """
    self._consensus = self._variations[0]
    self._start = self._consensus.getStart()
    self._end = self._consensus.getEnd()
    self._copyPos = self._consensus.getCopyPosition()
    self._copyReference = self._consensus.getCopyReference()

  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    """
    return False

  def __nextEvent(self, vtype):
    event = '%s%s' % (BreakendCluster.__svtype[vtype], BreakendCluster.__counter[vtype])
    BreakendCluster.__counter[vtype] += 1
    return event

  def __str__(self):
    """
    Print cluster in VCF format
    """
    if not self._consensus:
      return ""

    result = ""
    info = self._consensus.getInfo()
    info['event'] = self.__nextEvent(self._consensus.getType())
    cicopypos = info.get('cicopypos', None)
    cipos = info.get('cipos', None)
    ciend = info.get('ciend', None)
    cilen = info.get('cilen', None)
    del info['copypos']
    del info['cicopypos']
    del info['copyreference']

    if self._consensus.getType() == Variation.vtype.TRA: # translocation
      info['fulldepth'] = self._sample.getExactCoverages(self._rindex, self._start, self._end)
      info['svtype'] = 'DEL'
      result += "%s\t%s\t.\t%s\t%s\t.\t.\t%s\n" % (self._rname,
                                                   self._start,
                                                   self._consensus.getReferenceSequence(),
                                                   self._consensus.getSequence(),
                                                   self.infoString(info))

    # first record of breakend -------------------------------------------------
    del info['cipos']
    del info['ciend']
    del info['cilen']

    if cicopypos:
      info['cipos'] = cicopypos
    if cipos:
      info['cialtpos'] = cipos
    if cilen:
      info['cialtlen'] = cilen

    info['fulldepth'] = self._sample.getExactCoverage(self._rindex, self._copyPos)
    info['svtype'] = 'BND'
    refseq = self._sample.fetchReference(self._rindex, self._copyPos, self._copyPos + 2)
    result += "%s\t%s\t.\t%s\t%s[%s:%s[\t.\t.\t%s\n" % (self._copyReference,
                                                        self._copyPos,
                                                        refseq[0],
                                                        refseq[0],
                                                        self._rname,
                                                        self._start + 1,
                                                        self.infoString(info))

    # second record of breakend ------------------------------------------------
    if ciend:
      info['cialtpos'] = ciend
    else:
      del info['cialtpos']

    result += "%s\t%s\t.\t%s\t]%s:%s]%s\t.\t.\t%s" % (self._copyReference,
                                                      self._copyPos+1,
                                                      refseq[1],
                                                      self._rname,
                                                      self._end,
                                                      refseq[1],
                                                      self.infoString(info))
    return result
