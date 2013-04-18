#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "07.04. 2013"

class BaseFactory:
  """
  Factory with base operations with variations
  """

  def __init__(self, sample):
    """
    Initialize variables
    """
    self._sample = sample

  def _countMaxLength(self, tid, end):
    """
    Count maximal confidence for end position
    """
    return self._sample.getLengths()[tid] - end

  def _repairInfo(self, start, info):
    """
    Repair informations about variations (remove unused, count imprecise)
    """
    if info.get('cpos', -1) >= 0: # remove cpos
      del info['cpos']

    if info.get('cend', 1) <= 0: # remove cend
      del info['cend']

    if (start + info.get('cpos', 0)) < 0: # can't be negative
      info['cpos'] = -start

    if info.get('tracpos', -1) >= 0: # remove tracpos
      del info['tracpos']

    if info.get('tracend', 1) <= 0: # remove tracend
      del info['tracend']

    cilen = info.get('cilen', [-1, 1])

    if cilen[0] == cilen[1]: # remove cilen
      if cilen[0] != 0: # create svlen
        info['svlen'] = cilen[0]

      del info['cilen']

    if info.get('cpos', None) is None and info.get('max', -1) == start: # remove max
      del info['max']

    # count imprecise
    info['imprecise'] = 'cpos' in info or \
                        'cend' in info or \
                        'cilen' in info or \
                        'max' in info or \
                        'tracpos' in info or \
                        'tracend' in info
