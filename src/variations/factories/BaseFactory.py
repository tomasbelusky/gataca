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

    if info.get('svlen', 1) <= 0: # svlen can't be zero or negative
      del info['svlen']

    if (start + info.get('cpos', 0)) < 0: # can't be negative
      info['cpos'] = -start

    if info.get('tracpos', -1) >= 0: # remove tracpos
      del info['tracpos']

    if info.get('tracend', 1) <= 0: # remove tracend
      del info['tracend']

    end = info.get('end', info.get('max', start))
    length = end - start + 1
    maxStart = start + info.get('cpos', 0)
    maxEnd = end + info.get('cend', 0)
    maxLength = maxEnd - maxStart + 1

    if 'cilen' in info: # count cilen
      if start != end: # not inserted sequence
        if info['cilen'][1] < length: # check min length
          info['cilen'][0] = length
          info['cilen'][1] = length
        elif info['cilen'][0] < length:
          info['cilen'][0] = length
        elif maxLength < info['cilen'][0]: # check max length
          info['cilen'][0] = maxLength
          info['cilen'][1] = maxLength
        elif maxLength < info['cilen'][1]:
          info['cilen'][1] = maxLength

      if info['cilen'][0] == info['cilen'][1]: # remove cilen
        if info['cilen'][0] != 0: # create svlen
          info['svlen'] = info['cilen'][0]

        del info['cilen']
      elif info['cilen'][0] > info['cilen'][1]: # remove cilen
        info['svlen'] = info['cilen'][1]
        del info['cilen']

    if 'max' in info and info.get('cpos', info.get('cilen', None)) is None: # max -> end
      info['end'] = start + length - 1
      del info['max']

    # count imprecise
    info['imprecise'] = 'cpos' in info or \
                        'cend' in info or \
                        'cilen' in info or \
                        'max' in info or \
                        'tracpos' in info or \
                        'tracend' in info

    if not info['imprecise'] and 'svlen' not in info: # add svlen
      info['svlen'] = length
