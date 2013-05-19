#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
import sys
import math

def enum(**enums):
  """
  Enumeration type
  """
  return type('enum', (), enums)

def findMedian(items):
  """
  Return median from items list
  """
  length = len(items)
  half = length / 2.0

  if length % 2: # odd
    return items[int(half)]
  else: # even
    return (items[int(math.floor(half))] + items[int(math.ceil(half))]) / 2.0

def parseRegion(string):
  """
  Parse region of "chr:start-end" format
  """
  region = [None, None, None]

  if string: # parse region
    regionMatch = re.match(r'^([^:]*)(?::([0-9]*)(?:-([0-9]*))?)?$', string)

    if not regionMatch:
      sys.exit("Region has bad format")

    region = list(regionMatch.groups())

    for i in [1, 2]: # start and end of reference
      region[i] = int(region[i]) if region[i] else None

  return region
