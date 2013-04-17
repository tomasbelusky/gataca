#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import re
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

def parseInterval(value, interval):
  """
  Parse interval in MIN,MAX format
  """
  match = re.match(r'(?P<min>\d+),(?P<max>\d+)', interval)

  if not match:
    raise Exception("%s interval has bad format, should be MIN,MAX" % value)

  return int(match.group('min')), int(match.group('max'))
