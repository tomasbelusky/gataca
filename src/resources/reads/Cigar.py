#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "12.04. 2013"

from interface.interface import *

class Cigar:
  """
  Represents possible operations in reads
  """
  op = enum(# operations
             ALIGNMENT=0,
             INSERTION=1,
             DELETION=2,
             SKIPPED=3,
             SOFTCLIP=4,
             HARDCLIP=5,
             PADDING=6,
             MATCH=7,
             MISMATCH=8)
  sums = (op.ALIGNMENT, op.INSERTION, op.MATCH, op.MISMATCH)
  abbr = 'MIDNSHP=X' # abbreviations
