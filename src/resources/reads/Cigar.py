#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "12.04. 2013"

from src.interface.interface import *

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
  # sum of these operation's lengths give seqence length of read
  sums = (op.ALIGNMENT, op.INSERTION, op.SOFTCLIP, op.MATCH, op.MISMATCH)
  abbr = 'MIDNSHP=X' # abbreviations
