#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "17.04. 2013"

class Settings:
  """
  Command line settings
  """
  # region
  REFERENCE = None # reference of interest
  START = None # start position of specified region
  END = None # end position of specified region

  # reads
  POLICY = "fr" # how reads where sequenced
  MIN_QUALITY = 30 # minimal quality to not be filtered out
  MIN_PART_LENGTH = 10 # minimal length of split parts

  # coverage
  WINDOW_SIZE = 100 # length of window for getting coverage
  MIN_COVERAGE = 0 # minimal coverage
  MAX_COVERAGE = 0 # maximal coverage

  # insert size
  MIN_INSERT_SIZE = 0 # minimal insert size
  MAX_INSERT_SIZE = 0 # maximal insert size
  READS_NUM = 10000 # number of reads to estimate insert size
  CORE = 0.1 # core from <0,1> interval for getting insert size
  MIN_CORE_COUNT = 10 # minimal count of iterms in core
