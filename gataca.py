#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import optparse
import pysam
import re
import sys
import os

from src.interface.Settings import Settings
from src.resources.Sample import Sample
from src.resources.reads.Read import Read
from src.variations.Detector import Detector

def getParameters(argv):
  """
  Return parsed cmd parameters
  """
  parser = optparse.OptionParser(usage="Usage: %prog [OPTIONS] <sample.bam> <reference.fasta>",
                                 description="Detection of genome variations in mapped reads on reference genome.")

  inout = optparse.OptionGroup(parser, "Input/output")
  inout.add_option("-r", "--region",
                   help="specify region (chr:from-to) of your interest, default: whole genome",
                   metavar="STR")
  inout.add_option("-o", "--output",
                   help="name of output VCF file, default: standard output",
                   metavar="STR")
  parser.add_option_group(inout)

  read = optparse.OptionGroup(parser, "Reads")
  read.add_option("-p", "--policy",
                  help="set how reads were sequenced (fr, rf) [%default]",
                  metavar="STR", choices=('fr', 'rf'), default=Settings.POLICY)
  read.add_option("-q", "--min_quality",
                  help="minimal mapping Phred quality score of read [%default]",
                  type="int", metavar="INT", default=Settings.MIN_QUALITY)
  read.add_option("-l", "--min_length",
                  help="minimal length of split part [%default]",
                  type="int", metavar="INT", default=Settings.MIN_PART_LENGTH)
  parser.add_option_group(read)

  coverage = optparse.OptionGroup(parser, "Depth of coverage")
  coverage.add_option("-w", "--window_size",
                      help="size of window for getting coverage [%default]",
                      type="int", metavar="INT", default=Settings.WINDOW_SIZE)
  coverage.add_option("-c", "--coverage",
                      help="interval (min,max) of accepted coverage in windows, default: estimate from reads",
                      metavar="STR", default="%d,%d" % (Settings.MIN_COVERAGE, Settings.MAX_COVERAGE))
  coverage.add_option("-a", "--coverage_core",
                      help="core of windows from which minimal and maximal allowed coverage will be estimated [%default]",
                      type="float", metavar="FLOAT", default=Settings.INSERT_CORE)
  coverage.add_option("-u", "--min_coverage_count",
                      help="minimal number of windows in core [%default]",
                      type="int", metavar="INT", default=Settings.MIN_COVERAGE_COUNT)
  parser.add_option_group(coverage)

  insert = optparse.OptionGroup(parser, "Insert size")
  insert.add_option("-i", "--insert_size",
                    help="interval (min,max) of accepted size between reads, default: estimate from reads",
                    metavar="STR", default="%d,%d" % (Settings.MIN_INSERT, Settings.MAX_INSERT))
  insert.add_option("-n", "--insert_reads",
                    help="number of reads from which insert size will be estimated [%default]",
                    type="int", metavar="INT", default=Settings.INSERT_READS)
  insert.add_option("-e", "--insert_core",
                    help="core of reads from which minimal and maximal allowed insert size will be estimated [%default]",
                    type="float", metavar="FLOAT", default=Settings.INSERT_CORE)
  insert.add_option("-m", "--min_insert_count",
                    help="minimal number of reads in core [%default]",
                    type="int", metavar="INT", default=Settings.MIN_INSERT_COUNT)
  parser.add_option_group(insert)

  var = optparse.OptionGroup(parser, "Variations")
  var.add_option("-v", "--min_confidence",
                 help="minimal confidence about variation [%default]",
                 type="float", metavar="FLOAT", default=Settings.MIN_CONFIDENCE)
  parser.add_option_group(var)

  parser.set_defaults(region=None)
  return parser.parse_args(argv[1:])

def checkPositive(name, value):
  """
  Check if value is positive
  """
  if value <= 0:
    raise Exception("%s must be positive integer value" % name)

  return value

def checkInterval(name, interval, value, minValue, minInclude, maxValue, maxInclude):
  """
  Check if value is in interval
  """
  if (minInclude and value < minValue) or (not minInclude and value <= minValue) or \
     (maxInclude and maxValue < value) or (not maxInclude and maxValue <= value):
    raise Exception("%s must be in interval %s" % (name, interval))

  return value

def parseInterval(value, interval):
  """
  Parse interval in MIN,MAX format
  """
  match = re.match(r'(?P<min>\d+),(?P<max>\d+)', interval)

  if not match:
    raise Exception("%s interval has bad format, should be MIN,MAX" % value)

  minValue = int(match.group('min'))
  maxValue = int(match.group('max'))

  if max < min:
    Exception("%s interval cannot have minimal value bigger then maximal value" % value)

  return minValue, maxValue

def main(argv):
  """
  Main function
  """
  (options, args) = getParameters(argv)
  params = options.__dict__

  if len(args) < 2: # check required input files
    raise Exception("Please specify file with sample and file with reference genome")

  if not params['output']: # set output
    params['output'] = sys.stdout

  # parse region
  region = [None, None, None]

  if params['region']: # parse region
    regionMatch = re.match(r'^([^:]*)(?::([0-9]*)(?:-([0-9]*))?)?$', params['region'])

    if not regionMatch:
      raise Exception("Region has bad format")

    region = list(regionMatch.groups())

    for i in [1, 2]: # start and end of reference
      region[i] = int(region[i]) if region[i] else None

  # create tmp dir
  if not os.path.exists(Sample.TMP_PATH):
    os.makedirs(Sample.TMP_PATH)

  # set settings
  Settings.REFERENCE = region[0]
  Settings.START = region[1]
  Settings.END = region[2]
  Settings.POLICY = Read.ptype.FR if params['policy'] == "fr" else Read.ptype.RF
  Settings.MIN_QUALITY = params['min_quality']
  Settings.MIN_PART_LENGTH = checkPositive("Minimal length", params['min_length'])
  Settings.WINDOW_SIZE = checkPositive("Window size", params['window_size'])
  (Settings.MIN_COVERAGE, Settings.MAX_COVERAGE) = parseInterval("Coverage", params['coverage'])
  Settings.COVERAGE_CORE = checkInterval("Coverage core", "(0,1]", params['coverage_core'], 0, False, 1, True)
  Settings.MIN_COVERAGE_COUNT = checkPositive("Minimal coverage count", params['min_coverage_count'])
  (Settings.MIN_INSERT, Settings.MAX_INSERT) = parseInterval("Insert size", params['insert_size'])
  Settings.INSERT_READS = checkPositive("Insert reads", params['insert_reads'])
  Settings.INSERT_CORE = checkInterval("Insert core", "(0,1]", params['insert_core'], 0, False, 1, True)
  Settings.MIN_INSERT_COUNT = checkPositive("Minimal insert count", params['min_insert_count'])
  Settings.MIN_CONFIDENCE = checkInterval("Minimal confidence", "[0,1]", params['min_confidence'], 0, True, 1, True)

  # create objects and start
  refgenome = pysam.Fastafile(args[1])
  sample = Sample(args[0], refgenome)

  if Settings.REFERENCE and Settings.REFERENCE not in sample.getReferences(): # check reference name
    raise Exception("Unknown chromosome")

  detector = Detector(sample, refgenome, params['output'])
  detector.start()
  sample.close()

if __name__ == "__main__":
  #try:
    main(sys.argv)
  #except Exception, ex:
  #  sys.stderr.write("%s: error: %s\n" % (sys.argv[0], ex))
  #  sys.exit(1)

  #sys.exit(0)
