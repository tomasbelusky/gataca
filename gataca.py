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

  inoutGroup = optparse.OptionGroup(parser, "Input/output")
  inoutGroup.add_option("-r",
                        "--region",
                        help="specify region (chr:from-to) of your interest, default: whole genome",
                        metavar="STR")
  inoutGroup.add_option("-o",
                        "--output",
                        help="name of output VCF file, default: standard output",
                        metavar="STR")
  parser.add_option_group(inoutGroup)

  readGroup = optparse.OptionGroup(parser, "Reads")
  readGroup.add_option("-p",
                        "--policy",
                        help="set how reads were sequenced (fr, rf) [%default]",
                        metavar="STR",
                        choices=('fr', 'rf'),
                        default=Settings.POLICY)
  readGroup.add_option("-q",
                       "--min_quality",
                       help="minimal mapping Phred quality score of read [%default]",
                       type="int",
                       metavar="INT",
                       default=Settings.MIN_QUALITY)
  readGroup.add_option("-l",
                       "--min_length",
                       help="minimal length of split part [%default]",
                       type="int",
                       metavar="INT",
                       default=Settings.MIN_PART_LENGTH)
  parser.add_option_group(readGroup)

  coverageGroup = optparse.OptionGroup(parser, "Depth of coverage")
  coverageGroup.add_option("-d",
                           "--dont_use_coverage",
                           help="don't use coverage in variation detection",
                           action="store_false",
                           default=Settings.COUNT_COVERAGE)
  coverageGroup.add_option("-w",
                           "--window_size",
                           help="size of window for getting coverage [%default]",
                           type="int",
                           metavar="INT",
                           default=Settings.WINDOW_SIZE)
  parser.add_option_group(coverageGroup)

  insertSizeGroup = optparse.OptionGroup(parser, "Insert size")
  insertSizeGroup.add_option("-i",
                             "--interval",
                             help="interval (min-max) of accepted size between reads, default: estimate from reads",
                             metavar="STR",
                             default="%d-%d" % (Settings.MIN_INSERT_SIZE, Settings.MAX_INSERT_SIZE))
  insertSizeGroup.add_option("-n",
                             "--reads_num",
                             help="number of reads from which insert size will be estimated [%default]",
                             type="int",
                             metavar="INT",
                             default=Settings.READS_NUM)
  insertSizeGroup.add_option("-c",
                             "--core",
                             help="core of reads_num from which minimal and maximal insert size will be estimated [%default]",
                             type="float",
                             metavar="FLOAT",
                             default=Settings.CORE)
  insertSizeGroup.add_option("-m",
                             "--min_core_count",
                             help="minimal number of reads in core [%default]",
                             type="int",
                             metavar="INT",
                             default=Settings.MIN_CORE_COUNT)
  parser.add_option_group(insertSizeGroup)

  parser.set_defaults(region=None)
  return parser.parse_args(argv[1:])

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

  # parse min and max insert size
  insertSizeMatch = re.match(r'(?P<min>\d+)-(?P<max>\d+)', params['interval'])

  if not insertSizeMatch:
    raise Exception("Insert size interval has bad format")

  # create tmp dir
  if not os.path.exists(Sample.TMP_PATH):
    os.makedirs(Sample.TMP_PATH)

  # set settings
  Settings.REFERENCE = region[0]
  Settings.START = region[1]
  Settings.END = region[2]
  Settings.POLICY = Read.ptype.FR if params['policy'] == "fr" else Read.ptype.RF
  Settings.MIN_QUALITY = params['min_quality']
  Settings.MIN_PART_LENGTH = params['min_length']
  Settings.COUNT_COVERAGE = params['dont_use_coverage']
  Settings.WINDOW_SIZE = params['window_size']
  Settings.MIN_INSERT_SIZE = int(insertSizeMatch.group('min'))
  Settings.MAX_INSERT_SIZE = int(insertSizeMatch.group('max'))
  Settings.READS_NUM = params['reads_num']
  Settings.CORE = params['core']
  Settings.MIN_CORE_COUNT = params['min_core_count']

  # create objects, start and write output
  refgenome = pysam.Fastafile(args[1])
  sample = Sample(args[0], refgenome)

  if Settings.REFERENCE and Settings.REFERENCE not in sample.getReferences(): # check reference name
    raise Exception("Unknown reference")

  detector = Detector(sample, refgenome)

  if not detector.start():
    detector.write(params['output'])

if __name__ == "__main__":
  #try:
    main(sys.argv)
  #except Exception, ex:
  #  sys.stderr.write("%s: error: %s\n" % (sys.argv[0], ex))
  #  sys.exit(1)

  #sys.exit(0)
