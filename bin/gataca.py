#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import optparse
import pysam
import re
import sys

from resources.Sample import Sample
from variations.Detector import Detector

def getParameters(argv):
  """
  Return parsed cmd parameters
  """
  parser = optparse.OptionParser(description="Detection of genome variations from mapped reads on reference genome.")
  parser.add_option("-r", "--region", help="specify region of your interest, if not set variations will be find in whole genome", default=None, metavar="CHR:FROM-TO")
  parser.add_option("-p", "--policy", help="set how reads were sequenced (ff, fr, rf), default=%default", default="fr", choices=('fr', 'rf'))
  parser.add_option("-s", "--sample", help="file with aligned reads to reference genome [-f] in BAM format")
  parser.add_option("-o", "--output", help="Name of output VCF file, by default output will be print on standard output")
  parser.add_option("-f", "--reference", help="FASTA file with reference genome")
  parser.add_option("-i", "--insert_size", help="Interval of accepted size between reads", metavar="MIN-MAX", default="0-0")
  parser.add_option("-d", "--dont_use_coverage", help="Don't use coverage in variation detection", action="store_false", default=True)
  parser.set_defaults(region=None)
  return parser.parse_args(argv[1:])

def main(argv):
  """
  Main function
  """
  (options, args) = getParameters(argv)
  params = options.__dict__

  if not params['reference'] or not params['sample']: # check required params
    raise Exception("Please specify file with sample and file with reference genome")

  if not params['output']: # set output
    params['output'] = sys.stdout

  if params['policy'] == "fr": # set policy
    policy = Sample.ptype.FR
  else:
    policy = Sample.ptype.RF

  # set region
  region = [None, None, None]

  if params['region']: # parse region
    regionMatch = re.match(r'^([^:]*)(?::([0-9]*)(?:-([0-9]*))?)?$', params['region'])

    if not regionMatch:
      raise Exception("Region has bad format")

    region = list(regionMatch.groups())

    for i in [1, 2]: # start and end of reference
      region[i] = int(region[i]) if region[i] else None

  # set min and max insert size
  insertSizeMatch = re.match(r'(?P<min>\d+)-(?P<max>\d+)', params['insert_size'])

  if not insertSizeMatch:
    raise Exception("Insert size interval has bad format")

  # create objects, start and write output
  refgenome = pysam.Fastafile(params['reference'])
  sample = Sample(params['sample'], refgenome, policy, int(insertSizeMatch.group('min')), int(insertSizeMatch.group('max')), params['dont_use_coverage'])

  if region[0] and region[0] not in sample.getReferences(): # check reference name
    raise Exception("Unknown reference")

  detector = Detector(sample, refgenome, policy, region[0], region[1], region[2])

  if not detector.start():
    detector.write(params['output'])

if __name__ == "__main__":
  #try:
    main(sys.argv)
  #except Exception, ex:
  #  sys.stderr.write("%s: error: %s\n" % (sys.argv[0], ex))
  #  sys.exit(1)

  #sys.exit(0)
