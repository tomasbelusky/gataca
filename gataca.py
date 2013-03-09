#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import os
import sys
import re
import optparse
import pysam
from src.Sample import Sample
from src.Detector import Detector

def getParameters(argv):
  parser = optparse.OptionParser(description="Detection of genome variations from mapped reads on reference genome.")
  parser.add_option("-r", "--region", help="specify region of your interest, if not set variations will be find in whole genome", default=None, metavar="CHR:FROM-TO")
  parser.add_option("-p", "--policy", help="set how reads were sequenced (ff, fr, rf), default=%default", default="fr", choices=('ff', 'fr', 'rf'))
  parser.add_option("-s", "--sample", help="file with aligned reads to reference genome [-f] in BAM format")
  parser.add_option("-o", "--output", help="Name of output VCF file, by default output will be print on standard output")
  parser.add_option("-f", "--reference", help="FASTA file with reference genome")
  parser.set_defaults(region=None)
  return parser.parse_args(argv[1:])

def main(argv):
  (options, args) = getParameters(argv)
  params = options.__dict__

  if not params['reference'] or not params['sample']:
    return "Please specify file with sample and file with reference genome"

  if params['policy'] == "ff":
    policy = Sample.policyType.FF
  elif params['policy'] == "fr":
    policy = Sample.policyType.FR
  else:
    policy = Sample.policyType.RF

  #try:
  sample = Sample(params['sample'], params['reference'], policy)
  detector = Detector(sample, params['reference'], policy, params['region'])
  
  if detector.start():
    return

  if not params['output']:
    params['output'] = sys.stdout

  detector.write(params['output'])
  #except Exception, ex:
  #  sys.stderr.write("%s: error: %s\n" % (argv[0], ex))
  #  return 1

if __name__ == "__main__":
  sys.exit(main(sys.argv))
