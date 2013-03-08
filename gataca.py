#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import os
import sys
import re
import optparse
import pysam
from src.Detector import Detector

def getParameters(argv):
  parser = optparse.OptionParser(description="Detection of genome variations from mapped reads on reference genome.")
  parser.add_option("-r", "--region", help="specify region of your interest, if not set variations will be find in whole genome", default=None, metavar="CHR:FROM-TO")
  parser.add_option("-p", "--policy", help="set how reads were sequenced (ff, fr, rf), default=%default", default="ff", choices=('ff', 'fr', 'rf'))
  parser.add_option("-s", "--sample", help="BAM file with aligned reads to reference genome [-f]")
  parser.add_option("-f", "--reference", help="FASTA file with reference genome")
  parser.set_defaults(region=None)
  return parser.parse_args(argv[1:])

def main(argv):
  (options, args) = getParameters(argv)
  params = options.__dict__

  if not params['reference'] or not params['sample']:
    return "Please specify file with sample and file with reference genome"

  try:  
    detector = Detector(params['sample'], params['reference'], params['policy'], params['region'])
    return detector.start()
  except Exception, ex:
    sys.stderr.write("%s: error: %s\n" % (argv[0], ex))
    return 1

  """
  # open BAM and get the iterator
  samfile = pysam.Samfile(args[0], "rb")

  if options.region != None:
    it = samfile.fetch(region=options.region)
  else:
    it = samfile.fetch()

  # calculate the end position and print out BED
  take = (0, 1, 4, 7, 8)
  outfile = sys.stdout

  for read in it:
    if read.is_unmapped:
      continue

    # compute total length on reference
    t = sum([l for op, l in read.cigar if op in take])
    strand = "-" if read.is_reverse else "+"
    outfile.write("%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n" % (samfile.getrname(read.rname), read.pos, read.pos+t, read.tlen, read.qname, read.mapq, read.flag, read.seq, read.cigar, read.qual))
  """

if __name__ == "__main__":
  sys.exit(main(sys.argv))
