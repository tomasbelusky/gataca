#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "11.04. 2013"

import os
import subprocess

class Bwa:
  """
  Represents bwa tool
  """

  def index(self, filename):
    """
    Create index file for reference file
    """
    with open(os.devnull, 'wb') as devnull:
      subprocess.call(["bwa",
                       "index",
                       "%s.fasta" % filename],
                       stdout=devnull,
                       stderr=devnull)

  def align(self, filename):
    """
    Align reads in filename to reference genome stored in refname
    """
    with open(os.devnull, 'wb') as devnull:
      with open("%s.sai" % filename, "w") as sai:
        subprocess.call(["bwa",
                         "aln",
                         "%s.fasta" % filename,
                         "%s.fastq" % filename],
                         stdout=sai,
                         stderr=devnull)
      
      with open("%s.sam" % filename, "w") as sam:
        subprocess.call(["bwa",
                         "samse",
                         "%s.fasta" % filename,
                         "%s.sai" % filename,
                         "%s.fastq" % filename],
                         stdout=sam,
                         stderr=devnull)
