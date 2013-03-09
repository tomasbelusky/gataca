#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "05.03. 2013"

import sys
import pysam
from datetime import date
from VcfCreator import VcfCreator
from interface import *
from Sample import Sample
from Cluster import Cluster

class Detector:
  cigarOperations = enum(ALIGNMENT=0,
                         INSERTION=1,
                         DELETION=2,
                         SKIPPED=3,
                         SOFTCLIP=4,
                         HARDCLIP=5,
                         PADDING=6,
                         MATCH=7,
                         MISMATCH=8)

  def __init__(self, sample, reference, policy=Sample.policyType.FR, region=None):
    self.__sample = sample
    self.__reference = pysam.Fastafile(reference)
    self.__policy = policy
    self.__region = region
    self.__vcfCreator = VcfCreator(self.__sample, self.__reference)

  def start(self):
    self.__sample.preprocessing()
    return 0

  def write(self, filename):
    self.__vcfCreator.setHeader('fileDate', date.today().strftime('%Y%m%d'))
    self.__vcfCreator.setHeader('source', 'gataca')
    self.__vcfCreator.setInfo('SVTYPE', 1, 'String', 'Type of structural variant')
    self.__vcfCreator.setInfo('SVLEN', 1, 'Integer', 'Difference in length between REF and ALT')
    self.__vcfCreator.setInfo('EVENT', 1, 'String', 'ID of event associated to breakend')
    self.__vcfCreator.setInfo('BNDLEN', 1, 'Integer', 'Length of sequence between with breakends')
    self.__vcfCreator.setInfo('END', 1, 'Integer', 'End position of the variant')
    self.__vcfCreator.setInfo('IMPRECISE', 0, 'Flag', 'Imprecise structural variation')
    self.__vcfCreator.setInfo('CIPOS', 2, 'Integer', 'Confidence interval around POS')
    self.__vcfCreator.setInfo('CIEND', 2, 'Integer', 'Confidence interval around END')
    self.__vcfCreator.setInfo('ISEQ', 1, 'String', 'Imprecise inserted sequence')
    self.__vcfCreator.setInfo('CPREFIX', 2, 'Integer', 'Confidence prefix of inserted sequence')
    self.__vcfCreator.setInfo('CSUFFIX', 2, 'Integer', 'Confidence suffix of inserted sequence')
    self.__vcfCreator.setAlt('DEL', 'Deletion')
    self.__vcfCreator.setAlt('INS', 'Insertion')
    self.__vcfCreator.setAlt('INV', 'Inversion')
    self.__vcfCreator.setAlt('BND', 'Breakend')
    self.__vcfCreator.write(filename)
