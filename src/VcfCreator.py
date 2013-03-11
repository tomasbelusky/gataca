#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "09.03. 2013"

import re
import types
import codecs

class VcfCreator:
  def __init__(self, sample, reference):
    self.__sample = sample
    self.__reference = reference
    self.__headers = []
    self.__infos = []
    self.__alts = []
    self.__records = []

  def setHeader(self, key, value):
    self.__headers.append((key, value))

  def setInfo(self, iid, number, itype, description):
    self.__infos.append((iid, number, itype, description))

  def setAlt(self, aid, description):
    self.__alts.append((aid, description))

  def setRecord(self, record):
    self.__records.append(record)

  def __addAttribute(self, attributes, key, attribute, addComma, addApostrophe=False):
    result = ""

    if key in attributes:
      if addComma:
        result += ","

      addComma = True

      if addApostrophe:
        result += "%s=\"%s\"" % (attribute, attributes[key])
      else:
        result += "%s=%s" % (attribute, attributes[key])

    return result, addComma

  def write(self, name):
    if type(name) != types.FileType:
      output = open(name, 'w')
    else:
      output = name

    output.write("##fileformat=VCFv4.1\n")

    for key, value in self.__headers: # write header
      output.write("##%s=%s\n" % (key, value))

    output.write("##reference=%s\n" % self.__reference.filename)

    for contig in self.__sample.getRefSequences(): # write contigs
      (result1, addComma) = self.__addAttribute(contig, 'SN', 'ID', False)
      (result2, addComma) = self.__addAttribute(contig, 'LN', 'length', addComma)
      (result3, addComma) = self.__addAttribute(contig, 'AS', 'assembly', addComma)
      (result4, addComma) = self.__addAttribute(contig, 'M5', 'md5', addComma)
      (result5, addComma) = self.__addAttribute(contig, 'SP', 'species', addComma, True)
      (result6, addComma) = self.__addAttribute(contig, 'UR', 'URL', addComma)
      output.write("##contig=<%s%s%s%s%s%s>\n" % (result1, result2, result3, result4, result5, result6))

    for iid, number, itype, description in self.__infos: # write info
      output.write("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n" % (iid, number, itype, description))

    for aid, description in self.__alts: # write alt
      output.write("##ALT=<ID=%s,Description=\"%s\">\n" % (aid, description))

    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for rec in self.__records: # write records
      output.write("%s\n" % rec)

    if type(name) != types.FileType:
      output.close()
