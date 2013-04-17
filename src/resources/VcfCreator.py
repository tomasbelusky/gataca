#!/usr/bin/python2.7
# -*- encoding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "09.03. 2013"

import types

class VcfCreator:
  """
  Creator of VCF output
  """

  def __init__(self, filename):
    """
    Initialize variables
    """
    self.__filename = filename
    self.__headers = []
    self.__contigs = []
    self.__infos = []
    self.__alts = []
    self.__records = []

  def addHeader(self, key, value):
    """
    Add header
    """
    self.__headers.append((key, value))

  def __addContigAttribute(self, attributes, key, attribute, result, addApostrophe=False):
    """
    Add attribute with possible apostrophes into result string
    """
    if key in attributes: # attribute exists
      if result: # there is previous attribute
        result += ","

      if addApostrophe:
        result += "%s=\"%s\"" % (attribute, attributes[key])
      else:
        result += "%s=%s" % (attribute, attributes[key])

    return result

  def addContig(self, contig):
    """
    Add reference contig
    """
    result = self.__addContigAttribute(contig, 'SN', 'ID', "")
    result = self.__addContigAttribute(contig, 'LN', 'length', result)
    result = self.__addContigAttribute(contig, 'AS', 'assembly', result)
    result = self.__addContigAttribute(contig, 'M5', 'md5', result)
    result = self.__addContigAttribute(contig, 'SP', 'species', result, True)
    result = self.__addContigAttribute(contig, 'UR', 'URL', result)
    self.__contigs.append(result)

  def addInfo(self, iid, number, itype, description):
    """
    Add info
    """
    self.__infos.append((iid, number, itype, description))

  def addAlt(self, aid, description):
    """
    Add alternate
    """
    self.__alts.append((aid, description))

  def addRecord(self, record):
    """
    Add record with variation
    """
    if len(record.strip()):
      self.__records.append(record)

  def write(self, name):
    """
    Write VCF into stdout or file
    """
    if type(name) == types.FileType:
      output = name
    else:
      output = open(name, 'w')

    """
    INFO: don't print in test phase
    output.write("##fileformat=VCFv4.1\n")

    for key, value in self.__headers: # write header
      output.write("##%s=%s\n" % (key, value))

    output.write("##reference=%s\n" % self.__filename)

    for contig in self.__contigs: # write contigs
      output.write("##contig=<%s>\n" % contig)

    for iid, number, itype, description in self.__infos: # write info
      output.write("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n" % (iid, number, itype, description))

    for aid, description in self.__alts: # write alt
      output.write("##ALT=<ID=%s,Description=\"%s\">\n" % (aid, description))

    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    """

    for rec in self.__records: # write records
      output.write("%s\n" % rec)

    if type(name) != types.FileType:
      output.close()
