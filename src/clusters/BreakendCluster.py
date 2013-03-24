#!/usr/bin/python2.7
#-*- encoding: utf-8 -*-

_author_ = "Tomáš Beluský"
_date_ = "13.03. 2013"

from AbstractCluster import AbstractCluster
from src.Variation import Variation

class BreakendCluster:
  """
  Represents cluster with variations which are represents with breakends in VCF format
  """
  def _init_(self, reference, sample):
    """
    Initialize variables
    """
    AbstractCluster._init_(self, reference, sample)

  def _process(self):
    """
    Join variations together for printing them in VCF output
    TODO: whole method
    """
    pass

  def compare(self, variation):
    """
    Compare if variation fits into actual cluster
    TODO: whole method
    """
    return False

  def __str__(self):
    """
    Print cluster in VCF format
    TODO: whole method
    """
    return ""
