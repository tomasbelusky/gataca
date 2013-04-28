#!/usr/bin/env python2.5
#-*- coding: utf-8 -*-

__author__ = "Tomáš Beluský"
__date__ = "28.01.2012"

import cgi
import copy
import datetime
from types import *
import traceback

class Logger:
  """
  Log program's flow when exception has occured
  """
  LOG_PATH = "log"

  @staticmethod
  def addItem(key, value):
    """
    Add item into local variable's list
    """
    return "<li><span><a href=\"javascript:void(0)\" onclick=\"hide(this); return false;\">%s</a></span>: <span>%s</span></li>" % (key, value)

  @staticmethod
  def getLocalValue(local):
    """
    Return content of local variable
    """
    try:
      objType = type(local)

      if objType in (ClassType, InstanceType):
        local = local.__dict__
        objType = type(local)

      if objType is UnicodeType: # unicode
        return cgi.escape(local.encode("iso-8859-2"))
      elif objType is StringType: # string
        return cgi.escape(local)
      elif objType in(DictType, DictionaryType): # dictionary
        content = ""

        for key, value in sorted(local.items(), key=lambda l: l[0]):
          key = Logger.getLocalValue(key)
          value = Logger.getLocalValue(value)

          if key is not None and value is not None:
            content += Logger.addItem(key, value)

        return "<ul>%s</ul>" % content
      elif objType in(TupleType, ListType): # tuple or list
        content = ""

        for index, item in enumerate(local):
          value = Logger.getLocalValue(item)

          if value is not None:
            content += Logger.addItem(index, value)

        return "<ul>%s</ul>" % content
      else: # everything else
        return cgi.escape(str(local))
    except Exception:
      pass

    return None

  @staticmethod
  def save(tb):
    """
    Save traceback
    """
    fileContents = []
    excLines = traceback.format_exc().splitlines()
    count = 0

    while tb != None: # go through trace
      count += 1
      locals = copy.copy(tb.tb_frame.f_locals)
      trace = traceback.extract_tb(tb)
      block = "<a name=\"%d\"></a>" % count
      block += "<table class=\"info\">"
      block += "<tr><th>File:</th><td>%s</td></tr>" % trace[0][0]
      block += "<tr><th>Line:</th><td>%d</td></tr>" % trace[0][1]
      block += "<tr><th>Function:</th><td>%s</td></tr>" % trace[0][2]
      block += "<tr><th>Command:</th><td>%s</td></tr>" % trace[0][3]
      block += "<tr><th>Exception:</th><td>%s</td></tr>" % cgi.escape(excLines[-1])
      block += "</table>"
      block += "<table class=\"navigation\"><tr>"
      tb = tb.tb_next
      block += "<td>%s</td>" % (("<a href=\"#%d\">previous</a>" % (count + 1)) if tb != None else "&nbsp;")
      block += "<td>%s</td>" % (("<a href=\"#%d\">next</a>" % (count - 1)) if count != 1 else "&nbsp;")
      block += "</tr></table>"
      block += "<ul>"

      for i in sorted(locals): # add local variables
        value = Logger.getLocalValue(locals[i])

        if value is not None:
          block += Logger.addItem(i, value)

      block += "</ul>"
      fileContents.append(block)

    content = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">'
    content += '<html xmlns="http://www.w3.org/1999/xhtml">'
    content += '<head>'
    content += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'
    content += '<title>gataca log (%s)</title>' % excLines[-1]
    content += '<link rel="stylesheet" type="text/css" href="../interface/style.css" />'
    content += '<script type="text/javascript" src="../interface/script.js"></script>'
    content += '</head>'
    content += '<body>'

    for i in reversed(fileContents): # add trace in reversed order
      content += i

    content += '</body>'
    content += '</html>'

    now = datetime.datetime.now()
    filename = "%s/%s.html" % (Logger.LOG_PATH, now.strftime('%Y%m%d%H%M%S'))
    f = open(filename, 'w')
    f.write(content)
    f.close()
