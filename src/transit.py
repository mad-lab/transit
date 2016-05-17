#!/usr/bin/env python

# Copyright 2015.
#   Michael A. DeJesus, Chaitra Ambadipudi, and  Thomas R. Ioerger.
#
#
#    This file is part of TRANSIT.
#
#    TRANSIT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#
#    TRANSIT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import absolute_import

import sys
import wx
#Check if wx is the newest 3.0+ version:
try:
    from wx.lib.pubsub import pub
    pub.subscribe
    newWx = True
except AttributeError as e:
    from wx.lib.pubsub import Publisher as pub
    newWx = False

import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.transit_gui as transit_gui
import pytransit.analysis


method_wrap_width = 250
methods = pytransit.analysis.methods


wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"
transit_prefix = "[TRANSIT]"

if __name__ == "__main__":


    #If no arguments, show GUI:
    if len(sys.argv) == 1:
   
        transit_tools.transit_message("Running in GUI Mode")
         
        app = wx.App(False)

        #create an object of CalcFrame
        frame = transit_gui.TnSeekFrame(None)
        #show the frame
        frame.Show(True)
        #start the applications
        app.MainLoop()

    else:
        method_name = sys.argv[1]
        if method_name not in methods:
            print "Error: The '%s' method is unknown." % method_name
            print "Please use one of the known methods (or see documentation to add a new one):"
            for m in methods:
                print "\t - %s" % m
            print "Usage: python %s <method>" % sys.argv[0]
        else:
            methodobj = methods[method_name].method.fromconsole()
            methodobj.Run()            



