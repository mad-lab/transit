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

import pytransit.__main__




if __name__ == "__main__":
    pytransit.__main__.main()



