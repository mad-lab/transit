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

import sys
import glob
import os
import time
import math
import sys
import re
import shutil
import platform
import gzip

from tpp_tools import *
from tpp_gui import *


def main():
    # if -nowin is command-line arg, skip the GUI and set filenames in vars

    vars = Globals()
    #vars.version = "$Revision: 1.5 $".split()[1]

    if len(sys.argv) <= 1 and hasWx:
        app = wx.App(False)
        form = MyForm(vars)
        form.update_dataset_list()

        form.Show()
        app.MainLoop()

        # vars.action not defined, quit...

        if hasattr(vars, 'action'):
            if vars.action=="start":
                verify_inputs(vars)
                if vars.fq2=="": msg = 'running pre-processing on %s' % (vars.fq1)
                else: msg = 'running pre-processing on %s and %s' % (vars.fq1,vars.fq2)
                message(msg)
                message("transposon type: %s" % vars.transposon)
                save_config(vars)
                driver(vars)

        else:
            pass

    elif len(sys.argv) <= 1 and not hasWx:
        print "Please install wxPython to run in GUI Mode."
        print "To run in Console Mode please follow these instructions:"
        print ""
        show_help()

    else:
        flag = False
        initialize_globals(vars)
        for i in range(0, len(sys.argv)):
            if sys.argv[i] == '-help':
                show_help()
                sys.exit()
            if sys.argv[i] == '-tn5':
                vars.transposon = 'Tn5'
            if sys.argv[i] == '-reads1':
                vars.fq1 = sys.argv[i+1]
            elif sys.argv[i] == '-reads2':
                flag = True
                vars.fq2 = sys.argv[i+1]
            elif sys.argv[i] == '-bwa':
                vars.bwa = sys.argv[i+1]
            elif sys.argv[i] == '-ref':
                vars.ref = sys.argv[i+1]
            elif sys.argv[i] == '-maxreads':
                vars.maxreads = int(sys.argv[i+1])
            elif sys.argv[i] == '-prefix':
                vars.base = sys.argv[i+1]
            elif sys.argv[i] == '-mismatches':
                vars.mm1 = int(sys.argv[i+1])
        if flag==False: vars.fq2 = ""
        if vars.fq2=="": msg = 'running pre-processing on %s' % (vars.fq1)
        else: msg = 'running pre-processing on %s and %s' % (vars.fq1,vars.fq2)
        message(msg)
        message("transposon type: %s" % vars.transposon)
        verify_inputs(vars)
        save_config(vars)
        driver(vars)

if __name__ == "__main__":
    main()

