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


def run_main():
    (args, kwargs) = cleanargs(sys.argv[1:])
    main(*args, **kwargs)

def main(*args, **kwargs):
    vars = Globals()
    # Check for arguements
    if not args and not kwargs and hasWx:        
        app = wx.App(False)
        form = MyForm(vars)
        form.update_dataset_list()

        form.Show()
        form.Maximize(True)
        app.MainLoop()

        # vars.action not defined, quit...

        if hasattr(vars, 'action'):
            if vars.action=="start":
                verify_inputs(vars)
                if vars.fq2=="": msg = 'running pre-processing on %s' % (vars.fq1)
                else: msg = 'running pre-processing on %s and %s' % (vars.fq1,vars.fq2)
                message(msg)
                message("transposon type: %s" % vars.transposon)
                message("protocol: %s" % vars.protocol)
                save_config(vars)
                driver(vars)

        else:
            pass

    elif not args and not kwargs and not hasWx:
        print "Please install wxPython to run in GUI Mode."
        print "To run in Console Mode please follow these instructions:"
        print ""
        show_help()

    else:
        # Show help if needed
        if "help" in kwargs or "-help" in kwargs:
            show_help()
            sys.exit()

        # Check for strange flags
        known_flags = set(["tn5", "help", "himar1", "protocol", "primer", "reads1",
                           "reads2", "bwa", "ref", "maxreads", "output", "mismatches", "flags",
                           "barseq_catalog_in", "barseq_catalog_out",
                           "window-size", "bwa-alg", "replicon-ids","primer-start-window"])
        unknown_flags = set(kwargs.keys()) - known_flags
        if unknown_flags:
            print "error: unrecognized flags:", ", ".join(unknown_flags)
            show_help()
            sys.exit()

        # Initialize variables
        initialize_globals(vars, args, kwargs)

        # Check inputs make sense
        verify_inputs(vars)
        
        # Print some messages
        if vars.fq2:
            msg = 'running pre-processing on %s' % (vars.fq1)
        else:
            msg = 'running pre-processing on %s and %s' % (vars.fq1, vars.fq2)

        message(msg)
        message("protocol: %s" % vars.protocol)
        message("transposon type: %s" % vars.transposon)

        # Save configuration file
        save_config(vars)
        
        # Run TPP
        driver(vars)


if __name__ == "__main__":
    run_main()

