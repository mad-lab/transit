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
import os

try:
    import wx
    hasWx = True
    #Check if wx is the newest 3.0+ version:
    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False
except Exception as e:
    hasWx = False
    newWx = False

import math
import ntpath
import numpy
import scipy.optimize
import scipy.stats
import pytransit

import pytransit.tnseq_tools as tnseq_tools

def aton(aa):
    #TODO: Write docstring
    return(((aa-1)*3)+1)

def parseCoords(strand, aa_start, aa_end, start, end):
    #TODO: Write docstring
    if strand == "+":
        return((aton(aa_start) + start,  aton(aa_end) + start))
    # Coordinates are Reversed... to match with Trash FILE TA coordinates
    if strand == "-":
        return((end - aton(aa_end), end - aton(aa_start)))



def fetch_name(filepath):
    #TODO: Write docstring
    return os.path.splitext(ntpath.basename(filepath))[0]


def basename(filepath):
    #TODO: Write docstring
    return ntpath.basename(filepath)


def cleanargs(rawargs):
    #TODO: Write docstring
    args = []
    kwargs = {}
    count = 0
    while count < len(rawargs):
        if rawargs[count].startswith("-"):
            try:
                kwargs[rawargs[count][1:]] = rawargs[count+1]
                count += 1
            except IndexError as IE:
                kwargs[rawargs[count][1:]] = True
        else:
            args.append(rawargs[count])
        count += 1

    return (args, kwargs)


def getTabTableData(path, colnames):
    #TODO: Write docstring
    row = 0
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        tmp[-1] = tmp[-1].strip()
        rowdict = dict([(colnames[i], tmp[i]) for i in range(len(colnames))])
        data.append((row, rowdict))
        row+=1

    return data


def ShowMessage(MSG=""):
    #TODO: Write docstring
    wx.MessageBox(MSG, 'Info',
        wx.OK | wx.ICON_INFORMATION)

def ShowAskWarning(MSG=""):
    #TODO: Write docstring
    dial = wx.MessageDialog(None, MSG, 'Warning',
        wx.OK | wx.CANCEL | wx.ICON_EXCLAMATION)
    return dial.ShowModal()

def ShowError(MSG=""):
    #TODO: Write docstring
    dial = wx.MessageDialog(None, MSG, 'Error',
        wx.OK | wx.ICON_ERROR)
    dial.ShowModal()

def transit_message(msg=""):
    #TODO: Write docstring
    print pytransit.prefix, msg


def transit_error(text):
    #TODO: Write docstring
    transit_message(text)
    try:
        ShowError(text)
    except:
        pass


def validate_annotation(annotation):
    #TODO: Write docstring
    if not annotation:
        transit_error("Error: No annotation file selected!")
        return False
    return True

def validate_control_datasets(ctrldata):
    #TODO: Write docstring
    if len(ctrldata) == 0:
        transit_error("Error: No control datasets selected!")
        return False
    return True

def validate_both_datasets(ctrldata, expdata):
    #TODO: Write docstring
    if len(ctrldata) == 0 and len(expdata) == 0:
        transit_error("Error: No datasets selected!")
        return False
    elif len(ctrldata) == 0:
        transit_error("Error: No control datasets selected!")
        return False
    elif len(expdata) == 0:
        transit_error("Error: No experimental datasets selected!")
        return False
    else:
        return True
    

def validate_filetypes(datasets, transposons, justWarn=True):
    #TODO: Write docstring
    unknown = tnseq_tools.get_unknown_file_types(datasets, transposons)
    if unknown:
        if justWarn:
            answer = ShowAskWarning("Warning: Some of the selected datasets look like they were created using transposons that this method was not intended to work with: %s. Proceeding may lead to errors. Click OK to continue." % (",". join(unknown)))
            if answer == wx.ID_CANCEL:
                return False
            else:
                return True
        else:
            transit_error("Error: Some of the selected datasets look like they were created using transposons that this method was not intended to work with: %s." % (",". join(unknown)))
            return False
    return True





