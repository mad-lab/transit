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


#importing wx files

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

import os
import time
import datetime
import threading
import numpy
import matplotlib.pyplot as plt
import multiprocessing as mp
import math

from functools import partial

import traceback

# trash view stuff
import transit
import transit.analysis

import transit.trash as trash
import transit.transit_gui as transit_gui
import transit.transit_tools as transit_tools
import transit.fileDisplay as fileDisplay
import transit.qcDisplay as qcDisplay
import transit.images as images



method_wrap_width = 250
methods = transit.analysis.methods


wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"
transit_prefix = "[TRANSIT]"

#inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TnSeekFrame(transit_gui.MainFrame):
    #constructor
    def __init__(self,parent):
        #initialize parent class
        transit_gui.MainFrame.__init__(self,parent)


        self.workdir = os.getcwd()
        self.annotation = ""
        self.transposons = ["himar1", "tn5"]
        #import pkgutil
        #print [x for x in pkgutil.iter_modules(['transit/analysis'])]
        #print gumbel.Gumbel.__bases__

        self.logoImg.SetBitmap(images.TRANSIT.GetImage().ConvertToBitmap())

        self.index_ctrl = 0
        self.list_ctrl.InsertColumn(0, 'File', width=210)
        self.list_ctrl.InsertColumn(1, 'Total Reads', width=85)
        self.list_ctrl.InsertColumn(2, 'Density', width=85)
        self.list_ctrl.InsertColumn(3, 'Mean Count', width=90)
        self.list_ctrl.InsertColumn(4, 'Max Count', width=85)
        self.list_ctrl.InsertColumn(5, 'Full Path', width=403)


        self.index_exp = 0
        self.list_exp.InsertColumn(0, 'File', width=210)
        self.list_exp.InsertColumn(1, 'Total Reads', width=85)
        self.list_exp.InsertColumn(2, 'Density', width=85)
        self.list_exp.InsertColumn(3, 'Mean Count', width=90)
        self.list_exp.InsertColumn(4, 'Max Count', width=85)
        self.list_exp.InsertColumn(5, 'Full Path',width=403)


        self.index_file = 0
        self.list_files.InsertColumn(0, 'Name', width=350)
        self.list_files.InsertColumn(1, 'Type', width=100)
        self.list_files.InsertColumn(2, 'Date', width=230)
        self.list_files.InsertColumn(3, 'Full Path',width=403)


        self.verbose = True


        self.statusBar.SetStatusText("Welcome to TRANSIT")
        self.progress_count = 0
        pub.subscribe(self.setProgressRange, "progressrange")
        pub.subscribe(self.updateProgress, "progress")
        pub.subscribe(self.updateStatus, "status")
        pub.subscribe(self.addFile, "file")
        pub.subscribe(self.finishRun, "finish")
        pub.subscribe(self.saveHistogram, "histogram")   
 
        #self.outputDirPicker.SetPath(os.path.dirname(os.path.realpath(__file__)))

        methodChoiceChoices = [ "[Choose Method]"]
        for name in methods:
            methods[name].gui.definePanel(self)
            #methods[name].gui.panel.BackgroundColour = (0, 200, 20)
            self.methodSizer.Add(methods[name].gui.panel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
            methods[name].gui.Hide()

            if "himar1" in methods[name].transposons:
                tempMenuItem = wx.MenuItem( self.himar1MenuItem, wx.ID_ANY, methods[name].fullname(), wx.EmptyString, wx.ITEM_NORMAL )
                self.Bind( wx.EVT_MENU, partial(self.MethodSelectFunc,  methods[name].fullname()), tempMenuItem )
                self.himar1MenuItem.AppendItem( tempMenuItem )
            if "tn5" in methods[name].transposons:
                tempMenuItem = wx.MenuItem( self.tn5MenuItem, wx.ID_ANY, methods[name].fullname(), wx.EmptyString, wx.ITEM_NORMAL )
                self.Bind( wx.EVT_MENU, partial(self.MethodSelectFunc, methods[name].fullname()), tempMenuItem )
                self.tn5MenuItem.AppendItem( tempMenuItem )

        #progress
        self.progressPanel = wx.Panel( self.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        progressSizer = wx.BoxSizer( wx.VERTICAL )

        self.progressLabel = wx.StaticText( self.progressPanel, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.progressLabel.Wrap( -1 )
        progressSizer.Add( self.progressLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.progress = wx.Gauge( self.progressPanel, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
        progressSizer.Add( self.progress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        #self.progressPanel.BackgroundColour = (0, 0, 250)
        self.progressPanel.SetSizer( progressSizer )
        self.progressPanel.SetMaxSize(wx.Size(100, 100))
        self.progressPanel.Layout()

        #progressSizer.Fit( self.progressPanel )
        self.methodSizer.Add( self.progressPanel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        self.globalLabel = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"Florf", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.methodSizer.Add( self.globalLabel, 1, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        #self.methodSizer.Hide()
        self.progress.SetRange(50)
        #########
        
        self.optionsWindow.Fit()

        self.HideProgressSection()
        self.HideGlobalOptions()

    


    def Exit(self, event):
        """Exit Menu Item"""
        if self.verbose:
            print transit_prefix, "Exiting Transit"
        self.Close()


    def updateProgress(self, msg):
        """"""
        if newWx:
            method, count = msg
        else:
            method, count = msg.data
        self.progress_count = count
        try:
            self.progress.SetValue(self.progress_count)
        except:
            pass

    def setProgressRange(self, msg):
        """"""
        if newWx:
            count = msg
        else:
            count = msg.data
        try:
            self.progress.SetRange(count)
        except:
            pass

    

    def updateStatus(self, msg):
        """"""
        if newWx:
            method, text = msg 
        else:
            method, text = msg.data
        self.statusBar.SetStatusText(text)

    
    def saveHistogram(self, msg):
        if newWx:
            data, orf, path, delta = msg
        else:
            data, orf, path, delta = msg.data
        n, bins, patches = plt.hist(data, normed=1, facecolor='c', alpha=0.75, bins=100)
        plt.xlabel('Delta Sum')
        plt.ylabel('Probability')
        plt.title('%s - Histogram of Delta Sum' % orf)
        plt.axvline(delta, color='r', linestyle='dashed', linewidth=3)
        plt.grid(True)
        genePath = os.path.join(path, orf +".png")
        plt.savefig(genePath)
        plt.clf()



    def addFile(self, data):
        if not newWx:
            data = data.data
        fullpath = data["path"]
        name = transit_tools.basename(fullpath)
        type = data["type"]
        date = data["date"]
        self.list_files.InsertStringItem(self.index_file, name)
        self.list_files.SetStringItem(self.index_file, 1, "%s" % type)
        self.list_files.SetStringItem(self.index_file, 2, "%s" % (date))
        self.list_files.SetStringItem(self.index_file, 3, "%s" % (fullpath))
        self.index_file+=1
        

    def finishRun(self,msg):
        if not newWx: msg = msg.data
        try:
            #methods[msg].gui.Enable()
            pass
        except Exception as e:
            print transit_prefix, "Error:", e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def ResetProgress(self):
        self.progress_count = 0



    def HideAllOptions(self):
        self.HideGlobalOptions()
        self.HideProgressSection()
        for name in methods:
            methods[name].gui.Hide()


    def HideGlobalOptions(self):
        self.globalLabel.Hide()
        self.globalNTerminusLabel.Hide()
        self.globalCTerminusLabel.Hide()
        self.globalNTerminusText.Hide()
        self.globalCTerminusText.Hide()


    def ShowGlobalOptions(self):
        self.globalLabel.Show()
        self.globalNTerminusLabel.Show()
        self.globalCTerminusLabel.Show()
        self.globalNTerminusText.Show()
        self.globalCTerminusText.Show()

    def HideProgressSection(self):
        self.progressLabel.Hide()
        self.progress.Hide()


    def ShowProgressSection(self):
        self.progressLabel.Show()
        self.progress.Show()

    def onHimar1Checked(self, event):
        if self.methodCheckBoxHimar1.GetValue():
            self.transposons.append("himar1")
        else:
            self.transposons.remove("himar1")
        self.filterMethodsByTransposon()

    def onTn5Checked(self, event):
        if self.methodCheckBoxTn5.GetValue():
            self.transposons.append("tn5")
        else:
            self.transposons.remove("tn5")
        self.filterMethodsByTransposon()

    def filterMethodsByTransposon(self):
        newmethods = {}
        fullmethods = transit.analysis.methods
        goodTn = False
        for method in fullmethods:
            goodTn = False
            for tn in self.transposons:
                if tn in fullmethods[method].transposons:
                    goodTn = True
            if goodTn:
                newmethods[method] = fullmethods[method]
        
        methodChoiceChoices = [ "[Choose Method]"]
        for name in newmethods:
            methodChoiceChoices.append(methods[name].fullname())
        self.methodChoice.SetItems(methodChoiceChoices)
        self.methodChoice.SetSelection( 0 )



    def SaveFile(self, DIR=None, FILE="", WC=""):
        """
        Create and show the Save FileDialog
        """
        path = ""

        if not DIR:
            DIR = os.getcwd()

        dlg = wx.FileDialog(
            self, message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE, wildcard=WC, style=wx.SAVE|wx.OVERWRITE_PROMPT
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if self.verbose:
                print transit_prefix, "You chose the following output filename: %s" % path
        dlg.Destroy()
        return path

    def OpenFile(self, DIR=".", FILE="", WC=""):
        """
        Create and show the Open FileDialog
        """
        path = ""
        dlg = wx.FileDialog(
            self, message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE, wildcard=WC, style=wx.OPEN
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if self.verbose:
                print transit_prefix, "You chose the following file: %s" % path
        dlg.Destroy()
        return path



    def ShowMessage(self, MSG=""):
        wx.MessageBox(MSG, 'Info', 
            wx.OK | wx.ICON_INFORMATION)
        
    def ShowAskWarning(self, MSG=""):
        dial = wx.MessageDialog(None, MSG, 'Warning', 
            wx.OK | wx.CANCEL | wx.ICON_EXCLAMATION)
        return dial.ShowModal()
    
    def ShowError(self, MSG=""):
        dial = wx.MessageDialog(None, MSG, 'Error', 
            wx.OK | wx.ICON_ERROR)
        dial.ShowModal()
     
    def ctrlSelected(self, col=5):
        selected_ctrl = []
        current = -1
        while True:
            next = self.list_ctrl.GetNextSelected(current)
            if next == -1:
                break
            path = self.list_ctrl.GetItem(next, col).GetText()
            selected_ctrl.append(path)
            current = next
        return selected_ctrl


    def expSelected(self, col=5):
        selected_exp = []
        current = -1
        while True:
            next = self.list_exp.GetNextSelected(current)
            if next == -1:
                break
            path = self.list_exp.GetItem(next, col).GetText()
            selected_exp.append(path)
            current = next
        return selected_exp


    def allSelected(self, col=5):
        selected_all = self.ctrlSelected(col) + self.expSelected(col)
        return selected_all


    def ctrlAll(self, col=5):
        all_ctrl = []
        for i in range(self.list_ctrl.GetItemCount()):
            all_ctrl.append(self.list_ctrl.GetItem(i, col).GetText())
        return all_ctrl


    def expAll(self, col=5):
        all_exp = []
        for i in range(self.list_exp.GetItemCount()):
            all_exp.append(self.list_exp.GetItem(i, col).GetText())
        return all_exp


    def loadCtrlFileFunc(self, event):
        try:
        
            dlg = wx.FileDialog(
                self, message="Choose a file",
                defaultDir=self.workdir,
                defaultFile="",
                wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
                )
            if dlg.ShowModal() == wx.ID_OK:
                paths = dlg.GetPaths()
                print "You chose the following Control file(s):"
                for fullpath in paths:
                    print "\t%s" % fullpath
                    name = transit_tools.basename(fullpath)
                    (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = transit_tools.get_wig_stats(fullpath)
                    self.list_ctrl.InsertStringItem(self.index_ctrl, name)
                    self.list_ctrl.SetStringItem(self.index_ctrl, 1, "%1.1f" % (totalrd))
                    self.list_ctrl.SetStringItem(self.index_ctrl, 2, "%2.1f" % (density*100))
                    self.list_ctrl.SetStringItem(self.index_ctrl, 3, "%1.1f" % (meanrd))
                    self.list_ctrl.SetStringItem(self.index_ctrl, 4, "%d" % (maxrd))
                    self.list_ctrl.SetStringItem(self.index_ctrl, 5, "%s" % (fullpath))
                    self.index_ctrl+=1
            dlg.Destroy()
        except Exception as e:
            print transit_prefix, "Error:", e
            print "PATH", fullpath
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def loadExpFileFunc(self, event):
        try:

            dlg = wx.FileDialog(
                self, message="Choose a file",
                defaultDir=self.workdir,
                defaultFile="",
                wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
                )
            if dlg.ShowModal() == wx.ID_OK:
                paths = dlg.GetPaths()
                print "You chose the following Experimental file(s):"
                for fullpath in paths:
                    print "\t%s" % fullpath
                    name = transit_tools.basename(fullpath)
                    (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = transit_tools.get_wig_stats(fullpath)
                    self.list_exp.InsertStringItem(self.index_exp, name)
                    self.list_exp.SetStringItem(self.index_exp, 1, "%1.1f" % (totalrd))
                    self.list_exp.SetStringItem(self.index_exp, 2, "%2.1f" % (density*100))
                    self.list_exp.SetStringItem(self.index_exp, 3, "%1.1f" % (meanrd))
                    self.list_exp.SetStringItem(self.index_exp, 4, "%d" % (maxrd))
                    self.list_exp.SetStringItem(self.index_exp, 5, "%s" % (fullpath))
                    self.index_exp+=1
            dlg.Destroy()
        except Exception as e:
            print transit_prefix, "Error:", e
            print "PATH", fullpath
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def ctrlRemoveFunc(self, event):
        next = self.list_ctrl.GetNextSelected(-1)
        while next != -1:
            if self.verbose:
                print transit_prefix, "Removing control item (%d): %s" % (next, self.list_ctrl.GetItem(next, 0).GetText())
            self.list_ctrl.DeleteItem(next)
            next = self.list_ctrl.GetNextSelected(-1)
            self.index_ctrl-=1 

     


    def expRemoveFunc(self, event):
        next = self.list_exp.GetNextSelected(-1)
        while next != -1:
            if self.verbose:
                print transit_prefix, "Removing experimental item (%d): %s" % (next, self.list_exp.GetItem(next, 0).GetText())
            self.list_exp.DeleteItem(next)
            next = self.list_exp.GetNextSelected(-1)
            self.index_exp-=1


    def allViewFunc(self, event, gene=""):
        
        annotationpath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()

        if datasets and annotationpath:
            if self.verbose:
                print transit_prefix, "Visualizing counts for:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            viewWindow = trash.TrashFrame(self, datasets, annotationpath, gene)
            viewWindow.Show()
        elif not datasets:
            self.ShowError("Error: No datasets selected.")
            return
        else:
            self.ShowError("Error: No annotation file selected.")
            return





    def ctrlViewFunc(self, event, gene=""):
        annotationpath = self.annotation
        datasets = self.ctrlSelected()
        if datasets and annotationpath:
            if self.verbose:
                print transit_prefix, "Visualizing counts for:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            viewWindow = trash.TrashFrame(self, datasets, annotationpath, gene)
            viewWindow.Show()
        elif not datasets:
            if self.verbose:
                print transit_prefix, "No datasets selected to visualize!"
        else:
            if self.verbose:
                print transit_prefix, "No annotation file selected"



    def expViewFunc(self, event, gene=""):
        annotationpath = self.annotation
        datasets = self.expSelected()
        if datasets and annotationpath:
            if self.verbose:
                print transit_prefix, "Visualizing counts for:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            viewWindow = trash.TrashFrame(self, datasets, annotationpath, gene)
            viewWindow.Show()
        elif not datasets:
            if self.verbose:
                print transit_prefix, "No datasets selected to visualize!"
        else:
            if self.verbose:
                print transit_prefix, "No annotation file selected"



    def scatterFunc(self, event):
        """ """
        #annotationpath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()
        if len(datasets) == 2:
            if self.verbose:
                print transit_prefix, "Showing scatter plot for:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            (data, position) = transit_tools.get_data(datasets)
            X = data[0,:]
            Y = data[1,:]

            plt.plot(X,Y, "bo")
            plt.title('Scatter plot - Reads at TA sites')
            plt.xlabel(transit_tools.fetch_name(datasets[0]))
            plt.ylabel(transit_tools.fetch_name(datasets[1])) 
            plt.show()
        else:
            self.ShowError(MSG="Please make sure only two datasets are selected (across control and experimental datasets).")


    def qcFunc(self, event):
        datasets = self.allSelected()
        nfiles = len(datasets)

        if nfiles == 0:
            print transit_prefix, "You must select atleast one dataset control or experimental dataset."
            self.ShowError(MSG="You must select atleast one dataset control or experimental dataset.")
            
        elif nfiles >= 1:
            print transit_prefix, "Displaying results:", datasets[0]
            try:
                qcWindow = qcDisplay.qcFrame(self, datasets)
                qcWindow.Show()
            except Exception as e:
                print transit_prefix, "Error occured displaying file:", str(e)
                traceback.print_exc()



    def annotationFileFunc(self, event):

        wc = u"Prot Table (*.prot_table)|*.prot_table;|\nProt Table (*.txt)|*.txt;|\nProt Table (*.dat)|*.dat;|\nAll files (*.*)|*.*" 
        self.annotation = self.OpenFile(DIR=self.workdir, FILE="", WC=wc)
        if self.annotation:
            self.annotationFilePicker.SetLabel(transit_tools.basename(self.annotation))
            if self.verbose:
                print transit_prefix, "Annotation File Selected:", self.annotation
        else:
            self.annotationFilePicker.SetLabel("[Click to add Annotation File (.prot_table)]")
        
    def MethodSelectFunc(self, selected_name, test=""):
        #X = self.methodChoice.GetCurrentSelection()
        #selected_name = self.methodChoice.GetString(X)

        #If empty is selected
        if selected_name == "[Choose Method]":
            self.HideAllOptions()
            self.methodInfoText.SetLabel(u"Instructions")
            self.methodInstructions.Show()
            self.methodInstructions.SetLabel(self.instructions_text)
            self.methodInstructions.Wrap(method_wrap_width)
            self.methodShortText.Hide()
            self.methodLongText.Hide()
            self.methodDescText.Hide()

            self.method_choice = ""
        else:
            self.ShowGlobalOptions()
            #Show Selected Method and hide Others
            for name in methods:
                methods[name].gui.Hide()
                if methods[name].fullname() == selected_name:
                    text = """%s
                       
    %s 
                    """ % (methods[name].fullname(), methods[name].description)

                    self.methodInfoText.SetLabel("%s" % methods[name].short_name)
                    
                    #self.methodShortText.SetLabel("[%s]" % methods[name].short_name)
                    #self.methodShortText.Wrap(250)
                    #self.methodLongText.SetLabel(methods[name].long_name)
                    #self.methodLongText.Wrap(250)
                    self.methodDescText.SetLabel(methods[name].description)
                    self.methodDescText.Wrap(250)
                    self.methodDescText.Show()
                    self.methodInstructions.SetLabel("")
                    methods[name].gui.Show()
                else:
                    methods[name].gui.Hide()
            self.ShowProgressSection()
            self.method_choice = selected_name

        self.Layout()
        if self.verbose:
            print transit_prefix, "Selected Method: %s" % (selected_name)




    def displayFileFunc(self, event):
        next = self.list_files.GetNextSelected(-1)
        if next > -1:
            dataset = self.list_files.GetItem(next, 3).GetText()
            if self.verbose:
                print transit_prefix, "Displaying results:", self.list_files.GetItem(next, 0).GetText()

            try:
                #fileWindow = fileDisplay.FileFrame(self, dataset, self.list_files.GetItem(next, 1).GetText())
                fileWindow = fileDisplay.TransitGridFrame(self, dataset)
                fileWindow.Show()
            except Exception as e:
                print transit_prefix, "Error occurred displaying file:", str(e)
                traceback.print_exc()

        else:
            if self.verbose:
                print transit_prefix, "No results selected to display!"


    def graphFileFunc(self, event):
        # 0 - nothing
        # 1 - Volcano
        # 2 - Hist gene counts ratio
        
        plot_choice  = self.graphFileChoice.GetCurrentSelection()
        plot_name = self.graphFileChoice.GetString(plot_choice)
        if plot_choice == 0:
                return
        next = self.list_files.GetNextSelected(-1)
        if next > -1:
            dataset_path = self.list_files.GetItem(next, 3).GetText()
            dataset_name = self.list_files.GetItem(next, 0).GetText()
            dataset_type = self.list_files.GetItem(next, 1).GetText()
            
            if self.verbose:
                print transit_prefix, "Creating a", plot_name, " for dataset", dataset_name

            if plot_choice == 1:
                self.graphVolcanoPlot(dataset_name, dataset_type, dataset_path)
            elif plot_choice == 2:
                self.graphGeneCounts(dataset_name, dataset_type, dataset_path)

            self.graphFileChoice.SetSelection(0)
            
        else:
            self.ShowError(MSG="Please select a results file to plot!")
    
        


    def graphGeneCounts(self, dataset_name, dataset_type, dataset_path):
        try:
            if dataset_type == "Resampling":
                X = []
                for line in open(dataset_path):
                    if line.startswith("#"): continue
                    tmp = line.strip().split("\t")
                    try:
                        log2FC = float(tmp[8])
                    except:
                        log2FC = 0
                    X.append(log2FC)

                n, bins, patches = plt.hist(X, normed=1, facecolor='c', alpha=0.75, bins=100)
                plt.xlabel('log2 FC - Total Gene Counts')
                plt.ylabel('Probability')
                plt.title('Histogram of log2 Fold Change for Total Normalized Counts within Genes')
                plt.axvline(0, color='r', linestyle='dashed', linewidth=3)
                plt.grid(True)
                plt.show()
                plt.close()
            else:
               self.ShowError(MSG="Need to select a 'Resampling' results file for this type of plot.")

        except Exception as e:
            print transit_prefix, "Error occurred creating plot:", str(e)
            traceback.print_exc()





    def graphVolcanoPlot(self, dataset_name, dataset_type, dataset_path):
        try:
            if dataset_type == "Resampling":
                X = []; Y = [];
                for line in open(dataset_path):
                    if line.startswith("#"): continue
                    tmp = line.strip().split("\t")
                    try:
                        #log2FC = math.log(float(tmp[6])/float(tmp[5]),2)
                        log2FC = float(tmp[8])
                        log10qval = -math.log(float(tmp[-1].strip()), 10)
                    except:
                        log2FC = 0
                        log10qval = 0
                        #log2FC = 1
                        #log10qval = 1
                    X.append(log2FC)
                    Y.append(log10qval)

                plt.plot(X,Y, "bo")
                plt.xlabel("Log Fold Change (base 2)")
                plt.ylabel("-Log q-value (base 10)")
                plt.title("Resampling - Volcano plot")
                plt.show()
                plt.close()
            else:
               self.ShowError(MSG="Need to select a 'Resampling' results file for this type of plot.")
            
        except Exception as e:
            print "Error occurred creating plot:", str(e)
        



    def LoessPrevFunc(self,event):
        datasets_selected = self.ctrlSelected() + self.expSelected()
        if not datasets_selected:
            self.ShowError(MSG="Need to select at least one control or experimental dataset.")
            return

        data, position = transit_tools.get_data(datasets_selected)
        (K,N) = data.shape
        window = 100
        for j in range(K):

            size = len(position)/window + 1
            x_w = numpy.zeros(size)
            y_w = numpy.zeros(size)
            for i in range(size):
                x_w[i] = window*i
                y_w[i] = sum(data[j][window*i:window*(i+1)])

            y_smooth = transit_tools.loess(x_w, y_w, h=10000)
            plt.plot(x_w, y_w, "g+")
            plt.plot(x_w, y_smooth, "b-")
            plt.xlabel("Genomic Position")
            plt.ylabel("Reads per 100 insertion sites")

            plt.title("LOESS Fit - %s" % transit_tools.basename(datasets_selected[j]) )
            plt.show()
            plt.close()
    

    def addFileFunc(self, event):

        try:
            dlg = wx.FileDialog(
                self, message="Choose a file",
                defaultDir=self.workdir,
                defaultFile="",
                wildcard=u"Results Files (*.dat)|*.dat;|\nResults Files (*.txt)|*.txt;|\nAll files (*.*)|*.*",
                style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
                )
            if dlg.ShowModal() == wx.ID_OK:
                paths = dlg.GetPaths()
                print "You chose the following Results file(s):"
                for fullpath in paths:
                    print "\t%s" % fullpath
                    name = transit_tools.basename(fullpath)
                    line = open(fullpath).readline()
                    if line.startswith("#Gumbel"):
                        type = "Gumbel"
                    elif line.startswith("#Binomial"):
                        type = "Binomial"
                    elif line.startswith("#HMM - Sites"):
                        type = "HMM - Sites"
                    elif line.startswith("#HMM - Genes"):
                        type = "HMM - Genes"
                    elif line.startswith("#Resampling"):
                        type = "Resampling"
                    elif line.startswith("#DE-HMM - Sites"):
                        type = "DE-HMM - Sites"
                    elif line.startswith("#DE-HMM - Segments"):
                        type = "DE-HMM - Segments"
                    else:
                        type = "Unknown"
                    data = {"path":fullpath, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
                    if newWx:
                        wx.CallAfter(pub.sendMessage, "file", data=data)
                    else:
                        wx.CallAfter(pub.sendMessage, "file", data)
            dlg.Destroy()
        except Exception as e:
            print transit_prefix, "Error:", e
            print "PATH", fullpath
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def choseMethodsMenu(self, selected_name, event):
        if self.verbose:
            print transit_prefix, "Selected Method: %s" % (selected_name)
        self.MethodSelectFunc(selected_name) 


    def convertToIGV(self, dataset_list, annotationPath, path):
        fulldata = []
        position = []
        for j,dataset in enumerate(dataset_list):
            temp = []
            for line in open(dataset):
                if line.startswith("#"): continue
                if line.startswith("variable"): continue
                if line.startswith("location"): continue
                tmp = line.split()
                pos = int(tmp[0])
                read = int(tmp[1])
                temp.append(read)
                if j == 0:
                    position.append(pos)
            fulldata.append(temp)
                

        output = open(path, "w")
        output.write("#Converted to IGV with TRANSIT.\n")
        output.write("#Files:\n#%s\n" % "\n#".join(dataset_list))
        output.write("#Chromosome\tStart\tEnd\tFeature\t%s\tTAs\n" % ("\t".join([transit_tools.fetch_name(D) for D in dataset_list])))
        chrom = transit_tools.fetch_name(annotationPath)

        for i,pos in enumerate(position):
            output.write("%s\t%s\t%s\tTA%s\t%s\t1\n" % (chrom, position[i], position[i]+1, position[i], "\t".join(["%s" % fulldata[j][i] for j in range(len(fulldata))])))
        output.close()




    def ctrlToIGV(self, event):
        annotationPath = self.annotation
        datasets = self.ctrlSelected()
        defaultFile = "read_counts.igv"
        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        
        outputPath = self.SaveFile(defaultDir, defaultFile)
        if datasets and annotationPath and outputPath:
            if self.verbose:
                print transit_prefix, "Converting the following datasets to IGV format:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            self.convertToIGV(datasets, annotationPath, outputPath)
            if self.verbose:
                print transit_prefix, "Finished conversion"
        elif not datasets:
            self.ShowError("Error: No datasets selected to convert!")
        elif not annotationPath:
            self.ShowError("Error: No annotation file selected.")
        else:
            pass
        
        

    def expToIGV(self, event):
        annotationPath = self.annotation
        datasets = self.expSelected()
        defaultFile = "read_counts.igv"
        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        outputPath = self.SaveFile(defaultDir, defaultFile)
        if datasets and annotationPath and outputPath:
            if self.verbose:
                print transit_prefix, "Converting the following datasets to IGV format:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            self.convertToIGV(datasets, annotationPath, outputPath)
            if self.verbose:
                print transit_prefix, "Finished conversion"
        elif not datasets:
            self.ShowError("Error: No datasets selected to convert!")
        elif not annotationPath:
            self.ShowError("Error: No annotation file selected.")
        else:
            pass
       
 

    def allToIGV(self, event):
        annotationPath = self.annotation
        datasets = self.ctrlSelected() + self.expSelected()
        defaultFile = "read_counts.igv"
        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        outputPath = self.SaveFile(defaultDir, defaultFile)
        if datasets and annotationPath and outputPath:
            if self.verbose:
                print transit_prefix, "Converting the following datasets to IGV format:", ", ".join([transit_tools.fetch_name(d) for d in datasets])
            self.convertToIGV(datasets, annotationPath, outputPath)
            if self.verbose:
                print transit_pefix, "Finished conversion"
        elif not datasets:
            self.ShowError("Error: No datasets selected to convert!")
        elif not annotationPath:
            self.ShowError("Error: No annotation file selected.")
        else:
            pass


        
    def annotationPT_to_GFF3(self, event):
        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".gff3"
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        outputPath = self.SaveFile(defaultDir, defaultFile)

        ORGANISM = transit_tools.fetch_name(annotationpath)
        if not annotationpath:
            self.ShowError("Error: No annotation file selected.")

        elif outputPath:
            if self.verbose:
                print transit_prefix, "Converting annotation file from prot_table format to GFF3 format"
            year = time.localtime().tm_year
            month = time.localtime().tm_mon
            day = time.localtime().tm_mday
            
            output = open(outputPath, "w")
            output.write("##gff-version 3\n")
            output.write("##converted to IGV with TRANSIT.\n")
            output.write("##date %d-%d-%d\n" % (year, month, day))
            output.write("##Type DNA %s\n" % ORGANISM)
            
            for line in open(annotationpath):
                if line.startswith("#"): continue
                tmp = line.strip().split("\t")
                desc = tmp[0]; start = int(tmp[1]); end = int(tmp[2]); strand = tmp[3];
                length = tmp[4]; name = tmp[7]; orf = tmp[8];
                ID = name
                desc.replace("%", "%25").replace(";", "%3B").replace("=","%3D").replace(",","%2C")
                output.write("%s\tRefSeq\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;Alias=%s;locus_tag=%s;desc=%s\n" % (ORGANISM, start, end, strand, orf,ID, orf, orf,desc))
                
            output.close()
            if self.verbose:
                print transit_prefix, "Finished conversion"        
       



    def annotationPT_to_PTT(self, event):
 
        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".ptt.table"
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            self.ShowError("Error: No annotation file selected.")
        elif not datasets:
            self.ShowError("Error: Please add a .wig dataset, to determine TA sites.")            
        else:
            
            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath: return
            if self.verbose:
                print transit_prefix, "Converting annotation file from prot_table format to PTT format"
            (data, position) = transit_tools.get_data(datasets)
            orf2info = transit_tools.get_gene_info(annotationpath)
            hash = transit_tools.get_pos_hash(annotationpath)
            (orf2reads, orf2pos) = transit_tools.get_gene_reads(hash, data, position, orf2info)

            output = open(outputPath, "w")
            output.write("geneID\tstart\tend\tstrand\tTA coordinates\n")
            for line in open(annotationpath):
                if line.startswith("#"): continue
                tmp = line.strip().split("\t")
                orf = tmp[8]
                name = tmp[7]
                desc = tmp[0]
                start = int(tmp[1])
                end = int(tmp[2])
                strand = tmp[3]
                ta_str = "no TAs"
                if orf in orf2pos: ta_str = "\t".join([str(int(ta)) for ta in orf2pos[orf]])
                output.write("%s\t%s\t%s\t%s\t%s\n" % (orf, start, end, strand, ta_str))
            output.close()
            if self.verbose:
                print transit_prefix, "Finished conversion"
                


    def annotationPTT_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()
        
        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            self.ShowError("Error: No annotation file selected.")
        #elif not datasets:
        #    self.ShowError("Error: Please add a .wig dataset, to determine TA sites.")
        else:

            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath: return
            if self.verbose:
                print transit_prefix, "Converting annotation file from PTT format to prot_table format"
            #(data, position) = transit_tools.get_data(datasets)
            #orf2info = transit_tools.get_gene_info(annotationpath)
            #hash = transit_tools.get_pos_hash(annotationpath)
            #(orf2reads, orf2pos) = transit_tools.get_gene_reads(hash, data, position, orf2info)
            
            
            output = open(outputPath, "w")
            #output.write("geneID\tstart\tend\tstrand\tTA coordinates\n")
            for line in open(annotationpath):
                if line.startswith("#"): continue
                if line.startswith("geneID"): continue
                tmp = line.strip().split("\t")
                orf = tmp[0]
                if orf == "intergenic": continue
                name = "-"
                desc = "-"
                start = int(tmp[1])
                end = int(tmp[2])
                length = ((end-start+1)/3)-1
                strand = tmp[3]
                someID = "-"
                someID2 = "-"
                COG = "-"
                output.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" % (desc, start, end, strand, length, someID, someID2, name, orf, COG))
            output.close()
            if self.verbose:
                print transit_prefix, "Finished conversion"

        #geneID  start   end strand  TA coordinates
        


    def annotationGFF3_to_PT(self, event):

        annotationpath = self.annotation
        defaultFile = transit_tools.fetch_name(annotationpath) + ".prot_table"
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        datasets = self.ctrlSelected() + self.expSelected()
        if not annotationpath:
            self.ShowError("Error: No annotation file selected.")
        else:
            outputPath = self.SaveFile(defaultDir, defaultFile)
            if not outputPath: return
            if self.verbose:
                print transit_prefix, "Converting annotation file from GFF3 format to prot_table format"

            output = open(outputPath, "w")
            for line in open(annotationpath):
                if line.startswith("#"): continue
                tmp = line.strip().split("\t")
                chr = tmp[0]
                type = tmp[2]
                start = int(tmp[3])
                end = int(tmp[4])
                length = ((end-start+1)/3)-1
                strand = tmp[6]
                features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
                if "ID" not in features: continue
                orf = features["ID"]
                name = features.get("Name", "-")
                if name == "-": name = features.get("name", "-")
                
                desc = features.get("Description", "-")
                if desc == "-": desc = features.get("description", "-")
                if desc == "-": desc = features.get("Desc", "-")
                if desc == "-": desc = features.get("desc", "-")

                someID = "-"
                someID2 = "-"
                COG = "-"
                output.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" % (desc, start, end, strand, length, someID, someID2, name, orf, COG))
            output.close()
            if self.verbose:
                print transit_prefix, "Finished conversion"

                 

    def RunMethod(self, event):
        #FLORF
        #X = self.methodChoice.GetCurrentSelection()
        #selected_name = self.methodChoice.GetString(X)
        selected_name = self.method_choice
        for name in methods:
            if  methods[name].fullname() == selected_name:
                methodobj = methods[name].method
        try:
            M = methodobj.fromGUI(self)
            if M: 
                thread = threading.Thread(target=M.Run())
                thread.setDaemon(True)
                thread.start()
        except Exception as e:
            print transit_prefix, "Error:", str(e)
            traceback.print_exc()







def msg():
    return """transit.py {gumbel, hmm, resampling} [-h/--help]

    
Help:  
    """


if __name__ == "__main__":


    #If no arguments, show GUI:
    if len(sys.argv) == 1:
   
        print transit_prefix, "Running in GUI Mode"
         
        app = wx.App(False)

        #create an object of CalcFrame
        frame = TnSeekFrame(None)
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

        else:
            methodobj = methods[method_name].method.fromconsole()
            methodobj.Run()            



