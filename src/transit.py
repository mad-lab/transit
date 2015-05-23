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

if len(sys.argv) == 1:
    import wx
    from wx.lib.pubsub import pub
   
import os
import time
import datetime
import ntpath
import threading
import numpy
import matplotlib.pyplot as plt
import multiprocessing as mp
import math

# trash view stuff
import trash
import transit_gui
import transit_tools
import fileDisplay

import gumbelMH
import hmm_geom
import hmm_tools
import resampling

import imgTRANSIT

import argparse


wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"
transit_prefix = "[TRANSIT]"

#inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TnSeekFrame(transit_gui.MainFrame):
    #constructor
    def __init__(self,parent):
        #initialize parent class
        transit_gui.MainFrame.__init__(self,parent)


        self.logoImg.SetBitmap(imgTRANSIT.TRANSIT.GetImage().ConvertToBitmap())

        self.index_ctrl = 0
        self.list_ctrl.InsertColumn(0, 'File', width=210)
        self.list_ctrl.InsertColumn(1, 'Total Reads', width=85)
        self.list_ctrl.InsertColumn(2, 'Density', width=85)
        self.list_ctrl.InsertColumn(3, 'Mean Count', width=85)
        self.list_ctrl.InsertColumn(4, 'Max Count', width=85)
        self.list_ctrl.InsertColumn(5, 'Full Path', width=403)


        self.index_exp = 0
        self.list_exp.InsertColumn(0, 'File', width=210)
        self.list_exp.InsertColumn(1, 'Total Reads', width=85)
        self.list_exp.InsertColumn(2, 'Density', width=85)
        self.list_exp.InsertColumn(3, 'Mean Count', width=85)
        self.list_exp.InsertColumn(4, 'Max Count', width=85)
        self.list_exp.InsertColumn(5, 'Full Path',width=403)


        self.index_file = 0
        self.list_files.InsertColumn(0, 'Name', width=350)
        self.list_files.InsertColumn(1, 'Type', width=100)
        self.list_files.InsertColumn(2, 'Date', width=230)
        self.list_files.InsertColumn(3, 'Full Path',width=403)


        self.verbose = True


        self.progress_count = 0
        self.gumbel_count = 0
        self.hmm_count = 0
        self.resampling_count = 0
        pub.subscribe(self.updateGumbel, "gumbel")
        pub.subscribe(self.updateHMM, "hmm")
        pub.subscribe(self.updateResampling, "resampling")
        pub.subscribe(self.updateProgress, "update")
        pub.subscribe(self.addFile, "file")
        pub.subscribe(self.finishRun, "finish")
   
 
        #self.outputDirPicker.SetPath(os.path.dirname(os.path.realpath(__file__)))

        self.progress.Hide()
        self.progressLabel.Hide()
        self.HideGlobalOptions()
        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()

    def Exit(self, event):
        """Exit Menu Item"""
        if self.verbose:
            print transit_prefix, "Exiting Transit"
        self.Close()


    def updateGumbel(self, msg):
        self.gumbel_count+=1
        self.gumbelProgress.SetValue(self.gumbel_count)

        X = self.methodChoice.GetCurrentSelection()
        if X == 1: # if gumbel
            if self.gumbel_count > self.gumbelProgress.GetRange(): msg = ""
            self.statusBar.SetStatusText(msg)
    
    def updateHMM(self, msg):
        self.hmm_count+=1
        self.hmmProgress.SetValue(self.hmm_count)
        
        X = self.methodChoice.GetCurrentSelection()
        if X == 2: # if hmm
            if self.hmm_count > self.hmmProgress.GetRange(): msg = ""
            self.statusBar.SetStatusText(msg)


    def updateResampling(self, msg):
        self.resampling_count+=1
        self.resamplingProgress.SetValue(self.resampling_count)

        X = self.methodChoice.GetCurrentSelection()
        if X == 3: # if resampling
            if self.resampling_count > self.resamplingProgress.GetRange(): msg = ""
            self.statusBar.SetStatusText(msg)


    def updateProgress(self, msg):
        """"""
        self.progress_count += 1
        self.progress.SetValue(self.progress_count)
        self.statusBar.SetStatusText(msg)
        if self.progress_count > self.progress.GetRange():
            self.statusBar.SetStatusText("Finished!")
        #self.statusBar.SetStatusText("Running Gumbel Method... %2.0f%%" % (100.0*self.progress_count/self.progress.GetRange()))


    def addFile(self, data):
        fullpath = data["path"]
        name = ntpath.basename(fullpath)
        type = data["type"]
        date = data["date"]
        self.list_files.InsertStringItem(self.index_file, name)
        self.list_files.SetStringItem(self.index_file, 1, "%s" % type)
        self.list_files.SetStringItem(self.index_file, 2, "%s" % (date))
        self.list_files.SetStringItem(self.index_file, 3, "%s" % (fullpath))
        self.index_file+=1
        

    def finishRun(self,msg):
        if msg == "gumbel":
            self.gumbelButton.Enable()
        elif msg == "hmm":
            self.hmmButton.Enable()
        elif msg == "resampling":
            self.resamplingButton.Enable()


    def ResetProgress(self):
        self.progress_count = 0


    def SetProgressRange(self, X):
        self.progress.SetRange(X)


    def HideAllOptions(self):

        self.HideGlobalOptions()
        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()

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


    def HideGumbelOptions(self):
        self.gumbelLabel.Hide()
        self.gumbelInstructions.Hide()
        self.gumbelSampleLabel.Hide()
        self.gumbelSampleText.Hide()
        self.gumbelBurninLabel.Hide()
        self.gumbelBurninText.Hide()
        self.gumbelTrimLabel.Hide()
        self.gumbelTrimText.Hide()
        self.gumbelReadLabel.Hide()
        self.gumbelReadChoice.Hide()
        self.gumbelRepLabel.Hide()
        self.gumbelRepChoice.Hide()
        self.gumbelButton.Hide()
        self.gumbelProgress.Hide()


    def ShowGumbelOptions(self):
        self.gumbelLabel.Show()
        self.gumbelInstructions.Show()
        self.gumbelSampleLabel.Show()
        self.gumbelSampleText.Show()
        self.gumbelBurninLabel.Show()
        self.gumbelBurninText.Show()
        self.gumbelTrimLabel.Show()
        self.gumbelTrimText.Show()
        self.gumbelReadLabel.Show()
        self.gumbelReadChoice.Show()
        self.gumbelRepLabel.Show()
        self.gumbelRepChoice.Show()
        self.gumbelButton.Show()
        self.gumbelProgress.Show()


    def HideHMMOptions(self):
        self.hmmLabel.Hide()
        self.hmmInstructions.Hide()
        self.hmmRepLabel.Hide()
        self.hmmRepChoice.Hide()
        self.hmmButton.Hide()
        self.hmmProgress.Hide()

    
    def ShowHMMOptions(self):
        self.hmmLabel.Show()
        self.hmmInstructions.Show()
        self.hmmRepLabel.Show()
        self.hmmRepChoice.Show()
        self.hmmButton.Show()
        self.hmmProgress.Show()


    def HideResamplingOptions(self):
        self.resamplingLabel.Hide()
        self.resamplingInstructions.Hide()
        self.resamplingSampleLabel.Hide()
        self.resamplingSampleText.Hide()
        self.resamplingNormLabel.Hide()
        self.resamplingNormChoice.Hide()
        self.resamplingHistCheck.Hide()
        self.resamplingAdaptiveCheck.Hide()
        self.resamplingButton.Hide()
        self.resamplingProgress.Hide()
    

    def ShowResamplingOptions(self):
        self.resamplingLabel.Show()
        self.resamplingInstructions.Show()
        self.resamplingSampleLabel.Show()
        self.resamplingSampleText.Show()
        self.resamplingNormLabel.Show()
        self.resamplingNormChoice.Show()
        self.resamplingHistCheck.Show()
        self.resamplingAdaptiveCheck.Show()
        self.resamplingButton.Show()
        self.resamplingProgress.Show()
 

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
            fullpath = self.ctrlFilePicker.GetPath()
            name = ntpath.basename(fullpath)
            totalrd,density,meanrd,maxrd = transit_tools.get_wig_stats(fullpath)
            self.list_ctrl.InsertStringItem(self.index_ctrl, name)
            self.list_ctrl.SetStringItem(self.index_ctrl, 1, "%d" % (totalrd))
            self.list_ctrl.SetStringItem(self.index_ctrl, 2, "%2.1f" % (density*100))
            self.list_ctrl.SetStringItem(self.index_ctrl, 3, "%1.1f" % (meanrd))
            self.list_ctrl.SetStringItem(self.index_ctrl, 4, "%d" % (maxrd))
            self.list_ctrl.SetStringItem(self.index_ctrl, 5, "%s" % (fullpath))
            if self.verbose:
                print transit_prefix, "Adding control item (%d): %s" % (self.index_ctrl, name)
            self.index_ctrl+=1
        except e:
            print transit_prefix, "Error:", e

        


    def loadExpFileFunc(self, event):
        try:
            fullpath = self.expFilePicker.GetPath()
            name = ntpath.basename(fullpath)
            totalrd,density,meanrd,maxrd = transit_tools.get_wig_stats(fullpath)
            self.list_exp.InsertStringItem(self.index_exp, name)
            self.list_exp.SetStringItem(self.index_exp, 1, "%d" % (totalrd))
            self.list_exp.SetStringItem(self.index_exp, 2, "%2.1f" % (density*100))
            self.list_exp.SetStringItem(self.index_exp, 3, "%1.1f" % (meanrd))
            self.list_exp.SetStringItem(self.index_exp, 4, "%d" % (maxrd))
            self.list_exp.SetStringItem(self.index_exp, 5, "%s" % (fullpath))
            if self.verbose:
                print transit_prefix, "Adding experimental item (%d): %s" % (self.index_exp, name)
            self.index_exp+=1
        

        except e:
            print transit_prefix, "Error:", e 


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
        
        annotationpath = self.annotationFilePicker.GetPath()
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
        annotationpath = self.annotationFilePicker.GetPath()
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
        annotationpath = self.annotationFilePicker.GetPath()
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
        #annotationpath = self.annotationFilePicker.GetPath()
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



    def annotationFileFunc(self, event):
        if self.verbose:
            print transit_prefix, "Annotation File Selected:", self.annotationFilePicker.GetPath()
        
        
    def MethodSelectFunc(self, event):
        X = self.methodChoice.GetCurrentSelection()


        self.progressLabel.Show()
        if X == 0:
            self.HideAllOptions()
            self.progressLabel.Hide()
            self.mainInstructions.Show()
        elif X == 1:
            self.mainInstructions.Hide()
            self.ShowGlobalOptions()
            self.ShowGumbelOptions()
            self.HideHMMOptions()
            self.HideResamplingOptions()
        elif X == 2:
            self.mainInstructions.Hide()
            self.ShowGlobalOptions()
            self.ShowHMMOptions()
            self.HideGumbelOptions()
            self.HideResamplingOptions()
        elif X == 3:
            self.mainInstructions.Hide()
            self.ShowGlobalOptions()
            self.ShowResamplingOptions()
            self.HideGumbelOptions()
            self.HideHMMOptions()

        self.Layout()

        if self.verbose:
            print transit_prefix, "Selected Method (%d): %s" % (X, self.methodChoice.GetString(X))


    def displayFileFunc(self, event):
        next = self.list_files.GetNextSelected(-1)
        if next > -1:
            dataset = self.list_files.GetItem(next, 3).GetText()
            if self.verbose:
                print transit_prefix, "Displaying results:", self.list_files.GetItem(next, 0).GetText()

            try:
                fileWindow = fileDisplay.FileFrame(self, dataset, self.list_files.GetItem(next, 1).GetText())
                fileWindow.Show()
            except Exception as e:
                print "Error occurred displaying file", e
        else:
            if self.verbose:
                print transit_prefix, "No results selected to display!"


    def graphFileFunc(self, event):
        next = self.list_files.GetNextSelected(-1)
        if next > -1:
            dataset = self.list_files.GetItem(next, 3).GetText()
            if self.verbose:
                print transit_prefix, "Graphing results:", self.list_files.GetItem(next, 0).GetText()

            try:
                filetype = self.list_files.GetItem(next, 1).GetText()
                if filetype == "Resampling":
                    X = []; Y = [];
                    for line in open(dataset):
                        if line.startswith("#"): continue
                        tmp = line.strip().split("\t")
                        try:
                            log2FC = math.log(float(tmp[6])/float(tmp[5]),2)
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

                elif filetype == "Gumbel":
                    X = []; Y = [];
                    for line in open(dataset):
                        if line.startswith("#"): continue
                        tmp = line.strip().split("\t")
                        pprob = float(tmp[7])
                        if pprob < 0: continue
                        X.append(pprob)
                    
                    X.sort()
                    plt.plot(X,"bo")
                    plt.xlabel("Rank")
                    plt.ylabel("Zbar")
                    plt.title("Gumbel - Ranked Posterior Probability of Essentiality (Zbar)")
                    plt.show()



            except Exception as e:
                print "Error occurred displaying file:", e
        else:
            self.ShowError(MSG="No results selected to display!")




    def tempFileFunc(self, event):
        file1= "gumbel_H37Rv_Sassetti_glycerol_s100_b5_t1.dat"
        type="Gumbel"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        if self.verbose:
            print transit_prefix, "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "hmm_H37Rv_Sassetti_glycerol_sites.dat"
        type="HMM - Sites"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        if self.verbose:
            print transit_prefix, "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "hmm_H37Rv_Sassetti_glycerol_genes.dat"
        type="HMM - Genes"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        if self.verbose:
            print transit_prefix, "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "resampling_results_g0g1_g2g3.dat"
        #file1= "resampling_results_qval_test.dat"
        type="Resampling"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        if self.verbose:
            print transit_prefix, "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


    def addFileFunc(self, event):


        try:
            defaultDir = os.getcwd() 
            path = self.OpenFile(defaultDir, "")
            line = open(path).readline()
            if line.startswith("#Gumbel"):
                type = "Gumbel"
            elif line.startswith("#HMM - Sites"):
                type = "HMM - Sites"
            elif line.startswith("#HMM - Genes"):
                type = "HMM - Genes"
            elif line.startswith("#Resampling"):
                type = "Resampling"
            else:
                msg = """Please make sure the file has one of the following headers specifying the type of output as its first line:
    #Gumbel
    #HMM - Sites
    #HMM - Genes
    #Resampling
"""
                self.ShowError(msg)
                return
            
            data = {"path":path, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
            if self.verbose:
                print transit_prefix, "Adding File:", path
            wx.CallAfter(pub.sendMessage, "file", data=data)
        except:
            pass




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
        annotationPath = self.annotationFilePicker.GetPath()
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
        annotationPath = self.annotationFilePicker.GetPath()
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
        annotationPath = self.annotationFilePicker.GetPath()
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
        annotationpath = self.annotationFilePicker.GetPath()
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
 
        annotationpath = self.annotationFilePicker.GetPath()
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

        annotationpath = self.annotationFilePicker.GetPath()
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

        annotationpath = self.annotationFilePicker.GetPath()
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

                 

                
    def RunGumbelFunc(self, event):

        #Get options
        next = self.list_ctrl.GetNextSelected(-1)
        all_selected = self.ctrlSelected()
        if len(all_selected) ==0:
            self.ShowError("Error: No dataset selected.")
            return

        annotationPath = self.annotationFilePicker.GetPath()
        if not annotationPath:
            self.ShowError("Error: No annotation file selected.")
            return
        
        
        pathCol = self.list_ctrl.GetColumnCount() - 1
        readPath = self.list_ctrl.GetItem(next, pathCol).GetText()
        readPathList = all_selected
        name = ntpath.basename(readPath)
        min_read = int(self.gumbelReadChoice.GetString(self.gumbelReadChoice.GetCurrentSelection()))
        samples = int(self.gumbelSampleText.GetValue())
        burnin = int(self.gumbelBurninText.GetValue())
        trim = int(self.gumbelTrimText.GetValue())
        repchoice = self.gumbelRepChoice.GetString(self.gumbelRepChoice.GetCurrentSelection())
        ignoreCodon = True
        ignoreNTerm = float(self.globalNTerminusText.GetValue())
        ignoreCTerm = float(self.globalCTerminusText.GetValue())


        #Get Default file name
        defaultFile = "gumbel_%s_s%d_b%d_t%d.dat" % (".".join(name.split(".")[:-1]), samples, burnin, trim)

        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()


        #Ask user for output:        
        outputPath = self.SaveFile(defaultDir, defaultFile)
        
        if outputPath: 
            output = open(outputPath, "w")
        else:
            return

        self.statusBar.SetStatusText("Running Gumbel Method")
        self.gumbel_count = 0
        self.gumbelProgress.SetRange(samples+burnin)
        self.gumbelButton.Disable()

        kwargs = {}
        kwargs["readPathList"] = readPathList
        kwargs["annotationPath"] = annotationPath
        kwargs["min_read"] = min_read
        kwargs["samples"] = samples
        kwargs["burnin"] = burnin
        kwargs["trim"] = trim
        kwargs["repchoice"] = repchoice
        kwargs["ignoreCodon"] = ignoreCodon
        kwargs["ignoreNTerm"] = ignoreNTerm
        kwargs["ignoreCTerm"] = ignoreCTerm
        kwargs["output"] = output

        thread = threading.Thread(target=gumbelMH.runGumbel, args=(wx, pub.sendMessage), kwargs=kwargs)
        thread.setDaemon(True)
        thread.start()


    def RunHMMFunc(self, event):
        next = self.list_ctrl.GetNextSelected(-1)
        annotationPath = self.annotationFilePicker.GetPath()
        all_selected = self.ctrlSelected()
        if len(all_selected) == 0:
            self.ShowError("Error: No dataset selected.")
            return
        if not annotationPath:
            self.ShowError("Error: No annotation file selected.")
            return

        pathCol = self.list_ctrl.GetColumnCount() - 1
        readPath = self.list_ctrl.GetItem(next, pathCol).GetText()
        readPathList = all_selected
        name = ntpath.basename(readPath)
        repchoice = self.hmmRepChoice.GetString(self.hmmRepChoice.GetCurrentSelection())
        ignoreCodon = True
        ignoreNTerm = float(self.globalNTerminusText.GetValue())
        ignoreCTerm = float(self.globalCTerminusText.GetValue())


        #Get Default file name
        defaultFile = "hmm_%s_sites.dat" % (".".join(name.split(".")[:-1]))

        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        #Ask user for output:
        outputPath = self.SaveFile(defaultDir, defaultFile)

        if outputPath:
            output = open(outputPath, "w")
        else:
            return

        self.statusBar.SetStatusText("Running HMM Method")
        self.hmm_count = 0
        T = len([1 for line in open(readPath).readlines() if not line.startswith("#")])
        self.hmmProgress.SetRange(T*4 +1)

        self.hmmButton.Disable()

        kwargs = {}
        kwargs["readPathList"] = readPathList
        kwargs["annotationPath"] = annotationPath
        kwargs["repchoice"] = repchoice
        kwargs["ignoreCodon"] = ignoreCodon
        kwargs["ignoreNTerm"] = ignoreNTerm
        kwargs["ignoreCTerm"] = ignoreCTerm
        kwargs["output"] = output

        #HMMThread(readPathList, annotationPath, repchoice, ignoreCodon, ignoreNTerm, ignoreCTerm, output)
        thread = threading.Thread(target=hmm_geom.runHMM, args=(wx, pub.sendMessage), kwargs=kwargs)
        thread.setDaemon(True)
        thread.start()



    def RunResamplingFunc(self, event):
        selected_ctrl = self.ctrlSelected()
        selected_exp = self.expSelected()

        if len(selected_ctrl) == 0 and len(selected_exp) == 0:
            self.ShowError("Error: No datasets selected.")
            return

        elif len(selected_ctrl) == 0:
            self.ShowError("Error: No control datasets selected.")
            return

        elif len(selected_exp) == 0:
            self.ShowError("Error: No experimental datasets selected.")
            return
        annotationPath = self.annotationFilePicker.GetPath()
        if not annotationPath:
            self.ShowError("Error: No annotation file selected.")
            return


        #Get Default file name
        defaultFile = "resampling_results.dat"

        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()

        #Ask user for output:        
        outputPath = self.SaveFile(defaultDir, defaultFile)

        if outputPath:
            output = open(outputPath, "w")
        else:
            return

        #Check if user wants individual histograms
        if self.resamplingHistCheck.GetValue():
            histPath = os.path.join(ntpath.dirname(outputPath), transit_tools.fetch_name(outputPath))
            if not os.path.isdir(histPath):
                os.makedirs(histPath)
        else:
            histPath = ""

        doAdaptive= self.resamplingAdaptiveCheck.GetValue()

        ignoreCodon = True
        ignoreNTerm = float(self.globalNTerminusText.GetValue())
        ignoreCTerm = float(self.globalCTerminusText.GetValue())

        ctrlString = ",".join(selected_ctrl)
        expString = ",".join(selected_exp)

        if self.verbose:
            print transit_prefix, "Control String:", ctrlString
            print transit_prefix, "Experim String:", expString
            print transit_prefix, "outputPath:", outputPath
            print transit_prefix, "histPath:", histPath
         

        sampleSize = int(self.resamplingSampleText.GetValue())
 
        self.statusBar.SetStatusText("Running Resampling Method")
        self.resampling_count = 0
        T = len([1 for line in open(annotationPath).readlines() if not line.startswith("geneID")])
        self.resamplingProgress.SetRange(T+1)

        self.resamplingButton.Disable()

        kwargs = {}
        kwargs["ctrlList"] = ctrlString.split(",")
        kwargs["expList"] = expString.split(",")
        kwargs["annotationPath"] = annotationPath
        kwargs["sampleSize"] = sampleSize
        kwargs["histPath"] = histPath
        kwargs["doAdaptive"] = doAdaptive
        kwargs["ignoreCodon"] = ignoreCodon
        kwargs["ignoreNTerm"] = ignoreNTerm
        kwargs["ignoreCTerm"] = ignoreCTerm
        kwargs["doNormalize"] = True
        kwargs["output"] = output

        ##############
        #Instantiating pyplot outside thread to avoid 'main thread not in event loop error'
        n, bins, patches = plt.hist([0,1,2,3,5], normed=1, facecolor='c', alpha=0.75, bins=100)
        plt.xlabel('Delta Sum')
        plt.clf()
        #############


        thread = threading.Thread(target=resampling.runResampling, args=(wx, pub.sendMessage), kwargs=kwargs)
        thread.setDaemon(True)
        thread.start()

        #ResamplingThread(ctrlString, expString, annotationPath, sampleSize, histPath, doAdaptive, ignoreCodon, ignoreNTerm, ignoreCTerm, output)


if __name__ == "__main__":


    #If no arguments, show GUI:
    if len(sys.argv) == 1:
        #refer manual for details
        app = wx.App(False)

        #create an object of CalcFrame
        frame = TnSeekFrame(None)
        #show the frame
        frame.Show(True)
        #start the applications
        app.MainLoop()

    else:
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(help='Methods', dest='command')
        
        # A gumbel command
        gumbel_parser = subparsers.add_parser('gumbel', help='Gumbel method')
        gumbel_parser.add_argument("annotation", help="Path to the annotation file in .prot_table format.")
        gumbel_parser.add_argument("control_files", help="Comma separated list of paths for replicate CONTROL files.")
        gumbel_parser.add_argument("output_file", help="Output filename.")
        gumbel_parser.add_argument("-m", "--minread", help="Smallest read-count considered to be an insertion. Default is 1.", default=1, type=int)
        gumbel_parser.add_argument("-s", "--samples", help="Number of samples performed by the MH sampler. Default is 10000.", default=10000, type=int)
        gumbel_parser.add_argument("-b", "--burnin", help="Burn in period, Skips this number of samples before getting estimates. See documentation. Default is 500.", default=500, type=int)
        gumbel_parser.add_argument("-t", "--trim", help="Number of samples to trim. See documentation. Default is 1.", default=1, type=int)
        gumbel_parser.add_argument("-r", "--rep", help="How to handle replicates: 'Sum' or 'Mean'. Default is Sum.", default="Sum", type=str)
        gumbel_parser.add_argument("-iN", "--ignoreN", help="Ignore TAs occuring at X%% of the N terminus. Default is 0%.", default=0.0, type=float)
        gumbel_parser.add_argument("-iC", "--ignoreC", help="Ignore TAs occuring at X%% of the C terminus. Default is 0%.", default=0.0, type=float)
        
        # A hmm command
        hmm_parser = subparsers.add_parser('hmm', help='HMM method')
        hmm_parser.add_argument("annotation", help="Path to the annotation file in .prot_table format.")
        hmm_parser.add_argument("control_files", help="Comma separated list of paths for replicate CONTROL files.")
        hmm_parser.add_argument("output_file", help="Output filename.")
        hmm_parser.add_argument("-r", "--rep", help="How to handle replicates: 'Sum' or 'Mean'. Default is Sum.", default="Sum")
        hmm_parser.add_argument("-iN", "--ignoreN", help="Ignore TAs occuring at X%% of the N terminus. Default is 0%.", default=0.0, type=float)
        hmm_parser.add_argument("-iC", "--ignoreC", help="Ignore TAs occuring at X%% of the C terminus. Default is 0%.", default=0.0, type=float)
        
        # A resampling command
        resampling_parser = subparsers.add_parser('resampling', help='Resampling method')
        resampling_parser.add_argument("annotation", help="Path to the annotation file in .prot_table format.")
        resampling_parser.add_argument("control_files", help="Comma separated list of paths for replicate CONTROL files.")
        resampling_parser.add_argument("exp_files", help="Comma separated list of paths for replicate EXPERIMENTAL files.")
        resampling_parser.add_argument("output_file", help="Output filename.")
        resampling_parser.add_argument("-s", "--samples", help="Number of permutation samples obtained for each gene. Default is 10000.", default=10000, type=int)
        resampling_parser.add_argument("-H", "--hist", help="Number of samples to trim. See documentation. Default is to NOT create histograms.", action='store_true')
        resampling_parser.add_argument("-a", "--adaptive", help="Adaptive resampling; faster at the risk of lower accuracy. Default is off (i.e. regular resampling).", action='store_true')
        resampling_parser.add_argument("-dN", "--dontNormalize", help="Set this flag to NOT Normalize between conditions. Default is to normalize.", action='store_true')
        resampling_parser.add_argument("-iN", "--ignoreN", help="Ignore TAs occuring at X%% of the N terminus. Default is 0%.", default=0.0, type=float)
        resampling_parser.add_argument("-iC", "--ignoreC", help="Ignore TAs occuring at X%% of the C terminus. Default is 0%.", default=0.0, type=float)
        
        args = parser.parse_args()

        if args.command == "gumbel":
            #prepare kwargs
            kwargs = {}
            kwargs["readPathList"] = args.control_files.split(",")
            kwargs["annotationPath"] = args.annotation
            kwargs["min_read"] = args.minread
            kwargs["samples"] = args.samples
            kwargs["burnin"] = args.burnin
            kwargs["trim"] = args.trim
            kwargs["repchoice"] = args.rep
            kwargs["ignoreCodon"] = True
            kwargs["ignoreNTerm"] = args.ignoreN
            kwargs["ignoreCTerm"] = args.ignoreC
            kwargs["output"] = open(args.output_file, "w")
    
            #thread = threading.Thread(target=gumbelMH.runGumbel, args=(None, None), kwargs=kwargs)
            #thread.start()
            gumbelMH.runGumbel(None, None, **kwargs)
            
        elif args.command == "hmm":
            kwargs = {}
            kwargs["readPathList"] = args.control_files.split(",")
            kwargs["annotationPath"] = args.annotation
            kwargs["repchoice"] = args.rep
            kwargs["ignoreCodon"] = True
            kwargs["ignoreNTerm"] = args.ignoreN
            kwargs["ignoreCTerm"] = args.ignoreC
            kwargs["output"] = open(args.output_file, "w")

            #thread = threading.Thread(target=hmm_geom.runHMM, args=(None, None), kwargs=kwargs)
            #thread.start()
            hmm_geom.runHMM(None, None, **kwargs)

        elif args.command == "resampling":
            kwargs = {}
            kwargs["ctrlList"] = args.control_files.split(",")
            kwargs["expList"] = args.exp_files.split(",")
            kwargs["annotationPath"] = args.annotation
            kwargs["sampleSize"] = args.samples
           
            if args.hist:
                histPath = os.path.join(ntpath.dirname(args.output_file), transit_tools.fetch_name(args.output_file))
                if not os.path.isdir(histPath):
                    os.makedirs(histPath)
            else:
                histPath = ""
            kwargs["histPath"] = histPath

            kwargs["doAdaptive"] = args.adaptive
            kwargs["ignoreCodon"] = True
            kwargs["ignoreNTerm"] = args.ignoreN
            kwargs["ignoreCTerm"] = args.ignoreC
            kwargs["doNormalize"] = True
            kwargs["output"] = open(args.output_file, "w")
            
            #thread = threading.Thread(target=resampling.runResampling, args=(None, None), kwargs=kwargs)
            #thread.start()
            resampling.runResampling(None, None, **kwargs)

        else:
            print "Error: Command not recognized!"



