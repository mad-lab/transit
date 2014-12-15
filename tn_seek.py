#importing wx files
import wx
import sys
import os
import time
import datetime
import ntpath
import threading
import multiprocessing as mp
from math import *
from wx.lib.pubsub import pub

# trash view stuff
import trash
#import the newly created GUI file
import tn_seek_gui

import fileDisplay
#
import gumbelMH
import hmm_geom
import hmm_tools
import resampling

wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"


DEBUG = True

def get_wig_stats(path):
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        if line.startswith("variable"): continue
        if line.startswith("location"): continue
        tmp = line.split()
        pos = int(tmp[0])
        rd = int(tmp[1])
        data.append(rd)
    return (sum([1 for rd in data if rd > 0])/float(len(data)),  sum(data)/float(len(data)), max(data)) 


def fetch_name(filepath):
    return os.path.splitext(ntpath.basename(filepath))[0]




class GumbelThread(threading.Thread):
    """Gumbel Worker Thread Class."""
 
    #----------------------------------------------------------------------
    def __init__(self, readPath, annotationPath, min_read, samples, burnin, trim, output):
        """Init Worker Thread Class."""

        #parameters
        self.readPath = readPath
        self.annotationPath = annotationPath
        self.min_read = min_read
        self.samples = samples
        self.burnin = burnin
        self.trim = trim
        self.output = output

        #thread
        threading.Thread.__init__(self)
        self.start()    # start the thread
 
    #----------------------------------------------------------------------
    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread.

        print "Running Gumbel Method for %d Samples" % (self.samples)
        gumbelMH.runGumbel(self.readPath, self.annotationPath, self.min_read, self.samples, self.burnin, self.trim, self.output, wx, pub.sendMessage)
        print "Finished Gumbel Method"
        if not self.output.name.startswith("<"):
            data = {"path":self.output.name, "type":"Gumbel", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
            print "Adding File:", self.output.name
            wx.CallAfter(pub.sendMessage, "file", data=data)

        wx.CallAfter(pub.sendMessage, "gumbel", msg="Finished!")
        wx.CallAfter(pub.sendMessage,"finish", msg="gumbel")



class HMMThread(threading.Thread):
    """HMM Worker Thread Class."""

    #----------------------------------------------------------------------
    def __init__(self, readPath, annotationPath, output):
        """Init Worker Thread Class."""

        #parameters
        self.readPath = readPath
        self.annotationPath = annotationPath
        self.output = output

        #thread
        threading.Thread.__init__(self)
        self.start()    # start the thread

    #----------------------------------------------------------------------
    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread.

        print "Running HMM Method"

        hmm_geom.runHMM(self.readPath, self.annotationPath, self.output, wx, pub.sendMessage)
        print "Finished HMM Method"
        if not self.output.name.startswith("<"):
            data = {"path":self.output.name, "type":"HMM - Sites", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
            
            print "Adding File:", self.output.name
            wx.CallAfter(pub.sendMessage, "file", data=data)
            
            wx.CallAfter(pub.sendMessage, "hmm", msg="Creating HMM Sites output...")
            genes_path = ".".join(self.output.name.split(".")[:-1]) + "_genes." + self.output.name.split(".")[-1]
            hmm_tools.post_process_genes(self.output.name, self.annotationPath, output=open(genes_path,"w"))
            data["path"] =  genes_path
            data["type"] = "HMM - Genes"
            print "Adding File:", genes_path
            wx.CallAfter(pub.sendMessage, "file", data=data)
            wx.CallAfter(pub.sendMessage, "hmm", msg="Finished!")
            wx.CallAfter(pub.sendMessage,"finish", msg="hmm")



class ResamplingThread(threading.Thread):
    """HMM Worker Thread Class."""

    #----------------------------------------------------------------------
    def __init__(self, ctrlString, expString, annotationPath, output):
        """Init Worker Thread Class."""

        #parameters
        self.ctrlString = ctrlString
        self.expString = expString
        self.annotationPath = annotationPath
        self.output = output

        #thread
        threading.Thread.__init__(self)
        self.start()    # start the thread

    #----------------------------------------------------------------------
    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread.

        print "Running Resampling Method"
        resampling.runResampling(self.ctrlString, self.expString, self.annotationPath, self.output, wx, pub.sendMessage)
        print "Finished Resampling Method"
        if not self.output.name.startswith("<"):
            data = {"path":self.output.name, "type":"Resampling", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
            print "Adding File:", self.output.name
            wx.CallAfter(pub.sendMessage, "file", data=data)

        wx.CallAfter(pub.sendMessage, "resampling", msg="Finished!")
        wx.CallAfter(pub.sendMessage,"finish", msg="resampling")




#inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TnSeekFrame(tn_seek_gui.MainFrame):
    #constructor
    def __init__(self,parent):
        #initialize parent class
        tn_seek_gui.MainFrame.__init__(self,parent)


        self.index_ctrl = 0
        self.list_ctrl.InsertColumn(0, 'File', width=210)
        self.list_ctrl.InsertColumn(1, 'Density', width=85)
        self.list_ctrl.InsertColumn(2, 'Mean Read', width=85)
        self.list_ctrl.InsertColumn(3, 'Max Read', width=85)
        self.list_ctrl.InsertColumn(4, 'Full Path', width=403)


        self.index_exp = 0
        self.list_exp.InsertColumn(0, 'File', width=210)
        self.list_exp.InsertColumn(1, 'Density', width=85)
        self.list_exp.InsertColumn(2, 'Mean Read', width=85)
        self.list_exp.InsertColumn(3, 'Max Read', width=85)
        self.list_exp.InsertColumn(4, 'Full Path',width=403)


        self.index_file = 0
        self.list_files.InsertColumn(0, 'Name', width=350)
        self.list_files.InsertColumn(1, 'Type', width=100)
        self.list_files.InsertColumn(2, 'Date', width=230)
        self.list_files.InsertColumn(3, 'Full Path',width=403)



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
        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()



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
        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()


    def HideGumbelOptions(self):
        self.gumbelLabel.Hide()
        self.gumbelSampleLabel.Hide()
        self.gumbelSampleText.Hide()
        self.gumbelBurninLabel.Hide()
        self.gumbelBurninText.Hide()
        self.gumbelTrimLabel.Hide()
        self.gumbelTrimText.Hide()
        self.gumbelReadLabel.Hide()
        self.gumbelReadChoice.Hide()
        self.gumbelButton.Hide()
        self.gumbelProgress.Hide()


    def ShowGumbelOptions(self):
        self.gumbelLabel.Show()
        self.gumbelSampleLabel.Show()
        self.gumbelSampleText.Show()
        self.gumbelBurninLabel.Show()
        self.gumbelBurninText.Show()
        self.gumbelTrimLabel.Show()
        self.gumbelTrimText.Show()
        self.gumbelReadLabel.Show()
        self.gumbelReadChoice.Show()
        self.gumbelButton.Show()
        self.gumbelProgress.Show()


    def HideHMMOptions(self):
        self.hmmLabel.Hide()
        self.hmmButton.Hide()
        self.hmmProgress.Hide()

    
    def ShowHMMOptions(self):
        self.hmmLabel.Show()
        self.hmmButton.Show()
        self.hmmProgress.Show()


    def HideResamplingOptions(self):
        self.resamplingLabel.Hide()
        self.resamplingButton.Hide()
        self.resamplingProgress.Hide()
    

    def ShowResamplingOptions(self):
        self.resamplingLabel.Show()
        self.resamplingButton.Show()
        self.resamplingProgress.Show()
 

    def SaveFile(self, DIR=".", FILE="", WC=""):
        """
        Create and show the Save FileDialog
        """
        path = ""
        dlg = wx.FileDialog(
            self, message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE, wildcard=WC, style=wx.SAVE|wx.OVERWRITE_PROMPT
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print "You chose the following output filename: %s" % path
        dlg.Destroy()
        return path

        
     

    def ctrlSelected(self, col=4):
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


    def expSelected(self, col=4):
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



    def loadCtrlFileFunc(self, event):
        try:
            fullpath = self.ctrlFilePicker.GetPath()
            name = ntpath.basename(fullpath)
            density,meanrd,maxrd = get_wig_stats(fullpath)
            self.list_ctrl.InsertStringItem(self.index_ctrl, name)
            self.list_ctrl.SetStringItem(self.index_ctrl, 1, "%2.1f" % (density*100))
            self.list_ctrl.SetStringItem(self.index_ctrl, 2, "%1.1f" % (meanrd))
            self.list_ctrl.SetStringItem(self.index_ctrl, 3, "%d" % (maxrd))
            self.list_ctrl.SetStringItem(self.index_ctrl, 4, "%s" % (fullpath))
            if DEBUG:
                print "Adding control item (%d): %s" % (self.index_ctrl, name)
            self.index_ctrl+=1
            #self.ctrlFilePicker.SetPath("test.dat")
            #self.ctrlFilePicker.SetLabel("test.dat")
        except e:
            print "Error:", e

        


    def loadExpFileFunc(self, event):
        try:
            fullpath = self.expFilePicker.GetPath()
            name = ntpath.basename(fullpath)
            density,meanrd,maxrd = get_wig_stats(fullpath)
            self.list_exp.InsertStringItem(self.index_exp, name)
            self.list_exp.SetStringItem(self.index_exp, 1, "%2.1f" % (density*100))
            self.list_exp.SetStringItem(self.index_exp, 2, "%1.1f" % (meanrd))
            self.list_exp.SetStringItem(self.index_exp, 3, "%d" % (maxrd))
            self.list_exp.SetStringItem(self.index_exp, 4, "%s" % (fullpath))
            if DEBUG:
                print "Adding experimental item (%d): %s" % (self.index_exp, name)
            self.index_exp+=1
        

        except e:
            print "Error:", e 


    def ctrlRemoveFunc(self, event):
        next = self.list_ctrl.GetNextSelected(-1)
        while next != -1:
            if DEBUG:
                print "Removing control item (%d): %s" % (next, self.list_ctrl.GetItem(next, 0).GetText())
            self.list_ctrl.DeleteItem(next)
            next = self.list_ctrl.GetNextSelected(-1)
            self.index_ctrl-=1 

 
    def expRemoveFunc(self, event):
        next = self.list_exp.GetNextSelected(-1)
        while next != -1:
            if DEBUG:
                print "Removing experimental item (%d): %s" % (next, self.list_exp.GetItem(next, 0).GetText())
            self.list_exp.DeleteItem(next)
            next = self.list_exp.GetNextSelected(-1)
            self.index_exp-=1



    def ctrlViewFunc(self, event):
        annotationpath = self.annotationFilePicker.GetPath()
        datasets = self.ctrlSelected()
        if datasets and annotationpath:
            if DEBUG:
                print "Visualizing counts for:", ", ".join(datasets)
            viewWindow = trash.TrashFrame(self, datasets, annotationpath)
            viewWindow.Show()
        elif not datasets:
            if DEBUG:
                print "No datasets selected to visualize!"
        else:
            if DEBUG:
                print "No annotation file selected"



    def expViewFunc(self, event):
        annotationpath = self.annotationFilePicker.GetPath()
        datasets = self.expSelected()
        if datasets and annotationpath:
            if DEBUG:
                print "Visualizing counts for:", ", ".join(datasets)
            viewWindow = trash.TrashFrame(self, datasets, annotationpath)
            viewWindow.Show()
        elif not datasets:
            if DEBUG:
                print "No datasets selected to visualize!"
        else:
            if DEBUG:
                print "No annotation file selected"





    def annotationFileFunc(self, event):
        print "Annotation File Selected:", self.annotationFilePicker.GetPath()
        
        
    def MethodSelectFunc(self, event):
        X = self.methodChoice.GetCurrentSelection()


        self.progressLabel.Show()
        if X == 0:
            self.HideAllOptions()
            self.progressLabel.Hide()
        elif X == 1:
            self.ShowGumbelOptions()
            self.HideHMMOptions()
            self.HideResamplingOptions()
        elif X == 2:
            self.ShowHMMOptions()
            self.HideGumbelOptions()
            self.HideResamplingOptions()
        elif X == 3:
            self.ShowResamplingOptions()
            self.HideGumbelOptions()
            self.HideHMMOptions()

        self.Layout()

        if DEBUG:
            print "Selected Method (%d): %s" % (X, self.methodChoice.GetString(X))


    def displayFileFunc(self, event):
        next = self.list_files.GetNextSelected(-1)
        if next > -1:
            dataset = self.list_files.GetItem(next, 3).GetText()
            if DEBUG:
                print "Displaying results:", self.list_files.GetItem(next, 0).GetText()

            fileWindow = fileDisplay.FileFrame(self, dataset, self.list_files.GetItem(next, 1).GetText())
            fileWindow.Show()
        else:
            if DEBUG:
                print "No results selected to display!"



    def tempFileFunc(self, event):
        file1= "gumbel_H37Rv_Sassetti_glycerol_s100_b5_t1.dat"
        type="Gumbel"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "hmm_H37Rv_Sassetti_glycerol_sites.dat"
        type="HMM - Sites"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "hmm_H37Rv_Sassetti_glycerol_genes.dat"
        type="HMM - Genes"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)


        file1= "resampling_results.dat"
        type="Resampling"
        data = {"path":file1, "type":type, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print "Adding File:", file1
        wx.CallAfter(pub.sendMessage, "file", data=data)




    def RunGumbelFunc(self, event):

        #Get options
        next = self.list_ctrl.GetNextSelected(-1)
        readPath = self.list_ctrl.GetItem(next, 4).GetText()
        name = ntpath.basename(readPath)
        annotationPath = self.annotationFilePicker.GetPath()
        min_read = int(self.gumbelReadChoice.GetString(self.gumbelReadChoice.GetCurrentSelection()))
        samples = int(self.gumbelSampleText.GetValue())
        burnin = int(self.gumbelBurninText.GetValue())
        trim = int(self.gumbelTrimText.GetValue())
        

        #Get Default file name
        defaultFile = "gumbel_%s_s%d_b%d_t%d.dat" % (".".join(name.split(".")[:-1]), samples, burnin, trim)

        #Get Default directory
        defaultDir = os.path.dirname(os.path.realpath(__file__))

        #Ask user for output:        
        outputPath = self.SaveFile(defaultDir, defaultFile)
        
        if outputPath: 
            output = open(outputPath, "w")
        else:
            #output = sys.stdout
            return

        self.statusBar.SetStatusText("Running Gumbel Method")
        self.gumbel_count = 0
        self.gumbelProgress.SetRange(samples+burnin)

        self.gumbelButton.Disable()
        X = GumbelThread(readPath, annotationPath, min_read, samples, burnin, trim, output)



    def RunHMMFunc(self, event):
        next = self.list_ctrl.GetNextSelected(-1)
        readPath = self.list_ctrl.GetItem(next, 4).GetText()
        name = ntpath.basename(readPath)
        annotationPath = self.annotationFilePicker.GetPath()

        #Get Default file name
        defaultFile = "hmm_%s_sites.dat" % (".".join(name.split(".")[:-1]))

        #Get Default directory
        defaultDir = os.path.dirname(os.path.realpath(__file__))

        #Ask user for output:
        outputPath = self.SaveFile(defaultDir, defaultFile)

        if outputPath:
            output = open(outputPath, "w")
        else:
            output = sys.stdout

        self.statusBar.SetStatusText("Running HMM Method")
        self.hmm_count = 0
        T = len([1 for line in open(readPath).readlines() if not line.startswith("#")])
        self.hmmProgress.SetRange(T+T-1+1)

        self.hmmButton.Disable()
        HMMThread(readPath, annotationPath, output)


    def RunResamplingFunc(self, event):
        selected_ctrl = []
        current = -1
        while True:
            next = self.list_ctrl.GetNextSelected(current)
            if next == -1:
                break
            path = self.list_ctrl.GetItem(next, 4).GetText()
            selected_ctrl.append(path)
            current = next
    
        selected_exp = []
        current = -1
        while True:
            next = self.list_exp.GetNextSelected(current)
            if next == -1:
                break
            path = self.list_exp.GetItem(next, 4).GetText()
            selected_exp.append(path)
            current = next

        annotationPath = self.annotationFilePicker.GetPath()

        #Get Default file name
        defaultFile = "resampling_results.dat"

        #Get Default directory
        defaultDir = os.path.dirname(os.path.realpath(__file__))

        #Ask user for output:        
        outputPath = self.SaveFile(defaultDir, defaultFile)

        if outputPath:
            output = open(outputPath, "w")
        else:
            output = sys.stdout
 
        ctrlString = ",".join(selected_ctrl)
        expString = ",".join(selected_exp)

        print "Control String:", ctrlString
        print "Experim String:", expString
        
        self.statusBar.SetStatusText("Running Resampling Method")
        self.resampling_count = 0
        T = len([1 for line in open(annotationPath).readlines() if not line.startswith("geneID")])
        self.resamplingProgress.SetRange(T+1)

        self.resamplingButton.Disable()
        ResamplingThread(ctrlString, expString, annotationPath, output)


if __name__ == "__main__":
    #mandatory in wx, create an app, False stands for not deteriction stdin/stdout
    #refer manual for details
    app = wx.App(False)
     
    #create an object of CalcFrame
    frame = TnSeekFrame(None)
    #show the frame
    frame.Show(True)
    #start the applications
    app.MainLoop()


