#importing wx files
import wx
import os
import ntpath

# trash view stuff
import trash
 
#import the newly created GUI file
import tn_seek_gui

 
#importing * : to enable writing sin(13) instead of math.sin(13)
from math import *


import gumbelMH
import threading



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






class GumbelThread(Thread):
    """Gumbel Worker Thread Class."""
 
    #----------------------------------------------------------------------
    def __init__(self):
        """Init Worker Thread Class."""
        Thread.__init__(self)
        self.start()    # start the thread
 
    #----------------------------------------------------------------------
    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread.
        for i in range(20):
            time.sleep(1)
            wx.CallAfter(pub.sendMessage, "update", msg="")









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


        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()


    def HideAllOptions(self):
        self.HideGumbelOptions()
        self.HideHMMOptions()
        self.HideResamplingOptions()

    def HideGumbelOptions(self):
        self.gumbelLabel.Hide()
        self.gumbelSampleLabel.Hide()
        self.gumbelSampleText.Hide()
        self.gumbelReadLabel.Hide()
        self.gumbelReadChoice.Hide()
        self.gumbelButton.Hide()

    def ShowGumbelOptions(self):
        self.gumbelLabel.Show()
        self.gumbelSampleLabel.Show()
        self.gumbelSampleText.Show()
        self.gumbelReadLabel.Show()
        self.gumbelReadChoice.Show()
        self.gumbelButton.Show()


    def HideHMMOptions(self):
        self.hmmLabel.Hide()
        self.hmmButton.Hide()
    
    def ShowHMMOptions(self):
        self.hmmLabel.Show()
        self.hmmButton.Show()


    def HideResamplingOptions(self):
        self.resamplingLabel.Hide()
        self.resamplingButton.Hide()
    

    def ShowResamplingOptions(self):
        self.resamplingLabel.Show()
        self.resamplingButton.Show()
        

 
    #wx calls this function with and 'event' object
    def solveFunc(self,event):
        try:
            #evaluate the string in 'text' and put the answer back
            ans = eval(self.text.GetValue())
            self.text.SetValue (str(ans))
        except Exception:
            print 'error'
    #put a blank string in text when 'Clear' is clicked
    def clearFunc(self,event):
        self.text.SetValue(str(''))


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
        #secondWindow = window2(None)
        #secondWindow.Show()
        annotationpath = self.annotationFilePicker.GetPath()
        if DEBUG: print "Annotation path selected:", annotationpath
        if DEBUG: print "Focus:", self.FindFocus()

        next = self.list_ctrl.GetNextSelected(-1)
        if next > -1 and annotationpath:
            dataset = self.list_ctrl.GetItem(next, 4).GetText()
            if DEBUG:
                print "Visualizing counts for:", self.list_ctrl.GetItem(next, 0).GetText()
            viewWindow = trash.TrashFrame(self, dataset, annotationpath)
            viewWindow.Show()
        elif next == -1:
            if DEBUG:
                print "No datasets selected to visualize!"
        else:
            if DEBUG:
                print "No annotation file selected"

    def annotationFileFunc(self, event):
        pass
        
        
    def MethodSelectFunc(self, event):
        X = self.methodChoice.GetCurrentSelection()


        if X == 0:
            self.HideAllOptions()
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


    def RunGumbelFunc(self, event):
        next = self.list_ctrl.GetNextSelected(-1)
        readPath = self.list_ctrl.GetItem(next, 0).GetText()
        annotationPath = self.annotationFilePicker.GetPath()
        min_read = int(self.gumbelReadChoice.GetString(self.gumbelReadChoice.GetCurrentSelection()))
        samples = int(self.gumbelSampleText.GetValue())
        if DEBUG:
            print "Running Gumbel Method for %d Samples" % (samples)

        #gumbelMH.runGumbel(readPath, annotationPath, min_read, samples)

        #t1 = threading.Thread(target=func1 args=(arg1, arg2))
        #t1.start()
        #t1.join()
        t1 = threading.Thread(target=gumbelMH.runGumbel, args=(readPath, annotationPath, min_read, samples))
        t1.start()
        t1.join()

        if DEBUG:
            print "Finished Running Gumbel Method"

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


