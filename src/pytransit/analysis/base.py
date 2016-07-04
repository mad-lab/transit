#__all__ = []
import sys

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
    
import traceback
import datetime

import pytransit
import pytransit.tnseq_tools as tnseq_tools
import pytransit.transit_tools as transit_tools

file_prefix = "[FileDisplay]"

class TransitFile:
    #TODO write docstring

    def __init__(self, identifier="#Unknown", colnames=[]):
        #TODO write docstring
        self.identifier = identifier
        self.colnames = colnames

    def getData(self, path, colnames):
        #TODO write docstring
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

    def getHeader(self, path):
        #TODO write docstring
        return "Generic Transit File Type."


    def getMenus(self):
        menus = [("Display in Track View", self.displayInTrackView)]
        return menus

    def displayInTrackView(self, displayFrame, event):

        #print "Self:", self
        #print "Frame:", displayFrame
        #print "Event:", event
        #print "Frame parent:", displayFrame.parent
        try:
            gene = displayFrame.grid.GetCellValue(displayFrame.row, 0)
            displayFrame.parent.allViewFunc(displayFrame, gene)
        except Exception as e:
            print file_prefix, "Error occurred: %s" % e 

class AnalysisGUI:
    
    def __init__(self):
        self.wxobj = None
        self.panel = None

    def Hide(self):
        self.panel.Hide()

    def Show(self):
        self.panel.Show()

    def Enable(self):
        self.panel.Enable()



    def definePanel(self, wxobj):
        #TODO: write docstring
        
        self.wxobj = wxobj
        wPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        Section = wx.BoxSizer( wx.VERTICAL )

        Label = wx.StaticText(wPanel, id=wx.ID_ANY, label=str("Options"), pos=wx.DefaultPosition, size=wx.DefaultSize, style=0 )
        Label.Wrap( -1 )
        Section.Add( Label, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        Sizer1 = wx.BoxSizer( wx.HORIZONTAL )
        Section.Add( Sizer1, 1, wx.EXPAND, 5 )
        
        Button = wx.Button( wPanel, wx.ID_ANY, u"Run", wx.DefaultPosition, wx.DefaultSize, 0 )
        Section.Add( Button, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        wPanel.SetSizer( Section )
        wPanel.Layout()
        Section.Fit( wPanel )

        #Connect events
        Button.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )
        self.panel = wPanel



class AnalysisMethod:
    '''
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    '''
    
    def __init__(self, short_name, long_name, description, output, annotation_path, wxobj=None):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description
        self.output = output
        self.annotation_path = annotation_path

        self.newWx = newWx
        self.wxobj = wxobj


    @classmethod 
    def fromGUI(self, wxobj):
        #TODO: write docstring
        raise NotImplementedError

    @classmethod
    def fromargs(self, rawargs):
        #TODO: write docstring
        raise NotImplementedError

    @classmethod
    def fromconsole(self):
        #TODO: write docstring
        try:
            return self.fromargs(sys.argv[2:])
        except IndexError as e:
            traceback.print_exc()
            print self.usage_string()
        except TypeError as e:
            print "Error: %s" % str(e)
            traceback.print_exc()
            print self.usage_string()
        except ValueError as e:
            print "Error: %s" % str(e)
            traceback.print_exc()
            print self.usage_string()
        except Exception as e:
            print "Error: %s" % str(e)
            traceback.print_exc()
            print self.usage_string()
        sys.exit() 

    @classmethod
    def usage_string(self):
        #TODO: write docstring
        raise NotImplementedError

    def Run(self):
        #TODO write docstring
        raise NotImplementedError


    def print_members(self):
        #TODO: write docstring
        members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
        for m in members:
            print "%s = %s" % (m, getattr(self, m))


    def add_file(self, path=None, filetype=None):

        #TODO: write docstring
        if not path:
            path = self.output.name
        if not filetype:
            filetype = self.short_name

        data = {"path":path, "type":filetype, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}

        if self.wxobj:
            if newWx:
                wx.CallAfter(pub.sendMessage, "file", data=data)
            else:
                wx.CallAfter(pub.sendMessage, "file", data)

    def finish(self):
        #TODO: write docstring
        if self.wxobj:
            if newWx:
                wx.CallAfter(pub.sendMessage,"finish", msg=self.short_name.lower())
            else:
                wx.CallAfter(pub.sendMessage,"finish", self.short_name.lower())

    def progress_update(self, text, count):
        #TODO: write docstring
        if self.wxobj:
            if newWx:
                wx.CallAfter(pub.sendMessage, "progress", msg=(self.short_name, count))
            else:
                wx.CallAfter(pub.sendMessage, "progress", (self.short_name, count))
            wx.Yield()

    def progress_range(self, count):
        #TODO: write docstring
        if self.wxobj:
            if newWx:
                wx.CallAfter(pub.sendMessage, "progressrange", msg=count)
            else:
                wx.CallAfter(pub.sendMessage, "progressrange", count)
            wx.Yield()

        

    def status_message(self, text):
        #TODO: write docstring
        if self.wxobj:
            if newWx:
                wx.CallAfter(pub.sendMessage, "status", msg=(self.short_name, text))
            else:
                wx.CallAfter(pub.sendMessage, "status", (self.short_name, text))
            wx.Yield()

    def console_message(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

    def console_message_inplace(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text) )
        sys.stdout.flush()

    def transit_message(self, text):
        #TODO: write docstring
        self.console_message(text)
        self.status_message(text)

    def transit_message_inplace(self, text):
        #TODO: write docstring
        self.console_message_inplace(text)
        self.status_message(text)


    def transit_error(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowError(text)

    def transit_warning(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)    


class SingleConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine essentiality in a single condition (e.g. Gumbel, Binomial, HMM).
    '''

    def __init__(self, short_name, long_name, description, ctrldata, annotation_path, output, replicates="Sum", normalization=None, LOESS=False, ignoreCodon=True, NTerminus=0.0, CTerminus=0.0, wxobj=None):
        AnalysisMethod.__init__(self, short_name, long_name, description, output, annotation_path, wxobj)
        self.ctrldata = ctrldata
        self.replicates = replicates
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus



class DualConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine changes in essentiality between two conditions (e.g. Resampling, DEHMM).
    '''

    def __init__(self, short_name, long_name, description, ctrldata, expdata, annotation_path, output, normalization, replicates="Sum", LOESS=False, ignoreCodon=True, NTerminus=0.0, CTerminus=0.0, wxobj=None):
        AnalysisMethod.__init__(self, short_name, long_name, description, output, annotation_path, wxobj)
        self.ctrldata = ctrldata
        self.expdata = expdata
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus



class TransitAnalysis:
    def __init__(self, sn, ln, desc, tn, method_class=AnalysisMethod, gui_class=AnalysisGUI, filetypes=[TransitFile]):
        self.short_name = sn
        self.long_name = ln
        self.description = desc
        self.transposons = tn
        self.method = method_class
        self.gui = gui_class()
        self.filetypes = filetypes


    def __str__(self):
        return """Analysis Method:
    Short Name:  %s
    Long Name:   %s
    Description: %s
    Method:      %s
    GUI:         %s""" % (self.short_name, self.long_name, self.description, self.method, self.gui)


    def fullname(self):
        return "[%s]  -  %s" % (self.short_name, self.long_name)

    def getInstructionsText(self):
        return ""

    def getDescriptionText(self):
        return self.description 

    def getTransposonsText(self):
        if len(self.transposons) == 0:
            return "Tn attribute missing!"
        elif len(self.transposons) == 1:
            return "Intended for %s only" % self.transposons[0]
        elif len(self.transposons) == 2:
            return "Intended for %s && %s" % tuple(self.transposons)
        else:
            return "Intended for " + ", ".join(self.transposons[:-1]) + ", and " + self.transposons[-1]
    


if __name__ == "__main__":
    pass



