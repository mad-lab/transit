#__all__ = []
import sys
import wx
import datetime
#Check if wx is the newest 3.0+ version:
try:
    from wx.lib.pubsub import pub
    pub.subscribe
    newWx = True
except AttributeError as e:
    from wx.lib.pubsub import Publisher as pub
    newWx = False




class AnalysisGUI:
    
    def __init__(self, sn, ln, desc, wxobj):
        self.wxobj =  wxobj
        self.short_name = sn
        self.long_name = ln
        self.description = desc
        self.wxobj = wxobj
        self.panel = self.getPanel()
        self.wxobj.methodSizer.Add(self.panel, 1, wx.EXPAND |wx.ALL, 5 )


    def __str__(self):
        return """Method GUI:
    short name: %s
    long name: %s
    description: %s""" % (self.short_name, self.long_name, self.description)

    def fullname(self):
        return "[%s]  -  %s" % (self.short_name, self.long_name)

    def Hide(self):
        self.panel.Hide()

    def Show(self):
        self.panel.Show()

    def Enable(self):
        #self.panel.Button.Enable()
        self.panel.Enable()


    def getInstructions(self):
        return """ Description: 
    %s
                """ % (self.description)


    def getInstructions_old(self):
        return """Instruction:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



    def getPanel(self):
        #TODO: write docstring
        #raise NotImplementedError
        wPanel = wx.Panel( self.wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        Section = wx.BoxSizer( wx.VERTICAL )

        #print self.wxobj
        #print wx.ID_ANY
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


        return wPanel




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


    def __str__(self):
        #TODO: write docstring
        return "%s (%s): %s" % (self.long_name, self.short_name, self.description)

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
            print self.usage_string()
        except TypeError as e:
            print "Error: %s" % str(e)
            print self.usage_string()
        except ValueError as e:
            print "Error: %s" % str(e)
            print self.usage_string()
        except Exception as e:
            print "Error: %s" % str(e)
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


    def add_file(self, path=None):
        #TODO: write docstring
        if not path:
            data = {"path":self.output.name, "type":self.short_name, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        else:
            data = {"path":path, "type":self.short_name, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}

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





