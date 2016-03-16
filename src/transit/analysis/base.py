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


class AnalysisMethod:
    '''
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    '''
    
    def __init__(self,short_name, long_name, description, output, annotation_path, wxobj=None):
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
        return self.fromargs(sys.argv[1:])
    
    def Run(self):
        #TODO write docstring
        raise NotImplementedError


    def print_members(self):
        #TODO: write docstring
        members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
        for m in members:
            print "%s = %s" % (m, getattr(self, m))


    def add_file(self):
        #TODO: write docstring
        data = {"path":self.output.name, "type":self.short_name, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
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





