#__all__ = []
import sys

try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0

if hasWx:
    import wx.xrc
    from wx.lib.buttons import GenBitmapTextButton
    from pubsub import pub
    import wx.adv

import traceback
import datetime
import numpy
import pytransit.transit_tools as transit_tools

file_prefix = "[FileDisplay]"


class InvalidArgumentException(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        super(InvalidArgumentException, self).__init__(message)

#

if hasWx:
    class InfoIcon(wx.StaticBitmap):
        def __init__(self, panel, flag, bmp=None, tooltip=""):
            if not bmp:
                bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (16, 16))
            wx.StaticBitmap.__init__(self, panel, flag, bmp)
            tp = wx.ToolTip(tooltip)
            self.SetToolTip(tp)


class TransitGUIBase:

    def __init__(self):
        self.wxobj = None
        self.short_name = "TRANSIT"
        self.long_name = "TRANSIT"
        self.short_desc = "TRANSIT - Short Description"
        self.long_desc =  "TRANSIT - Long Description"

#

    def status_message(self, text, time=-1):
        #TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "status", msg=(self.short_name, text, time))
            wx.Yield()

#

    def console_message(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

#

    def console_message_inplace(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text) )
        sys.stdout.flush()

#

    def transit_message(self, text):
        #TODO: write docstring
        self.console_message(text)
        self.status_message(text)

#

    def transit_message_inplace(self, text):
        #TODO: write docstring
        self.console_message_inplace(text)
        self.status_message(text)

#

    def transit_error(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowError(text)

#

    def transit_warning(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)

#

class TransitFile (TransitGUIBase):
    #TODO write docstring

#

    def __init__(self, identifier="#Unknown", colnames=[]):
        #TODO write docstring
        TransitGUIBase.__init__(self)
        self.identifier = identifier
        self.colnames = colnames

#

    def getData(self, path, colnames):
        #TODO write docstring
        row = 0
        data = []
        shownError = False
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.split("\t")
            tmp[-1] = tmp[-1].strip()
            #print colnames
            #print  len(colnames), len(tmp)
            try:
                rowdict = dict([(colnames[i], tmp[i]) for i in range(len(colnames))])
            except Exception as e:
                if not shownError:
                    self.transit_warning("Error reading data! This may be caused by trying to load a old results file, when the format has changed.")
                    shownError = True
                rowdict = dict([(colnames[i], tmp[i]) for i in range(min(len(colnames), len(tmp)))])
            data.append((row, rowdict))
            row+=1
        return data

#

    def getHeader(self, path):
        #TODO write docstring
        return "Generic Transit File Type."

#

    def getMenus(self):
        menus = [("Display in Track View", self.displayInTrackView)]
        return menus

#

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

#

class AnalysisGUI:

    def __init__(self):
        self.wxobj = None
        self.panel = None
        self.LABELSIZE = (100,-1)
        self.WIDGETSIZE = (100,-1)

#

    def Hide(self):
        self.panel.Hide()

#

    def Show(self):
        self.panel.Show()

#

    def Enable(self):
        self.panel.Enable()

#

    def GlobalEnable(self):
        pass

#

    def GlobalHide(self):
        pass

#

    def GlobalShow(self):
        pass

#

    def GlobalDisable(self):
        pass

#

    def definePanel(self, wxobj):
        #TODO: write docstring

        self.wxobj = wxobj
        wPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        Section = wx.BoxSizer( wx.VERTICAL )

        Label = wx.StaticText(wPanel, id=wx.ID_ANY, label=str("Method Options"), pos=wx.DefaultPosition, size=(130, -1), style=0 )
        Label.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
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

#

    def defineTextBox(self, panel, labelText="", widgetText="", tooltipText="", labSize=None, widgetSize=None):
        if not labSize: labSize = self.LABELSIZE
        if not widgetSize: widgetSize = self.WIDGETSIZE

        sizer = wx.BoxSizer( wx.HORIZONTAL )
        label = wx.StaticText(panel, wx.ID_ANY, labelText, wx.DefaultPosition, labSize, 0)
        label.Wrap( -1 )
        textBox = wx.TextCtrl( panel, wx.ID_ANY, widgetText, wx.DefaultPosition, widgetSize, 0 )
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        sizer.Add(textBox, 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        sizer.Add(InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText), 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        return (label, textBox, sizer)

#

    def defineChoiceBox(self, panel, labelText="", widgetChoice=[""], tooltipText="", labSize=None, widgetSize=None):
        if not labSize: labSize = self.LABELSIZE
        if not widgetSize: widgetSize = self.WIDGETSIZE

        sizer = wx.BoxSizer( wx.HORIZONTAL )
        label = wx.StaticText(panel, wx.ID_ANY, labelText, wx.DefaultPosition, labSize, 0)
        label.Wrap( -1 )
        choiceBox = wx.Choice( panel, wx.ID_ANY, wx.DefaultPosition, widgetSize, widgetChoice, 0 )
        choiceBox.SetSelection(0)
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        sizer.Add(choiceBox, 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        sizer.Add(InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText), 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        return (label, choiceBox, sizer)

#

    def defineCheckBox(self, panel, labelText="", widgetCheck=False, tooltipText="", widgetSize=None):
        if not widgetSize: widgetSize = self.WIDGETSIZE
        sizer = wx.BoxSizer( wx.HORIZONTAL )
        checkBox = wx.CheckBox(panel, label = labelText, size=widgetSize)
        checkBox.SetValue(widgetCheck)
        sizer.Add(checkBox, 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        sizer.Add(InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText), 0, wx.ALIGN_CENTER_VERTICAL, 5 )
        return (checkBox, sizer)

#

class AnalysisMethod:
    '''
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    '''

    def __init__(self, short_name, long_name, short_desc, long_desc, output, annotation_path, wxobj=None):
        self.short_name = short_name
        self.long_name = long_name
        self.short_desc = short_desc
        self.long_desc = long_desc
        self.output = output
        self.annotation_path = annotation_path

        self.WX_VERSION = WX_VERSION
        self.wxobj = wxobj

#

    @classmethod
    def fromGUI(self, wxobj):
        #TODO: write docstring
        raise NotImplementedError

#

    @classmethod
    def fromargs(self, rawargs):
        #TODO: write docstring
        raise NotImplementedError

#

    @classmethod
    def fromconsole(self):
        #TODO: write docstring
        try:
            return self.fromargs(sys.argv[2:])
        except InvalidArgumentException as e:
            print "Error: %s" % str(e)
            print self.usage_string()
        except IndexError as e:
            print "Error: %s" % str(e)
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

#

    @classmethod
    def usage_string(self):
        #TODO: write docstring
        raise NotImplementedError

#

    def Run(self):
        #TODO write docstring
        raise NotImplementedError

#

    def print_members(self):
        #TODO: write docstring
        members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
        for m in members:
            print "%s = %s" % (m, getattr(self, m))

#

    def add_file(self, path=None, filetype=None):

        #TODO: write docstring
        if not path:
            path = self.output.name
        if not filetype:
            filetype = self.short_name

        data = {"path":path, "type":filetype, "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}

        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "file", data=data)

#

    def finish(self):
        #TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage,"finish", msg=self.short_name.lower())

#

    def progress_update(self, text, count):
        #TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "progress", msg=(self.short_name, count))
            wx.Yield()

        self.transit_message_inplace(text)
#

    def progress_range(self, count):
        #TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "progressrange", msg=count)
            wx.Yield()

#

    def status_message(self, text, time=-1):
        #TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "status", msg=(self.short_name, text, time))
            wx.Yield()

#

    def console_message(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

#

    def console_message_inplace(self, text):
        #TODO: write docstring
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text) )
        sys.stdout.flush()

#

    def transit_message(self, text):
        #TODO: write docstring
        self.console_message(text)
        self.status_message(text)

#

    def transit_message_inplace(self, text):
        #TODO: write docstring
        self.console_message_inplace(text)
        self.status_message(text)

#

    def transit_error(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowError(text)

#

    def transit_warning(self,text):
        self.transit_message(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)


class SingleConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine essentiality in a single condition (e.g. Gumbel, Binomial, HMM).
    '''

    def __init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output, replicates="Sum", normalization=None, LOESS=False, ignoreCodon=True, NTerminus=0.0, CTerminus=0.0, wxobj=None):
        AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, output,
            annotation_path, wxobj)
        self.ctrldata = ctrldata
        self.replicates = replicates
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus

#

class DualConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine changes in essentiality between two conditions (e.g. Resampling, DEHMM).
    '''

    def __init__(self, short_name, long_name, short_desc, long_desc, ctrldata, expdata, annotation_path, output, normalization, replicates="Sum", LOESS=False, ignoreCodon=True, NTerminus=0.0, CTerminus=0.0, wxobj=None):
        AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc,
            output, annotation_path, wxobj)
        self.ctrldata = ctrldata
        self.expdata = expdata
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus

#

class QuadConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine changes in essentiality between four conditions (e.g. GI).
    '''

    def __init__(self, short_name, long_name, short_desc, long_desc, ctrldataA, ctrldataB, expdataA, expdataB, annotation_path, output, normalization, replicates="Sum", LOESS=False, ignoreCodon=True, NTerminus=0.0, CTerminus=0.0, wxobj=None):
        AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc,
            output, annotation_path, wxobj)
        self.ctrldataA = ctrldataA
        self.ctrldataB = ctrldataB
        self.expdataA = expdataA
        self.expdataB = expdataB
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus
#

class MultiConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that compare essentiality between multiple conditions (e.g Anova).
    '''

    def __init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation_path, output, normalization=None, LOESS=False, ignoreCodon=True, wxobj=None, ignored_conditions=[], included_conditions=[], nterm=0.0, cterm=0.0):
        AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, output,
                annotation_path, wxobj)
        self.combined_wig = combined_wig
        self.metadata = metadata
        self.normalization = normalization
        self.NTerminus = nterm
        self.CTerminus = cterm
        self.unknown_cond_flag = "FLAG-UNMAPPED-CONDITION-IN-WIG"
        self.ignored_conditions = ignored_conditions
        self.included_conditions = included_conditions

    def filter_wigs_by_conditions(self, data, conditions, covariates = [], interactions = [], ignored_conditions = [], included_conditions = []):
        """
            Filters conditions that are ignored/included.
            ([[Wigdata]], [Condition], [[Covar]], [Condition], [Condition]) -> Tuple([[Wigdata]], [Condition])
        """
        ignored_conditions, included_conditions = (set(ignored_conditions), set(included_conditions))
        d_filtered, cond_filtered, filtered_indexes = [], [], [];
        if len(ignored_conditions) > 0 and len(included_conditions) > 0:
            self.transit_error("Both ignored and included conditions have len > 0", ignored_conditions, included_conditions)
            sys.exit(0)
        elif (len(ignored_conditions) > 0):
            self.transit_message("conditions ignored: {0}".format(ignored_conditions))
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag) and (c not in ignored_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)
        elif (len(included_conditions) > 0):
            self.transit_message("conditions included: {0}".format(included_conditions))
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag) and (c in included_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)
        else:
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)

        covariates_filtered = [[c[i] for i in filtered_indexes] for c in covariates]
        interactions_filtered = [[c[i] for i in filtered_indexes] for c in interactions]

        return (numpy.array(d_filtered),
                numpy.array(cond_filtered),
                numpy.array(covariates_filtered),
                numpy.array(interactions_filtered))

#

class TransitAnalysis:
    def __init__(self, sn, ln, short_desc, long_desc, tn, method_class=AnalysisMethod, gui_class=AnalysisGUI, filetypes=[TransitFile]):
        self.short_name = sn
        self.long_name = ln
        self.short_desc = short_desc
        self.long_desc = long_desc
        self.transposons = tn
        self.method = method_class
        self.gui = gui_class()
        self.filetypes = filetypes

#

    def __str__(self):
        return """Analysis Method:
    Short Name:  %s
    Long Name:   %s
    Short Desc:  %s
    Long Desc:   %s
    Method:      %s
    GUI:         %s""" % (self.short_name, self.long_name, self.short_desc, self.long_desc, self.method, self.gui)

#

    def fullname(self):
        return "[%s]  -  %s" % (self.short_name, self.short_desc)

#

    def getInstructionsText(self):
        return ""

#

    def getDescriptionText(self):
        return self.long_desc

#

    def getTransposonsText(self):
        if len(self.transposons) == 0:
            return "Tn attribute missing!"
        elif len(self.transposons) == 1:
            return "Intended for %s only" % self.transposons[0]
        elif len(self.transposons) == 2:
            return "Intended for %s or %s" % tuple(self.transposons)
        else:
            return "Intended for " + ", ".join(self.transposons[:-1]) + ", and " + self.transposons[-1]


#

if __name__ == "__main__":
    pass



