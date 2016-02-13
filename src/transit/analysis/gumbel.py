import base
import wx


method_name = "Gumbel"

def Hide(wxobj):
    wxobj.gumbelPanel.Hide()

def Show(wxobj):
    wxobj.gumbelPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run Gumbel" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.gumbelPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    wxobj.gumbelPanel.SetMinSize( wx.Size( 50,1 ) )
    wxobj.gumbelPanel.SetMaxSize( wx.Size( 250,-1 ) )

    gumbelSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelLabel.Wrap( -1 )
    gumbelSection.Add( wxobj.gumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    bSizer14 = wx.BoxSizer( wx.HORIZONTAL )

    bSizer15 = wx.BoxSizer( wx.HORIZONTAL )

    bSizer16 = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelSampleLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelSampleLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelBurninLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Burn-In", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelBurninLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelBurninLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelTrimLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Trim", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelTrimLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelTrimLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelReadLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelReadLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelReadLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelRepLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelRepLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelRepLabel, 1, wx.ALL, 5 )

    bSizer15.Add( bSizer16, 1, wx.EXPAND, 5 )

    bSizer17 = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelSampleText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    wxobj.gumbelBurninText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"500", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelBurninText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    wxobj.gumbelTrimText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"1", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelTrimText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    gumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
    wxobj.gumbelReadChoice = wx.Choice( wxobj.gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelReadChoiceChoices, 0 )
    wxobj.gumbelReadChoice.SetSelection( 0 )
    bSizer17.Add( wxobj.gumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    gumbelRepChoiceChoices = [ u"Sum", u"Mean" ]
    wxobj.gumbelRepChoice = wx.Choice( wxobj.gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelRepChoiceChoices, 0 )
    wxobj.gumbelRepChoice.SetSelection( 0 )
    bSizer17.Add( wxobj.gumbelRepChoice, 0, wx.ALL|wx.EXPAND, 5 )


    bSizer15.Add( bSizer17, 1, wx.EXPAND, 5 )


    bSizer14.Add( bSizer15, 1, wx.EXPAND, 5 )
        
    gumbelSection.Add( bSizer14, 1, wx.EXPAND, 5 )

    wxobj.gumbelButton = wx.Button( wxobj.gumbelPanel, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
    gumbelSection.Add( wxobj.gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.progressLabel1 = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.progressLabel1.Wrap( -1 )
    gumbelSection.Add( wxobj.progressLabel1, 0, wx.ALL, 5 )

    wxobj.gumbelProgress = wx.Gauge( wxobj.gumbelPanel, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
    wxobj.gumbelProgress.SetValue( 0 )
    gumbelSection.Add( wxobj.gumbelProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


    wxobj.gumbelPanel.SetSizer( gumbelSection )
    wxobj.gumbelPanel.Layout()
    gumbelSection.Fit( wxobj.gumbelPanel )


    #Connect events
    wxobj.gumbelButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.gumbelPanel


def updateProgressBar(wxobj, count):
    wxobj.gumbelProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.gumbelProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.gumbelButton.Enable()





########## CLASS #######################
class Gumbel(base.SingleConditionMethod):
    """   
    Gumbel
 
    """
    def __init__(self,
                ctrldata,
                outpath,
                samples=10000,
                burnin=500,
                trim=1,
                minread=1,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0):

        base.SingleConditionMethod.__init__(self, short_name, long_name, description, ctrldata, replicates="Sum", normalization=None, LOESS=False, NTerminus=0.0, CTerminus=0.0)
        self.samples = samples
        self.burnin = burnin
        self.trim = trim
        self.minread = minread
        self.short_name = "Gumbel"
        self.long_name = "Gumbel Method"
        self.description = "Gumbel method from DeJesus et al. (Bioinformatics, 2013)"



    @classmethod
    def fromGUI(self, wxobj):
        #Get options
        next = wxobj.list_ctrl.GetNextSelected(-1)
        all_selected = wxobj.ctrlSelected()
        if len(all_selected) ==0:
            wxobj.ShowError("Error: No dataset selected.")
            return

        annotationPath = wxobj.annotationFilePicker.GetPath()
        if not annotationPath:
            wxobj.ShowError("Error: No annotation file selected.")
            return


        pathCol = wxobj.list_ctrl.GetColumnCount() - 1
        readPath = wxobj.list_ctrl.GetItem(next, pathCol).GetText()
        readPathList = all_selected
        name = transit_tools.basename(readPath)
        minread = int(wxobj.gumbelReadChoice.GetString(wxobj.gumbelReadChoice.GetCurrentSelection()))
        samples = int(wxobj.gumbelSampleText.GetValue())
        burnin = int(wxobj.gumbelBurninText.GetValue())
        trim = int(wxobj.gumbelTrimText.GetValue())
        replicates = wxobj.gumbelRepChoice.GetString(wxobj.gumbelRepChoice.GetCurrentSelection())
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        normalization=None,
        LOESS=False

        return Gumbel(ctrldata,
                outpath,
                samples,
                burnin,
                trim,
                minread,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus):


        #Get Default file name
        defaultFile = "gumbel_%s_s%d_b%d_t%d.dat" % (".".join(name.split(".")[:-1]), samples, burnin, trim)

        #Get Default directory
        #defaultDir = os.path.dirname(os.path.realpath(__file__))
        defaultDir = os.getcwd()


        #Ask user for output:        
        outputPath = wxobj.SaveFile(defaultDir, defaultFile)

        if outputPath:
            output = open(outputPath, "w")
        else:
            return

        wxobj.statusBar.SetStatusText("Running Gumbel Method")
        wxobj.gumbel_count = 0
        wxobj.gumbelProgress.SetRange(samples+burnin)
        wxobj.gumbelButton.Disable()



    @classmethod
    def fromconsole(self, kwargs):
        pass



    def Run(self):
        print "FLORF: Run"
        pass





if __name__ == "__main__":


    G = Gumbel("Gumbel", "Gumbel", "Gumbel Bayesian method", [])


    G.message("Printing the member variables:")   
    G.print_members()
    print ""
    G.message("Instructions:")   
    print G.getInstructions()



