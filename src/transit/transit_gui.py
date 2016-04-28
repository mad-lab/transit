# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Jun  6 2014)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc
import transit.analysis.gumbel
from wx.lib.buttons import GenBitmapTextButton




###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
    
    def __init__( self, parent ):
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"TRANSIT", pos = wx.DefaultPosition, size = wx.Size( 1300,975 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
        

        #Define variables

        #Define Text
        self.instructions_text = """
1. Choose the annotation file ("prot table") that corresponds to the datasets to be analyzed.
2. Add the desired Control and Experimental datasets.
3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the "View" button.
4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.
"""


        # Define ART
        bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))



        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        
        bSizer1 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.mainWindow = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1,-1 ), wx.HSCROLL|wx.VSCROLL )
        self.mainWindow.SetScrollRate( 5, 5 )
        self.mainWindow.SetMinSize( wx.Size( 700,-1 ) )
        
        bSizer4 = wx.BoxSizer( wx.VERTICAL )

        # Organism
        orgSizer = wx.StaticBoxSizer( wx.StaticBox( self.mainWindow, wx.ID_ANY, u"Organism" ), wx.VERTICAL )
        
        bSizer10 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.m_staticText5 = wx.StaticText( self.mainWindow, wx.ID_ANY, u"Annotation File:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText5.Wrap( -1 )
        bSizer10.Add( self.m_staticText5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        

        self.annotationFilePicker = GenBitmapTextButton(self.mainWindow, 1, bmp, '[Click to add Annotation File (.prot_table)]', size= wx.Size(400, 30))


        bSizer10.Add( self.annotationFilePicker, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        
        self.m_panel2 = wx.Panel( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel2.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel2.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer10.Add( self.m_panel2, 1, wx.EXPAND |wx.ALL, 5 )
        
        orgSizer.Add( bSizer10, 1, wx.EXPAND, 5 )
        
        
        bSizer4.Add( orgSizer, 0, wx.EXPAND, 5 )
        
        ctrlSizer = wx.StaticBoxSizer( wx.StaticBox( self.mainWindow, wx.ID_ANY, u"Control Samples" ), wx.VERTICAL )
        
        bSizer2 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.ctrlRemoveButton = wx.Button( self.mainWindow, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer2.Add( self.ctrlRemoveButton, 0, wx.ALL, 5 )
        
        self.ctrlView = wx.Button( self.mainWindow, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.ctrlView.Hide()
        
        bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )
        
        self.ctrlScatter = wx.Button( self.mainWindow, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.ctrlScatter.Hide()
        
        bSizer2.Add( self.ctrlScatter, 0, wx.ALL, 5 )
        
        self.ctrlFilePicker = GenBitmapTextButton(self.mainWindow, 1, bmp, '[Click to add Control Dataset(s)]', size= wx.Size(400, 30))
        bSizer2.Add( self.ctrlFilePicker, 1, wx.ALL, 5 )
        
        self.m_panel21 = wx.Panel( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel21.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel21.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer2.Add( self.m_panel21, 1, wx.EXPAND |wx.ALL, 5 )
        
        
        ctrlSizer.Add( bSizer2, 0, wx.EXPAND, 5 )
        
        self.list_ctrl = wx.ListCtrl( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )

        self.list_ctrl.SetMaxSize(wx.Size(940,200))
        ctrlSizer.Add( self.list_ctrl, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( ctrlSizer, 1, wx.EXPAND, 5 )
        
        expSizer1 = wx.StaticBoxSizer( wx.StaticBox( self.mainWindow, wx.ID_ANY, u"Experimental Samples" ), wx.VERTICAL )
        
        bSizer3 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.expSizer = wx.Button( self.mainWindow, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer3.Add( self.expSizer, 0, wx.ALL, 5 )
        
        self.expView = wx.Button( self.mainWindow, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.expView.Hide()
        
        bSizer3.Add( self.expView, 0, wx.ALL, 5 )
        
        self.expScatter = wx.Button( self.mainWindow, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.expScatter.Hide()
        
        bSizer3.Add( self.expScatter, 0, wx.ALL, 5 )
       

        self.expFilePicker = GenBitmapTextButton(self.mainWindow, 1, bmp, '[Click to add Experimental Dataset(s)]', size= wx.Size(400, 30))
        bSizer3.Add( self.expFilePicker, 1, wx.ALL, 5 )
        
        self.m_panel22 = wx.Panel( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel22.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel22.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer3.Add( self.m_panel22, 1, wx.EXPAND |wx.ALL, 5 )
        
        
        expSizer1.Add( bSizer3, 0, wx.EXPAND, 5 )
        
        self.list_exp = wx.ListCtrl( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        self.list_exp.SetMaxSize(wx.Size(940, 200))
        expSizer1.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( expSizer1, 1, wx.EXPAND, 5 )
        
        filesSizer = wx.StaticBoxSizer( wx.StaticBox( self.mainWindow, wx.ID_ANY, u"Results Files" ), wx.VERTICAL )
        
        bSizer141 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.displayButton = wx.Button( self.mainWindow, wx.ID_ANY, u"Display Table", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer141.Add( self.displayButton, 0, wx.ALL, 5 )
        
        self.graphFileButton = wx.Button( self.mainWindow, wx.ID_ANY, u"Display Graph", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.graphFileButton.Hide()
        
        bSizer141.Add( self.graphFileButton, 0, wx.ALL, 5 )
        
        self.addFileButton = GenBitmapTextButton(self.mainWindow, 1 , bmp, 'Add Results File', size= wx.Size(150, 30))


        bSizer141.Add( self.addFileButton, 0, wx.ALL, 5 )
        
        graphFileChoiceChoices = [ u"[Visualization Options]", u"Volcano Plot", u"Histogram of Total Gene Counts" ]
        self.graphFileChoice = wx.Choice( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, graphFileChoiceChoices, 0 )
        self.graphFileChoice.SetSelection( 0 )
        bSizer141.Add( self.graphFileChoice, 0, wx.ALL, 5 )
        
        
        filesSizer.Add( bSizer141, 0, 0, 5 )
        
        self.list_files = wx.ListCtrl( self.mainWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        self.list_files.SetMaxSize(wx.Size(940,200))
        filesSizer.Add( self.list_files, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( filesSizer, 1, wx.EXPAND, 5 )
        
        
        self.mainWindow.SetSizer( bSizer4 )
        self.mainWindow.Layout()
        bSizer4.Fit( self.mainWindow )
        bSizer1.Add( self.mainWindow, 1, wx.ALL|wx.EXPAND, 5 )
        
        self.m_panel5 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.DOUBLE_BORDER|wx.TAB_TRAVERSAL )
        self.m_panel5.SetMaxSize( wx.Size( 2,-1 ) )
        
        bSizer1.Add( self.m_panel5, 0, wx.ALL|wx.EXPAND, 5 )
        
        self.optionsWindow = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1,-1 ), wx.HSCROLL|wx.VSCROLL |wx.EXPAND)
        self.optionsWindow.SetScrollRate( 5, 5 )
        self.optionsWindow.SetMinSize( wx.Size( 310,1000 ) )
        #self.optionsWindow.SetMaxSize( wx.Size( 310,1000 ) )
        #self.optionsWindow.BackgroundColour = (200, 0, 20) 
        

        self.optionsSizer = wx.BoxSizer( wx.VERTICAL )

        # Logo Section
        self.logoImg = wx.StaticBitmap( self.optionsWindow, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
        self.optionsSizer.Add( self.logoImg, 0, wx.ALL|wx.ALIGN_CENTER, 5 )

        self.versionLabel = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
        self.versionLabel.Wrap( -1 )
        self.versionLabel.SetFont( wx.Font( 10, 74, 90, 92, False, "Sans" ) )

        self.optionsSizer.Add( self.versionLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        # Method Information 
        self.methodInfoText = wx.StaticBox( self.optionsWindow, wx.ID_ANY, u"Instructions" )
        self.methodInfoSizer = wx.StaticBoxSizer( self.methodInfoText, wx.VERTICAL )

        self.methodShortText = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.methodShortText.Wrap(250)
        self.methodShortText.Hide()
        self.methodInfoSizer.Add( self.methodShortText, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        self.methodLongText = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.methodLongText.Wrap(250)
        self.methodLongText.Hide()
        self.methodInfoSizer.Add( self.methodLongText, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        self.methodDescText = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.methodDescText.Wrap(250)
        self.methodDescText.Hide()
        self.methodInfoSizer.Add( self.methodDescText, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.methodInstructions = wx.StaticText( self.optionsWindow, wx.ID_ANY, self.instructions_text, wx.DefaultPosition, wx.DefaultSize, 0 )
        self.methodInstructions.Wrap( 250 )
        self.methodInfoSizer.Add( self.methodInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        self.optionsSizer.Add( self.methodInfoSizer, 0, wx.ALL|wx.EXPAND, 5 )




        # Method Options 
        self.methodSizerText = wx.StaticBox( self.optionsWindow, wx.ID_ANY, u"Method Options" )
        self.methodSizer = wx.StaticBoxSizer( self.methodSizerText, wx.VERTICAL )

        #self.methodSizerText.Hide()
        #self.methodSizer.SetMinSize( wx.Size( 250,-1 ) ) 
        
        self.m_panel1 = wx.Panel( self.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel1.SetMinSize( wx.Size( 50,1 ) )
        
        self.methodSizer.Add( self.m_panel1, 0, wx.ALL, 5 )
        
        self.globalLabel = wx.StaticText( self.optionsWindow, wx.ID_ANY, u"Global Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.globalLabel.Wrap( -1 )
        self.methodSizer.Add( self.globalLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        self.globalPanel = wx.Panel( self.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.globalPanel.SetMinSize( wx.Size( 50,90 ) )
        self.globalPanel.SetMaxSize( wx.Size( 250,-1) )
       
        #self.globalPanel.BackgroundColour = (200, 230, 250) 
        bSizer1431 = wx.BoxSizer( wx.VERTICAL )
        
        bSizer1521 = wx.BoxSizer( wx.HORIZONTAL )
        
        bSizer1621 = wx.BoxSizer( wx.VERTICAL )
        
        self.globalNTerminusLabel = wx.StaticText( self.globalPanel, wx.ID_ANY, u"Ignore N-Terminus %:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.globalNTerminusLabel.Wrap( -1 )
        bSizer1621.Add( self.globalNTerminusLabel, 0, wx.ALL, 5 )
        
        self.globalCTerminusLabel = wx.StaticText( self.globalPanel, wx.ID_ANY, u"Ignore C-Terminus %:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.globalCTerminusLabel.Wrap( -1 )
        bSizer1621.Add( self.globalCTerminusLabel, 0, wx.ALL, 5 )
        
        
        bSizer1521.Add( bSizer1621, 1, wx.EXPAND, 5 )
        
        bSizer1721 = wx.BoxSizer( wx.VERTICAL )
        
        self.globalNTerminusText = wx.TextCtrl( self.globalPanel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer1721.Add( self.globalNTerminusText, 0, wx.ALL, 5 )
        
        self.globalCTerminusText = wx.TextCtrl( self.globalPanel, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer1721.Add( self.globalCTerminusText, 0, wx.ALL, 5 )
        
        
        bSizer1521.Add( bSizer1721, 1, wx.EXPAND, 5 )
        
        
        bSizer1431.Add( bSizer1521, 1, wx.EXPAND, 5 )
        
        
        self.globalPanel.SetSizer( bSizer1431 )
        self.globalPanel.Layout()
        bSizer1431.Fit( self.globalPanel )
        self.methodSizer.Add( self.globalPanel, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        #--------------------#

        self.optionsSizer.Add( self.methodSizer, 0, wx.EXPAND, 5 )

        self.optionsWindow.SetSizer( self.optionsSizer )
        self.optionsWindow.Layout()

        self.optionsWindow.Fit()
        
        #self.optionsSizer.Fit( self.optionsWindow )
        bSizer1.Add( self.optionsWindow, 0, wx.ALL, 5 )
        
        #--------------------#        

        self.SetSizer( bSizer1 )
        self.Layout()
        self.m_menubar1 = wx.MenuBar( 0 )
        self.fileMenuItem = wx.Menu()
        self.exportMenuItem = wx.Menu()
        self.ctrlExportMenuItem = wx.Menu()
        self.ctrlExportIGVMenuItem = wx.MenuItem( self.ctrlExportMenuItem, wx.ID_ANY, u"to IGV", wx.EmptyString, wx.ITEM_NORMAL )
        self.ctrlExportMenuItem.AppendItem( self.ctrlExportIGVMenuItem )
        
        self.exportMenuItem.AppendSubMenu( self.ctrlExportMenuItem, u"Control Samples" )
        
        self.expExportMenuItem = wx.Menu()
        self.expExportIGVMenuItem = wx.MenuItem( self.expExportMenuItem, wx.ID_ANY, u"to IGV", wx.EmptyString, wx.ITEM_NORMAL )
        self.expExportMenuItem.AppendItem( self.expExportIGVMenuItem )
        
        self.exportMenuItem.AppendSubMenu( self.expExportMenuItem, u"Experimental Samples" )
        
        self.allExportMenuItem = wx.Menu()
        self.allExportIGVMenuItem = wx.MenuItem( self.allExportMenuItem, wx.ID_ANY, u"to IGV", wx.EmptyString, wx.ITEM_NORMAL )
        self.allExportMenuItem.AppendItem( self.allExportIGVMenuItem )
        
        self.exportMenuItem.AppendSubMenu( self.allExportMenuItem, u"All Samples" )
        
        self.fileMenuItem.AppendSubMenu( self.exportMenuItem, u"Export" )
        
        self.convertMenuItem = wx.Menu()
        self.annotationConvertPTToPTTMenu = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, u"prot_table to PTT", wx.EmptyString, wx.ITEM_NORMAL )
        self.convertMenuItem.AppendItem( self.annotationConvertPTToPTTMenu )
        
        self.annotationConvertPTToGFF3Menu = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, u"prot_table to GFF3", wx.EmptyString, wx.ITEM_NORMAL )
        self.convertMenuItem.AppendItem( self.annotationConvertPTToGFF3Menu )
        
        self.annotationConvertPTTToPT = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, u"PTT to prot_table", wx.EmptyString, wx.ITEM_NORMAL )
        self.convertMenuItem.AppendItem( self.annotationConvertPTTToPT )
        
        self.annotationConvertGFF3ToPT = wx.MenuItem( self.convertMenuItem, wx.ID_ANY, u"GFF3 to prot_table", wx.EmptyString, wx.ITEM_NORMAL )
        self.convertMenuItem.AppendItem( self.annotationConvertGFF3ToPT )
        
        self.fileMenuItem.AppendSubMenu( self.convertMenuItem, u"Convert" )
        
        self.fileExitMenuItem = wx.MenuItem( self.fileMenuItem, wx.ID_ANY, u"Exit", wx.EmptyString, wx.ITEM_NORMAL )
        self.fileMenuItem.AppendItem( self.fileExitMenuItem )
        
        self.m_menubar1.Append( self.fileMenuItem, u"&File" ) 
        
        self.viewMenuItem = wx.Menu()
        self.scatterMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Scatter Plot", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.scatterMenuItem )
        
        self.trackMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Track View", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.trackMenuItem )
        
        self.qcMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Quality Control", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.qcMenuItem )
        
        self.m_menubar1.Append( self.viewMenuItem, u"&View" ) 
       

        # 
        self.methodsMenuItem = wx.Menu()
        self.himar1MenuItem = wx.Menu()
        self.tn5MenuItem = wx.Menu()
        
        
        self.methodsMenuItem.AppendSubMenu(self.himar1MenuItem, "Himar1 Methods")
        self.methodsMenuItem.AppendSubMenu(self.tn5MenuItem, "Tn5 Methods")
        self.m_menubar1.Append( self.methodsMenuItem, u"&Analysis" )
        
 
        self.SetMenuBar( self.m_menubar1 )

        
        self.helpMenuItem = wx.Menu()
        self.documentationMenuItem = wx.MenuItem(self.helpMenuItem, wx.ID_ANY, u"Documentation", wx.EmptyString, wx.ITEM_NORMAL)
        self.helpMenuItem.AppendItem(self.documentationMenuItem)
        self.aboutMenuItem = wx.MenuItem(self.helpMenuItem, wx.ID_ANY, u"About", wx.EmptyString, wx.ITEM_NORMAL)
        self.helpMenuItem.AppendItem(self.aboutMenuItem)
        self.m_menubar1.Append( self.helpMenuItem, u"&Help" )
        
        
        self.statusBar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
        



        self.Centre( wx.BOTH )
        
        # Connect Events

        self.annotationFilePicker.Bind( wx.EVT_BUTTON, self.annotationFileFunc )
        self.ctrlRemoveButton.Bind( wx.EVT_BUTTON, self.ctrlRemoveFunc )
        self.ctrlView.Bind( wx.EVT_BUTTON, self.allViewFunc )
        self.ctrlScatter.Bind( wx.EVT_BUTTON, self.scatterFunc )
        self.ctrlFilePicker.Bind( wx.EVT_BUTTON, self.loadCtrlFileFunc )
        self.expSizer.Bind( wx.EVT_BUTTON, self.expRemoveFunc )
        self.expView.Bind( wx.EVT_BUTTON, self.allViewFunc )
        self.expScatter.Bind( wx.EVT_BUTTON, self.scatterFunc )
        self.expFilePicker.Bind( wx.EVT_BUTTON, self.loadExpFileFunc )
        self.displayButton.Bind( wx.EVT_BUTTON, self.displayFileFunc )
        self.graphFileButton.Bind( wx.EVT_BUTTON, self.graphFileFunc )
        self.addFileButton.Bind( wx.EVT_BUTTON, self.addFileFunc )
        self.graphFileChoice.Bind( wx.EVT_CHOICE, self.graphFileFunc )
        self.Bind( wx.EVT_MENU, self.ctrlToIGV, id = self.ctrlExportIGVMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.expToIGV, id = self.expExportIGVMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.allToIGV, id = self.allExportIGVMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.annotationPT_to_PTT, id = self.annotationConvertPTToPTTMenu.GetId() )
        self.Bind( wx.EVT_MENU, self.annotationPT_to_GFF3, id = self.annotationConvertPTToGFF3Menu.GetId() )
        self.Bind( wx.EVT_MENU, self.annotationPTT_to_PT, id = self.annotationConvertPTTToPT.GetId() )
        self.Bind( wx.EVT_MENU, self.annotationGFF3_to_PT, id = self.annotationConvertGFF3ToPT.GetId() )
        self.Bind( wx.EVT_MENU, self.Exit, id = self.fileExitMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.scatterFunc, id = self.scatterMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.allViewFunc, id = self.trackMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.qcFunc, id = self.qcMenuItem.GetId() )

        self.Bind( wx.EVT_MENU, self.aboutFunc, id = self.aboutMenuItem.GetId() )
        self.Bind( wx.EVT_MENU, self.documentationFunc, id = self.documentationMenuItem.GetId() )
    
    def __del__( self ):
        pass
    
    
    # Virtual event handlers, overide them in your derived class

    def onHimar1Checked(self, event):
        event.Skip()

    def onTn5Checked(self, event):
        event.Skip()

    def annotationFileFunc( self, event ):
        event.Skip()
    
    def ctrlRemoveFunc( self, event ):
        event.Skip()
    
    def allViewFunc( self, event ):
        event.Skip()
    
    def scatterFunc( self, event ):
        event.Skip()
    
    def loadCtrlFileFunc( self, event ):
        event.Skip()
    
    def expRemoveFunc( self, event ):
        event.Skip()
    
    def loadExpFileFunc( self, event ):
        event.Skip()
    
    def displayFileFunc( self, event ):
        event.Skip()
    
    def graphFileFunc( self, event ):
        event.Skip()
    
    def addFileFunc( self, event ):
        event.Skip()
    
    
    def MethodSelectFunc( self, event ):
        event.Skip()
    
#    def RunGumbelFunc( self, event ):
#        event.Skip()
    
    def RunBinomialFunc( self, event ):
        event.Skip()
    
    def LoessPrevFunc( self, event ):
        event.Skip()
    
    def RunHMMFunc( self, event ):
        event.Skip()
    
    
    def RunResamplingFunc( self, event ):
        event.Skip()
    
    
    def RunDEHMMFunc( self, event ):
        event.Skip()
    
    def ctrlToIGV( self, event ):
        event.Skip()
    
    def expToIGV( self, event ):
        event.Skip()
    
    def allToIGV( self, event ):
        event.Skip()
    
    def annotationPT_to_PTT( self, event ):
        event.Skip()
    
    def annotationPT_to_GFF3( self, event ):
        event.Skip()
    
    def annotationPTT_to_PT( self, event ):
        event.Skip()
    
    def annotationGFF3_to_PT( self, event ):
        event.Skip()
    
    def Exit( self, event ):
        event.Skip()
    
    
    def qcFunc( self, event ):
        event.Skip()
    
    def aboutFunc( self, event ):
        event.Skip()

    def documentationFunc( self, event ):
        event.Skip()
    

