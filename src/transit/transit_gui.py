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
        
        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        
        bSizer1 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.m_scrolledWindow2 = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1,-1 ), wx.HSCROLL|wx.VSCROLL )
        self.m_scrolledWindow2.SetScrollRate( 5, 5 )
        self.m_scrolledWindow2.SetMinSize( wx.Size( 700,-1 ) )
        
        bSizer4 = wx.BoxSizer( wx.VERTICAL )
       

        #Method choice
        choiceSizer = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow2, wx.ID_ANY, u"Method" ), wx.VERTICAL )

        choiceSizer_H = wx.BoxSizer( wx.HORIZONTAL )


        self.methodCheckBoxHimar1 = wx.CheckBox(self.m_scrolledWindow2, label = 'Himar1',pos = (10,10))
        self.methodCheckBoxHimar1.SetValue(True)
        self.methodCheckBoxTn5 = wx.CheckBox(self.m_scrolledWindow2, label = 'Tn5',pos = (10,10))
        self.methodCheckBoxTn5.SetValue(True)
        choiceSizer_H.Add(self.methodCheckBoxHimar1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        choiceSizer_H.Add(self.methodCheckBoxTn5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        self.methodChoiceStaticText = wx.StaticText( self.m_scrolledWindow2, wx.ID_ANY, u"Method Choice:", wx.DefaultPosition, wx.DefaultSize, 0 ) 
        self.methodChoiceStaticText.Wrap( -1 )

        choiceSizer_H.Add(self.methodChoiceStaticText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )


        methodChoiceChoices = [ u"[Choose Method]"]
        self.methodChoice = wx.Choice( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, methodChoiceChoices, 0 )
        self.methodChoice.SetSelection( 0 )
        choiceSizer_H.Add( self.methodChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        
        
        choiceSizer.Add(choiceSizer_H, 1, wx.EXPAND, 5 )
        bSizer4.Add( choiceSizer, 0, wx.EXPAND, 5 )


        # Organism
        orgSizer = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow2, wx.ID_ANY, u"Organism" ), wx.VERTICAL )
        
        bSizer10 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.m_staticText5 = wx.StaticText( self.m_scrolledWindow2, wx.ID_ANY, u"Annotation File:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText5.Wrap( -1 )
        bSizer10.Add( self.m_staticText5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        
        self.annotationFilePicker = wx.FilePickerCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Prot Table (*.prot_table)|*.prot_table;|\nProt Table (*.txt)|*.txt;|\nProt Table (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
        bSizer10.Add( self.annotationFilePicker, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
        
        self.m_panel2 = wx.Panel( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel2.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel2.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer10.Add( self.m_panel2, 1, wx.EXPAND |wx.ALL, 5 )
        
        orgSizer.Add( bSizer10, 1, wx.EXPAND, 5 )
        
        
        bSizer4.Add( orgSizer, 0, wx.EXPAND, 5 )
        
        ctrlSizer = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow2, wx.ID_ANY, u"Control Samples" ), wx.VERTICAL )
        
        bSizer2 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.ctrlRemoveButton = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer2.Add( self.ctrlRemoveButton, 0, wx.ALL, 5 )
        
        self.ctrlView = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.ctrlView.Hide()
        
        bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )
        
        self.ctrlScatter = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.ctrlScatter.Hide()
        
        bSizer2.Add( self.ctrlScatter, 0, wx.ALL, 5 )
        
        bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))
        self.ctrlFilePicker = GenBitmapTextButton(self.m_scrolledWindow2, 1, bmp, '[Click to add Control Dataset(s)]', size= wx.Size(400, 30))
        #self.ctrlFilePicker = wx.FilePickerCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )

        #self.ctrlFilePicker = wx.FilePickerCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FD_MULTIPLE )
        #self.ctrlFilePicker = wx.FileDialog( self.m_scrolledWindow2, message=u"Select a file", wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", style=wx.FD_MULTIPLE )

        bSizer2.Add( self.ctrlFilePicker, 1, wx.ALL, 5 )
        
        self.m_panel21 = wx.Panel( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel21.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel21.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer2.Add( self.m_panel21, 1, wx.EXPAND |wx.ALL, 5 )
        
        
        ctrlSizer.Add( bSizer2, 0, wx.EXPAND, 5 )
        
        self.list_ctrl = wx.ListCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        ctrlSizer.Add( self.list_ctrl, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( ctrlSizer, 1, wx.EXPAND, 5 )
        
        expSizer1 = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow2, wx.ID_ANY, u"Experimental Samples" ), wx.VERTICAL )
        
        bSizer3 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.expSizer = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer3.Add( self.expSizer, 0, wx.ALL, 5 )
        
        self.expView = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.expView.Hide()
        
        bSizer3.Add( self.expView, 0, wx.ALL, 5 )
        
        self.expScatter = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.expScatter.Hide()
        
        bSizer3.Add( self.expScatter, 0, wx.ALL, 5 )
       

        bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))
        self.expFilePicker = GenBitmapTextButton(self.m_scrolledWindow2, 1, bmp, '[Click to add Experimental Dataset(s)]', size= wx.Size(400, 30))
        #self.expFilePicker = wx.FilePickerCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.EmptyString, u"Select a .wig file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
        bSizer3.Add( self.expFilePicker, 1, wx.ALL, 5 )
        
        self.m_panel22 = wx.Panel( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel22.SetMinSize( wx.Size( 100,-1 ) )
        self.m_panel22.SetMaxSize( wx.Size( 150,-1 ) )
        
        bSizer3.Add( self.m_panel22, 1, wx.EXPAND |wx.ALL, 5 )
        
        
        expSizer1.Add( bSizer3, 0, wx.EXPAND, 5 )
        
        self.list_exp = wx.ListCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        expSizer1.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( expSizer1, 1, wx.EXPAND, 5 )
        
        filesSizer = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow2, wx.ID_ANY, u"Results Files" ), wx.VERTICAL )
        
        bSizer141 = wx.BoxSizer( wx.HORIZONTAL )
        
        self.displayButton = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Display Table", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer141.Add( self.displayButton, 0, wx.ALL, 5 )
        
        self.graphFileButton = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Display Graph", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.graphFileButton.Hide()
        
        bSizer141.Add( self.graphFileButton, 0, wx.ALL, 5 )
        
        #self.addFileButton = wx.Button( self.m_scrolledWindow2, wx.ID_ANY, u"Add Results File", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.addFileButton = GenBitmapTextButton(self.m_scrolledWindow2, 1 , bmp, 'Add Results File', size= wx.Size(150, 30))


        bSizer141.Add( self.addFileButton, 0, wx.ALL, 5 )
        
        graphFileChoiceChoices = [ u"[Visualization Options]", u"Volcano Plot", u"Histogram of Total Gene Counts" ]
        self.graphFileChoice = wx.Choice( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, graphFileChoiceChoices, 0 )
        self.graphFileChoice.SetSelection( 0 )
        bSizer141.Add( self.graphFileChoice, 0, wx.ALL, 5 )
        
        
        filesSizer.Add( bSizer141, 0, 0, 5 )
        
        self.list_files = wx.ListCtrl( self.m_scrolledWindow2, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        filesSizer.Add( self.list_files, 1, wx.ALL|wx.EXPAND, 5 )
        
        
        bSizer4.Add( filesSizer, 1, wx.EXPAND, 5 )
        
        
        self.m_scrolledWindow2.SetSizer( bSizer4 )
        self.m_scrolledWindow2.Layout()
        bSizer4.Fit( self.m_scrolledWindow2 )
        bSizer1.Add( self.m_scrolledWindow2, 1, wx.ALL|wx.EXPAND, 5 )
        
        self.m_panel5 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.DOUBLE_BORDER|wx.SUNKEN_BORDER|wx.TAB_TRAVERSAL )
        self.m_panel5.SetMaxSize( wx.Size( 2,-1 ) )
        
        bSizer1.Add( self.m_panel5, 0, wx.ALL|wx.EXPAND, 5 )
        
        self.m_scrolledWindow1 = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1,-1 ), wx.HSCROLL|wx.VSCROLL )
        self.m_scrolledWindow1.SetScrollRate( 5, 5 )
        self.m_scrolledWindow1.SetMinSize( wx.Size( 310,1000 ) )
        
        self.methodSizer = wx.StaticBoxSizer( wx.StaticBox( self.m_scrolledWindow1, wx.ID_ANY, u"Methods" ), wx.VERTICAL )
        
        self.methodSizer.SetMinSize( wx.Size( 250,-1 ) ) 
        self.m_staticText13 = wx.StaticText( self.m_scrolledWindow1, wx.ID_ANY, u"Tn-Seq\nAnalysis", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText13.Wrap( -1 )
        self.m_staticText13.SetFont( wx.Font( 20, 74, 90, 92, False, "Sans" ) )
        self.m_staticText13.Hide()
        
        self.methodSizer.Add( self.m_staticText13, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        self.logoImg = wx.StaticBitmap( self.m_scrolledWindow1, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
        self.methodSizer.Add( self.logoImg, 0, wx.ALL, 5 )
        
        self.versionLabel = wx.StaticText( self.m_scrolledWindow1, wx.ID_ANY, u"v1.4.3", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
        self.versionLabel.Wrap( -1 )
        self.versionLabel.SetFont( wx.Font( 10, 74, 90, 92, False, "Sans" ) )
        
        self.methodSizer.Add( self.versionLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        #methodChoiceChoices = [ u"[Choose Method]"]
        #self.methodChoice = wx.Choice( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, methodChoiceChoices, 0 )
        #self.methodChoice.SetSelection( 0 )
        #methodSizer.Add( self.methodChoice, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
       
        self.mainInstructions = wx.StaticText( self.m_scrolledWindow1, wx.ID_ANY, u"Instructions:\n\n1. Choose the annotation file (\"prot table\") that corresponds to the datasets to be analyzed.\n2. Add the desired Control and Experimental datasets.\n3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the \"View\" button.\n4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.\n", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.mainInstructions.Wrap( 250 )
        self.methodSizer.Add( self.mainInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        self.m_panel1 = wx.Panel( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel1.SetMinSize( wx.Size( 50,1 ) )
        
        self.methodSizer.Add( self.m_panel1, 0, wx.ALL, 5 )
        
        self.globalLabel = wx.StaticText( self.m_scrolledWindow1, wx.ID_ANY, u"Global Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.globalLabel.Wrap( -1 )
        self.methodSizer.Add( self.globalLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
        
        self.globalPanel = wx.Panel( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.globalPanel.SetMinSize( wx.Size( 50,90 ) )
        self.globalPanel.SetMaxSize( wx.Size( 250,-1) )
        
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
        self.methodSizer.Add( self.globalPanel, 0, wx.ALL|wx.EXPAND, 5 )


        #print self.gui
        #gui = transit.analysis.defineGUI()
        #Add Methods
        #for method in transit.analysis.methods:
        #    methodSizer.Add( transit.analysis.methods[method]["module"].getPanel(self), 1, wx.EXPAND |wx.ALL, 5 )
        


        #-------------------#
        # Progress

        #self.progressPanel = wx.Panel( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        #self.progressSizer = wx.BoxSizer( wx.VERTICAL )

        #self.progressLabel = wx.StaticText( self.progressPanel, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.progressLabel.Wrap( -1 )
        #self.progressSizer.Add( self.progressLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        #self.progress = wx.Gauge( self.progressPanel, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
        #self.progressSizer.Add( self.progress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        #self.progressPanel.SetSizer( self.progressSizer )
        #self.progressPanel.Layout()
        #self.progressSizer.Fit( self.progressPanel )
        #self.methodSizer.Add( self.progressPanel, 1, wx.EXPAND |wx.ALL, 5 )

        #--------------------#

        self.m_scrolledWindow1.SetSizer( self.methodSizer )
        self.m_scrolledWindow1.Layout()
        self.methodSizer.Fit( self.m_scrolledWindow1 )
        bSizer1.Add( self.m_scrolledWindow1, 0, wx.ALL, 5 )
        
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
        
        self.m_menubar1.Append( self.fileMenuItem, u"File" ) 
        
        self.viewMenuItem = wx.Menu()
        self.scatterMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Scatter Plot", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.scatterMenuItem )
        
        self.trackMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Track View", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.trackMenuItem )
        
        self.qcMenuItem = wx.MenuItem( self.viewMenuItem, wx.ID_ANY, u"Quality Control", wx.EmptyString, wx.ITEM_NORMAL )
        self.viewMenuItem.AppendItem( self.qcMenuItem )
        
        self.m_menubar1.Append( self.viewMenuItem, u"View" ) 
        
        self.SetMenuBar( self.m_menubar1 )
        
        self.statusBar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
        
        self.Centre( wx.BOTH )
        
        # Connect Events

        self.methodCheckBoxHimar1.Bind(wx.EVT_CHECKBOX,self.onHimar1Checked)
        self.methodCheckBoxTn5.Bind(wx.EVT_CHECKBOX,self.onTn5Checked)
        self.annotationFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.annotationFileFunc )
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
        self.methodChoice.Bind( wx.EVT_CHOICE, self.MethodSelectFunc )
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
    

