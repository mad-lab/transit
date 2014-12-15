# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Jun  6 2014)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc

###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Tn-Seq Analysis", pos = wx.DefaultPosition, size = wx.Size( 1092,852 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer1 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer4 = wx.BoxSizer( wx.VERTICAL )
		
		orgSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Organism" ), wx.VERTICAL )
		
		bSizer10 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.m_staticText5 = wx.StaticText( self, wx.ID_ANY, u"Annotation File:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText5.Wrap( -1 )
		bSizer10.Add( self.m_staticText5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.annotationFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Prot Table (*.prot_table)|*.prot_table;|\nProt Table (*.txt)|*.txt;|\nProt Table (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		bSizer10.Add( self.annotationFilePicker, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		
		orgSizer.Add( bSizer10, 1, wx.EXPAND, 5 )
		
		
		bSizer4.Add( orgSizer, 0, wx.EXPAND, 5 )
		
		ctrlSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Control Samples" ), wx.VERTICAL )
		
		bSizer2 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.ctrlRemoveButton = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.ctrlRemoveButton, 0, wx.ALL, 5 )
		
		self.ctrlView = wx.Button( self, wx.ID_ANY, u"View", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )
		
		self.ctrlFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		bSizer2.Add( self.ctrlFilePicker, 1, wx.ALL, 5 )
		
		
		ctrlSizer.Add( bSizer2, 0, wx.EXPAND, 5 )
		
		self.list_ctrl = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		ctrlSizer.Add( self.list_ctrl, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( ctrlSizer, 1, wx.EXPAND, 5 )
		
		expSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Experimental Samples" ), wx.VERTICAL )
		
		bSizer3 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.expSizer = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer3.Add( self.expSizer, 0, wx.ALL, 5 )
		
		self.expView = wx.Button( self, wx.ID_ANY, u"View", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer3.Add( self.expView, 0, wx.ALL, 5 )
		
		self.expFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a .wig file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		bSizer3.Add( self.expFilePicker, 1, wx.ALL, 5 )
		
		
		expSizer.Add( bSizer3, 0, wx.EXPAND, 5 )
		
		self.list_exp = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		expSizer.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( expSizer, 1, wx.EXPAND, 5 )
		
		filesSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Files" ), wx.VERTICAL )
		
		bSizer141 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.displayButton = wx.Button( self, wx.ID_ANY, u"Display", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer141.Add( self.displayButton, 0, wx.ALL, 5 )
		
		self.tempButton = wx.Button( self, wx.ID_ANY, u"Add Test Files", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.tempButton.Hide()
		
		bSizer141.Add( self.tempButton, 0, wx.ALL, 5 )
		
		
		filesSizer.Add( bSizer141, 0, 0, 5 )
		
		self.list_files = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		filesSizer.Add( self.list_files, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( filesSizer, 1, wx.EXPAND, 5 )
		
		
		bSizer1.Add( bSizer4, 1, wx.EXPAND, 5 )
		
		methodSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Methods" ), wx.VERTICAL )
		
		methodSizer.SetMinSize( wx.Size( 200,-1 ) ) 
		self.m_staticText13 = wx.StaticText( self, wx.ID_ANY, u"Tn-Seq\nAnalysis", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText13.Wrap( -1 )
		self.m_staticText13.SetFont( wx.Font( 20, 74, 90, 92, False, "Sans" ) )
		
		methodSizer.Add( self.m_staticText13, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		methodChoiceChoices = [ u"[Choose Method]", u"Gumbel", u"HMM", u"Resampling" ]
		self.methodChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, methodChoiceChoices, 0 )
		self.methodChoice.SetSelection( 0 )
		methodSizer.Add( self.methodChoice, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.progressLabel = wx.StaticText( self, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.progressLabel.Wrap( -1 )
		methodSizer.Add( self.progressLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.progress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.progress.SetValue( 0 ) 
		methodSizer.Add( self.progress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.gumbelProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.gumbelProgress.SetValue( 0 ) 
		methodSizer.Add( self.gumbelProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.hmmProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.hmmProgress.SetValue( 0 ) 
		methodSizer.Add( self.hmmProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.resamplingProgress.SetValue( 0 ) 
		methodSizer.Add( self.resamplingProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.m_panel1.SetMinSize( wx.Size( 100,10 ) )
		
		methodSizer.Add( self.m_panel1, 0, wx.ALL, 5 )
		
		gumbelSection = wx.BoxSizer( wx.VERTICAL )
		
		self.gumbelInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have one control sample selected.\n2. Modify the options as desired.\n3. Click on the \"Run Gumbel\" button.\n4. Choose a name for the output file.\n5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.\n\n\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelInstructions.Wrap( 200 )
		gumbelSection.Add( self.gumbelInstructions, 0, wx.ALL, 5 )
		
		self.gumbelLabel = wx.StaticText( self, wx.ID_ANY, u"Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelLabel.Wrap( -1 )
		gumbelSection.Add( self.gumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer14 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer15 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer16 = wx.BoxSizer( wx.VERTICAL )
		
		self.gumbelSampleLabel = wx.StaticText( self, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelSampleLabel.Wrap( -1 )
		bSizer16.Add( self.gumbelSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.gumbelBurninLabel = wx.StaticText( self, wx.ID_ANY, u"Burn-In", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelBurninLabel.Wrap( -1 )
		bSizer16.Add( self.gumbelBurninLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.gumbelTrimLabel = wx.StaticText( self, wx.ID_ANY, u"Trim", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelTrimLabel.Wrap( -1 )
		bSizer16.Add( self.gumbelTrimLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.gumbelReadLabel = wx.StaticText( self, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelReadLabel.Wrap( -1 )
		bSizer16.Add( self.gumbelReadLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		
		bSizer15.Add( bSizer16, 1, wx.EXPAND, 5 )
		
		bSizer17 = wx.BoxSizer( wx.VERTICAL )
		
		self.gumbelSampleText = wx.TextCtrl( self, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer17.Add( self.gumbelSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
		
		self.gumbelBurninText = wx.TextCtrl( self, wx.ID_ANY, u"500", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer17.Add( self.gumbelBurninText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
		
		self.gumbelTrimText = wx.TextCtrl( self, wx.ID_ANY, u"1", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer17.Add( self.gumbelTrimText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
		
		gumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
		self.gumbelReadChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelReadChoiceChoices, 0 )
		self.gumbelReadChoice.SetSelection( 0 )
		bSizer17.Add( self.gumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
		
		
		bSizer15.Add( bSizer17, 1, wx.EXPAND, 5 )
		
		
		bSizer14.Add( bSizer15, 1, wx.EXPAND, 5 )
		
		
		gumbelSection.Add( bSizer14, 1, wx.EXPAND, 5 )
		
		self.gumbelButton = wx.Button( self, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
		gumbelSection.Add( self.gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( gumbelSection, 0, 0, 5 )
		
		hmmSection = wx.BoxSizer( wx.VERTICAL )
		
		self.hmmInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have one control sample selected.\n2. Click on the \"Run HMM\" button.\n3. Choose a name for the output file. Note: An additional file with the gene-level analysis will also be created.\n4. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.\n\n\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmInstructions.Wrap( 200 )
		hmmSection.Add( self.hmmInstructions, 0, wx.ALL, 5 )
		
		self.hmmLabel = wx.StaticText( self, wx.ID_ANY, u"HMM Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmLabel.Wrap( -1 )
		hmmSection.Add( self.hmmLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.hmmButton = wx.Button( self, wx.ID_ANY, u"Run HMM", wx.DefaultPosition, wx.DefaultSize, 0 )
		hmmSection.Add( self.hmmButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( hmmSection, 0, 0, 5 )
		
		resamplingSection = wx.BoxSizer( wx.VERTICAL )
		
		self.resamplingInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have selected at least one control sample and one experimental sample.\n2. Click on the \"Run Resampling\" button.\n3. Choose a name for the output file.\n4. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.\n\n\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingInstructions.Wrap( 200 )
		resamplingSection.Add( self.resamplingInstructions, 0, wx.ALL, 5 )
		
		self.resamplingLabel = wx.StaticText( self, wx.ID_ANY, u"Resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingLabel.Wrap( -1 )
		resamplingSection.Add( self.resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingButton = wx.Button( self, wx.ID_ANY, u"Run Resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( resamplingSection, 0, 0, 5 )
		
		
		bSizer1.Add( methodSizer, 0, wx.EXPAND, 5 )
		
		
		self.SetSizer( bSizer1 )
		self.Layout()
		self.statusBar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.annotationFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.annotationFileFunc )
		self.ctrlRemoveButton.Bind( wx.EVT_BUTTON, self.ctrlRemoveFunc )
		self.ctrlView.Bind( wx.EVT_BUTTON, self.ctrlViewFunc )
		self.ctrlFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadCtrlFileFunc )
		self.expSizer.Bind( wx.EVT_BUTTON, self.expRemoveFunc )
		self.expView.Bind( wx.EVT_BUTTON, self.expViewFunc )
		self.expFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadExpFileFunc )
		self.displayButton.Bind( wx.EVT_BUTTON, self.displayFileFunc )
		self.tempButton.Bind( wx.EVT_BUTTON, self.tempFileFunc )
		self.methodChoice.Bind( wx.EVT_CHOICE, self.MethodSelectFunc )
		self.gumbelButton.Bind( wx.EVT_BUTTON, self.RunGumbelFunc )
		self.hmmButton.Bind( wx.EVT_BUTTON, self.RunHMMFunc )
		self.resamplingButton.Bind( wx.EVT_BUTTON, self.RunResamplingFunc )
	
	def __del__( self ):
		pass
	
	
	# Virtual event handlers, overide them in your derived class
	def annotationFileFunc( self, event ):
		event.Skip()
	
	def ctrlRemoveFunc( self, event ):
		event.Skip()
	
	def ctrlViewFunc( self, event ):
		event.Skip()
	
	def loadCtrlFileFunc( self, event ):
		event.Skip()
	
	def expRemoveFunc( self, event ):
		event.Skip()
	
	def expViewFunc( self, event ):
		event.Skip()
	
	def loadExpFileFunc( self, event ):
		event.Skip()
	
	def displayFileFunc( self, event ):
		event.Skip()
	
	def tempFileFunc( self, event ):
		event.Skip()
	
	def MethodSelectFunc( self, event ):
		event.Skip()
	
	def RunGumbelFunc( self, event ):
		event.Skip()
	
	def RunHMMFunc( self, event ):
		event.Skip()
	
	def RunResamplingFunc( self, event ):
		event.Skip()
	

