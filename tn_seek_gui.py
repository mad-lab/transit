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
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Tn-Seek", pos = wx.DefaultPosition, size = wx.Size( 1092,594 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer1 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer4 = wx.BoxSizer( wx.VERTICAL )
		
		orgSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Organism" ), wx.HORIZONTAL )
		
		self.m_staticText5 = wx.StaticText( self, wx.ID_ANY, u"Annotation File:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText5.Wrap( -1 )
		orgSizer.Add( self.m_staticText5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.annotationFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		orgSizer.Add( self.annotationFilePicker, 1, wx.ALL, 5 )
		
		
		bSizer4.Add( orgSizer, 1, wx.EXPAND, 5 )
		
		ctrlSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Control Samples" ), wx.VERTICAL )
		
		bSizer2 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.ctrlRemoveButton = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.ctrlRemoveButton, 0, wx.ALL, 5 )
		
		self.ctrlView = wx.Button( self, wx.ID_ANY, u"View", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )
		
		self.ctrlFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"Read Files (*.wig; *.txt;*.dat)|*.wig;*.txt;*.dat|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
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
		
		self.expFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a .wig file", u"Read Files (*.wig; *.txt;*.dat)|*.wig;*.txt;*.dat|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		bSizer3.Add( self.expFilePicker, 1, wx.ALL, 5 )
		
		
		expSizer.Add( bSizer3, 0, wx.EXPAND, 5 )
		
		self.list_exp = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		expSizer.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( expSizer, 1, wx.EXPAND, 5 )
		
		
		bSizer1.Add( bSizer4, 1, wx.EXPAND, 5 )
		
		methodSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Methods" ), wx.VERTICAL )
		
		methodSizer.SetMinSize( wx.Size( 200,-1 ) ) 
		methodChoiceChoices = [ wx.EmptyString, u"Gumbel", u"HMM", u"Resampling" ]
		self.methodChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, methodChoiceChoices, 0 )
		self.methodChoice.SetSelection( 0 )
		methodSizer.Add( self.methodChoice, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.m_panel1.SetMinSize( wx.Size( 100,50 ) )
		
		methodSizer.Add( self.m_panel1, 0, wx.ALL, 5 )
		
		gumbelSection = wx.BoxSizer( wx.VERTICAL )
		
		self.gumbelLabel = wx.StaticText( self, wx.ID_ANY, u"Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelLabel.Wrap( -1 )
		gumbelSection.Add( self.gumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer5 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.gumbelSampleLabel = wx.StaticText( self, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelSampleLabel.Wrap( -1 )
		bSizer5.Add( self.gumbelSampleLabel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.gumbelSampleText = wx.TextCtrl( self, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer5.Add( self.gumbelSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		
		gumbelSection.Add( bSizer5, 0, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer6 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.gumbelReadLabel = wx.StaticText( self, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelReadLabel.Wrap( -1 )
		bSizer6.Add( self.gumbelReadLabel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		gumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
		self.gumbelReadChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelReadChoiceChoices, 0 )
		self.gumbelReadChoice.SetSelection( 0 )
		bSizer6.Add( self.gumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		
		gumbelSection.Add( bSizer6, 0, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.gumbelButton = wx.Button( self, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
		gumbelSection.Add( self.gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( gumbelSection, 0, wx.EXPAND, 5 )
		
		hmmSection = wx.BoxSizer( wx.VERTICAL )
		
		self.hmmLabel = wx.StaticText( self, wx.ID_ANY, u"HMM Section", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmLabel.Wrap( -1 )
		hmmSection.Add( self.hmmLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.hmmButton = wx.Button( self, wx.ID_ANY, u"Run HMM", wx.DefaultPosition, wx.DefaultSize, 0 )
		hmmSection.Add( self.hmmButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( hmmSection, 0, wx.EXPAND, 5 )
		
		resamplingSection = wx.BoxSizer( wx.VERTICAL )
		
		self.resamplingLabel = wx.StaticText( self, wx.ID_ANY, u"Resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingLabel.Wrap( -1 )
		resamplingSection.Add( self.resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingButton = wx.Button( self, wx.ID_ANY, u"Run Resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( resamplingSection, 0, wx.EXPAND, 5 )
		
		
		bSizer1.Add( methodSizer, 0, wx.EXPAND, 5 )
		
		
		self.SetSizer( bSizer1 )
		self.Layout()
		self.m_statusBar1 = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.annotationFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.annotationFileFunc )
		self.ctrlRemoveButton.Bind( wx.EVT_BUTTON, self.ctrlRemoveFunc )
		self.ctrlView.Bind( wx.EVT_BUTTON, self.ctrlViewFunc )
		self.ctrlFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadCtrlFileFunc )
		self.expSizer.Bind( wx.EVT_BUTTON, self.expRemoveFunc )
		self.expView.Bind( wx.EVT_BUTTON, self.expViewFunc )
		self.expFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadExpFileFunc )
		self.methodChoice.Bind( wx.EVT_CHOICE, self.MethodSelectFunc )
		self.gumbelButton.Bind( wx.EVT_BUTTON, self.RunGumbelFunc )
	
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
	
	def MethodSelectFunc( self, event ):
		event.Skip()
	
	def RunGumbelFunc( self, event ):
		event.Skip()
	

