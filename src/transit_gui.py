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
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"TRANSIT", pos = wx.DefaultPosition, size = wx.Size( 1200,975 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
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
		
		self.ctrlView = wx.Button( self, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.ctrlView.Hide()
		
		bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )
		
		self.ctrlScatter = wx.Button( self, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.ctrlScatter.Hide()
		
		bSizer2.Add( self.ctrlScatter, 0, wx.ALL, 5 )
		
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
		
		self.expView = wx.Button( self, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.expView.Hide()
		
		bSizer3.Add( self.expView, 0, wx.ALL, 5 )
		
		self.expScatter = wx.Button( self, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.expScatter.Hide()
		
		bSizer3.Add( self.expScatter, 0, wx.ALL, 5 )
		
		self.expFilePicker = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a .wig file", u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		bSizer3.Add( self.expFilePicker, 1, wx.ALL, 5 )
		
		
		expSizer.Add( bSizer3, 0, wx.EXPAND, 5 )
		
		self.list_exp = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		expSizer.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( expSizer, 1, wx.EXPAND, 5 )
		
		filesSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Results Files" ), wx.VERTICAL )
		
		bSizer141 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.displayButton = wx.Button( self, wx.ID_ANY, u"Display Table", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer141.Add( self.displayButton, 0, wx.ALL, 5 )
		
		self.graphFileButton = wx.Button( self, wx.ID_ANY, u"Display Graph", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer141.Add( self.graphFileButton, 0, wx.ALL, 5 )
		
		self.addFileButton = wx.Button( self, wx.ID_ANY, u"Add Results File", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer141.Add( self.addFileButton, 0, wx.ALL, 5 )
		
		graphFileChoiceChoices = [ u"[Choose Figure to Plot]", u"Volcano Plot", u"Histogram of gene logFC" ]
		self.graphFileChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, graphFileChoiceChoices, 0 )
		self.graphFileChoice.SetSelection( 0 )
		bSizer141.Add( self.graphFileChoice, 0, wx.ALL, 5 )
		
		
		filesSizer.Add( bSizer141, 0, 0, 5 )
		
		self.list_files = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
		filesSizer.Add( self.list_files, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer4.Add( filesSizer, 1, wx.EXPAND, 5 )
		
		
		bSizer1.Add( bSizer4, 1, wx.EXPAND, 5 )
		
		methodSizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Methods" ), wx.VERTICAL )
		
		methodSizer.SetMinSize( wx.Size( 250,-1 ) ) 
		self.m_staticText13 = wx.StaticText( self, wx.ID_ANY, u"Tn-Seq\nAnalysis", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText13.Wrap( -1 )
		self.m_staticText13.SetFont( wx.Font( 20, 74, 90, 92, False, "Sans" ) )
		self.m_staticText13.Hide()
		
		methodSizer.Add( self.m_staticText13, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.logoImg = wx.StaticBitmap( self, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		methodSizer.Add( self.logoImg, 0, wx.ALL, 5 )
		
		self.versionLabel = wx.StaticText( self, wx.ID_ANY, u"v1.4.2", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.versionLabel.Wrap( -1 )
		self.versionLabel.SetFont( wx.Font( 10, 74, 90, 92, False, "Sans" ) )
		
		methodSizer.Add( self.versionLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		methodChoiceChoices = [ u"[Choose Method]", u"Gumbel", u"HMM", u"Resampling" ]
		self.methodChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, methodChoiceChoices, 0 )
		self.methodChoice.SetSelection( 0 )
		methodSizer.Add( self.methodChoice, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.progress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.progress.SetValue( 0 ) 
		methodSizer.Add( self.progress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.mainInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Choose the annotation file (\"prot table\") that corresponds to the datasets to be analyzed.\n2. Add the desired Control and Experimental datasets.\n3. (Optional) If you wish to visualize their read counts, select the desired datasets and click on the \"View\" button.\n4. Select the desired analysis method from the dropdown menu on the top-right of the window, and follow its instructions.\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.mainInstructions.Wrap( 250 )
		methodSizer.Add( self.mainInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.m_panel1.SetMinSize( wx.Size( 50,1 ) )
		
		methodSizer.Add( self.m_panel1, 0, wx.ALL, 5 )
		
		self.resamplingInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have selected at least one control sample and one experimental sample.\n2. Click on the \"Run Resampling\" button.\n3. Choose a name for the output file.\n4. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingInstructions.Wrap( 250 )
		methodSizer.Add( self.resamplingInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.hmmInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have one control sample selected.\n2. Click on the \"Run HMM\" button.\n3. Choose a name for the output file. Note: An additional file with the gene-level analysis will also be created.\n4. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.\n\n\n", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmInstructions.Wrap( 250 )
		methodSizer.Add( self.hmmInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.gumbelInstructions = wx.StaticText( self, wx.ID_ANY, u"Instructions:\n\n1. Make sure you have one control sample selected.\n2. Modify the options as desired.\n3. Click on the \"Run Gumbel\" button.\n4. Choose a name for the output file.\n5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelInstructions.Wrap( 250 )
		methodSizer.Add( self.gumbelInstructions, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.globalLabel = wx.StaticText( self, wx.ID_ANY, u"Global Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.globalLabel.Wrap( -1 )
		methodSizer.Add( self.globalLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer1431 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer1521 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer1621 = wx.BoxSizer( wx.VERTICAL )
		
		self.globalNTerminusLabel = wx.StaticText( self, wx.ID_ANY, u"Ignore N-Terminus %:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.globalNTerminusLabel.Wrap( -1 )
		bSizer1621.Add( self.globalNTerminusLabel, 0, wx.ALL, 5 )
		
		self.globalCTerminusLabel = wx.StaticText( self, wx.ID_ANY, u"Ignore C-Terminus %:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.globalCTerminusLabel.Wrap( -1 )
		bSizer1621.Add( self.globalCTerminusLabel, 0, wx.ALL, 5 )
		
		
		bSizer1521.Add( bSizer1621, 1, wx.EXPAND, 5 )
		
		bSizer1721 = wx.BoxSizer( wx.VERTICAL )
		
		self.globalNTerminusText = wx.TextCtrl( self, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1721.Add( self.globalNTerminusText, 0, wx.ALL, 5 )
		
		self.globalCTerminusText = wx.TextCtrl( self, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1721.Add( self.globalCTerminusText, 0, wx.ALL, 5 )
		
		
		bSizer1521.Add( bSizer1721, 1, wx.EXPAND, 5 )
		
		
		bSizer1431.Add( bSizer1521, 1, wx.EXPAND, 5 )
		
		
		methodSizer.Add( bSizer1431, 0, 0, 5 )
		
		gumbelSection = wx.BoxSizer( wx.VERTICAL )
		
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
		
		self.gumbelRepLabel = wx.StaticText( self, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.gumbelRepLabel.Wrap( -1 )
		bSizer16.Add( self.gumbelRepLabel, 1, wx.ALL, 5 )
		
		
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
		
		gumbelRepChoiceChoices = [ u"Sum", u"Mean" ]
		self.gumbelRepChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelRepChoiceChoices, 0 )
		self.gumbelRepChoice.SetSelection( 0 )
		bSizer17.Add( self.gumbelRepChoice, 0, wx.ALL|wx.EXPAND, 5 )
		
		
		bSizer15.Add( bSizer17, 1, wx.EXPAND, 5 )
		
		
		bSizer14.Add( bSizer15, 1, wx.EXPAND, 5 )
		
		
		gumbelSection.Add( bSizer14, 1, wx.EXPAND, 5 )
		
		self.gumbelButton = wx.Button( self, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
		gumbelSection.Add( self.gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( gumbelSection, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		hmmSection = wx.BoxSizer( wx.VERTICAL )
		
		self.hmmLabel = wx.StaticText( self, wx.ID_ANY, u"HMM Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmLabel.Wrap( -1 )
		hmmSection.Add( self.hmmLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer143 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer152 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer162 = wx.BoxSizer( wx.VERTICAL )
		
		self.hmmRepLabel = wx.StaticText( self, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.hmmRepLabel.Wrap( -1 )
		bSizer162.Add( self.hmmRepLabel, 1, wx.ALL, 5 )
		
		self.hmmLoessCheck = wx.CheckBox( self, wx.ID_ANY, u"LOESS Correction", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer162.Add( self.hmmLoessCheck, 0, wx.ALL, 5 )
		
		
		bSizer152.Add( bSizer162, 1, wx.EXPAND, 5 )
		
		bSizer172 = wx.BoxSizer( wx.VERTICAL )
		
		hmmRepChoiceChoices = [ u"Sum", u"Mean" ]
		self.hmmRepChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, hmmRepChoiceChoices, 0 )
		self.hmmRepChoice.SetSelection( 0 )
		bSizer172.Add( self.hmmRepChoice, 0, wx.ALL|wx.EXPAND, 5 )
		
		self.hmmLoessPrev = wx.Button( self, wx.ID_ANY, u"Plot LOESS", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer172.Add( self.hmmLoessPrev, 0, wx.ALL, 5 )
		
		
		bSizer152.Add( bSizer172, 1, wx.EXPAND, 5 )
		
		
		bSizer143.Add( bSizer152, 1, wx.EXPAND, 5 )
		
		
		hmmSection.Add( bSizer143, 1, wx.EXPAND, 5 )
		
		self.hmmButton = wx.Button( self, wx.ID_ANY, u"Run HMM", wx.DefaultPosition, wx.DefaultSize, 0 )
		hmmSection.Add( self.hmmButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( hmmSection, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		resamplingSection = wx.BoxSizer( wx.VERTICAL )
		
		self.resamplingLabel = wx.StaticText( self, wx.ID_ANY, u"Resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingLabel.Wrap( -1 )
		resamplingSection.Add( self.resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		bSizer142 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer151 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer161 = wx.BoxSizer( wx.VERTICAL )
		
		self.resamplingSampleLabel = wx.StaticText( self, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingSampleLabel.Wrap( -1 )
		bSizer161.Add( self.resamplingSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.resamplingNormLabel = wx.StaticText( self, wx.ID_ANY, u"Normalization", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingNormLabel.Wrap( -1 )
		bSizer161.Add( self.resamplingNormLabel, 0, wx.ALL, 5 )
		
		
		bSizer151.Add( bSizer161, 1, wx.EXPAND, 5 )
		
		bSizer171 = wx.BoxSizer( wx.VERTICAL )
		
		self.resamplingSampleText = wx.TextCtrl( self, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer171.Add( self.resamplingSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
		
		resamplingNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u"zinfnb", u"quantile", u"betageom", u"nonorm" ]
		self.resamplingNormChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, resamplingNormChoiceChoices, 0 )
		self.resamplingNormChoice.SetSelection( 0 )
		bSizer171.Add( self.resamplingNormChoice, 0, wx.ALL, 5 )
		
		
		bSizer151.Add( bSizer171, 1, wx.EXPAND, 5 )
		
		
		bSizer142.Add( bSizer151, 1, wx.EXPAND, 5 )
		
		
		resamplingSection.Add( bSizer142, 1, wx.EXPAND, 5 )
		
		self.resamplingLoessCheck = wx.CheckBox( self, wx.ID_ANY, u"Correction for Genome Positional Bias", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingLoessCheck, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingLoessPrev = wx.Button( self, wx.ID_ANY, u"Plot LOESS fit", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingLoessPrev, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingHistCheck = wx.CheckBox( self, wx.ID_ANY, u"Generate Resampling Histograms", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.resamplingHistCheck.SetValue(True) 
		resamplingSection.Add( self.resamplingHistCheck, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingAdaptiveCheck = wx.CheckBox( self, wx.ID_ANY, u"Adaptive Resampling (faster)", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingAdaptiveCheck, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingButton = wx.Button( self, wx.ID_ANY, u"Run Resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
		resamplingSection.Add( self.resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		methodSizer.Add( resamplingSection, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.progressLabel = wx.StaticText( self, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.progressLabel.Wrap( -1 )
		methodSizer.Add( self.progressLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.gumbelProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.gumbelProgress.SetValue( 0 ) 
		methodSizer.Add( self.gumbelProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.hmmProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.hmmProgress.SetValue( 0 ) 
		methodSizer.Add( self.hmmProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.resamplingProgress = wx.Gauge( self, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
		self.resamplingProgress.SetValue( 0 ) 
		methodSizer.Add( self.resamplingProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		
		bSizer1.Add( methodSizer, 0, wx.EXPAND, 5 )
		
		
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
		
		self.m_menubar1.Append( self.viewMenuItem, u"View" ) 
		
		self.SetMenuBar( self.m_menubar1 )
		
		self.statusBar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.annotationFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.annotationFileFunc )
		self.ctrlRemoveButton.Bind( wx.EVT_BUTTON, self.ctrlRemoveFunc )
		self.ctrlView.Bind( wx.EVT_BUTTON, self.allViewFunc )
		self.ctrlScatter.Bind( wx.EVT_BUTTON, self.scatterFunc )
		self.ctrlFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadCtrlFileFunc )
		self.expSizer.Bind( wx.EVT_BUTTON, self.expRemoveFunc )
		self.expView.Bind( wx.EVT_BUTTON, self.allViewFunc )
		self.expScatter.Bind( wx.EVT_BUTTON, self.scatterFunc )
		self.expFilePicker.Bind( wx.EVT_FILEPICKER_CHANGED, self.loadExpFileFunc )
		self.displayButton.Bind( wx.EVT_BUTTON, self.displayFileFunc )
		self.graphFileButton.Bind( wx.EVT_BUTTON, self.graphFileFunc )
		self.addFileButton.Bind( wx.EVT_BUTTON, self.addFileFunc )
		self.methodChoice.Bind( wx.EVT_CHOICE, self.MethodSelectFunc )
		self.gumbelButton.Bind( wx.EVT_BUTTON, self.RunGumbelFunc )
		self.hmmLoessPrev.Bind( wx.EVT_BUTTON, self.LoessPrevFunc )
		self.hmmButton.Bind( wx.EVT_BUTTON, self.RunHMMFunc )
		self.resamplingLoessPrev.Bind( wx.EVT_BUTTON, self.LoessPrevFunc )
		self.resamplingButton.Bind( wx.EVT_BUTTON, self.RunResamplingFunc )
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
	
	def __del__( self ):
		pass
	
	
	# Virtual event handlers, overide them in your derived class
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
	
	def RunGumbelFunc( self, event ):
		event.Skip()
	
	def LoessPrevFunc( self, event ):
		event.Skip()
	
	def RunHMMFunc( self, event ):
		event.Skip()
	
	
	def RunResamplingFunc( self, event ):
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
	
	
	

