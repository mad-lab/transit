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
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.DefaultPosition, size = wx.Size( 500,300 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer1 = wx.BoxSizer( wx.VERTICAL )
		
		self.text = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1.Add( self.text, 1, wx.ALL|wx.EXPAND, 5 )
		
		self.solveButton = wx.Button( self, wx.ID_ANY, u"Solve", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1.Add( self.solveButton, 1, wx.ALL|wx.EXPAND, 5 )
		
		self.clearButton = wx.Button( self, wx.ID_ANY, u"Clear", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1.Add( self.clearButton, 1, wx.ALL|wx.EXPAND, 5 )
		
		
		self.SetSizer( bSizer1 )
		self.Layout()
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.solveButton.Bind( wx.EVT_BUTTON, self.solveFunc )
		self.clearButton.Bind( wx.EVT_BUTTON, self.clearFunc )
	
	def __del__( self ):
		pass
	
	
	# Virtual event handlers, overide them in your derived class
	def solveFunc( self, event ):
		event.Skip()
	
	def clearFunc( self, event ):
		event.Skip()
	

