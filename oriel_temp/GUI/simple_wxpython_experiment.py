''' 
This is an experiment using wxPython
'''

import wx


# create a class (itself an object) "Example", which inherits from the base class wx.Frame:
class Example(wx.Frame):

	#def __init__(self, parent, title):
	def __init__(self, *args, **kwargs):
		#super(Example, self).__init__(parent, title=title, size=(450,450))
		super(Example, self).__init__(*args, **kwargs)

		#self.Move((800,200))
		#self.Centre()
		#self.Show()
		self.InitUI()

	def InitUI(self):

		menubar = wx.MenuBar()
		fileMenu = wx.Menu()
		fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
		menubar.Append(fileMenu, '&File')	# the & means the F is underlined and is an accelerator key (i.e. so that Alt+F is a shortcut)
		self.SetMenuBar(menubar)

		self.Bind(wx.EVT_MENU, self.OnQuit, fitem)

		self.SetSize((400,400))
		self.SetTitle('Simple menu')
		self.Centre()
		self.Show(True)

	def OnQuit(self, e):
		self.Close()



def main():
	app = wx.App()
	Example(None)
	app.MainLoop()


if __name__ == '__main__':

	# create an application object
	#app = wx.App()

	# create wx.Frame object:  a widget (a special container widget)
	#frame = wx.Frame(None, -1, 'simple wxpython app')
	# call Show() method to actually display frame on the screen
	#frame.Show()

	# create an instance of the Example class
	#Example(None, title='Size')


	# enter the main GUI loop
	#app.MainLoop()

	main()


