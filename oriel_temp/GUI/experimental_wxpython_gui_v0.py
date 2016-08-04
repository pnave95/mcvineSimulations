# experiment in wxpython GUI building

import wx

try:
	from wx import glcanvas
except ImportError:
	raise ImportError, "Required dependency wx.glcanvas not present"


try:
	from OpenGL.GL import *
except ImportError:
	raise ImportError, "Required dependency OpenGL not present"



class GLFrame(wx.Frame):
	'''A simple class for using OpenGL with wxPython.'''

	def __init__(self, parent, id, title, pos=wx.DefaultPosition, size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE, name='frame'):
		#
		# Forcing a specific style on the window
		#  Should this include styles passed?  (what does this mean?)
		style = wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE

		super(GLFrame, self).__init__(parent, id, title, pos, size, style, name)

		self.GLinitialized = False # meaning?
		attribList = (glcanvas.WX_GL_RGBA,  # RGBA
					  glcanvas.WX_GL_DOUBLEBUFFER,	# Double Buffered
					  glcanvas.WX_GL_DEPTH_SIZE, 24)	# 24 bit  (but what is depth size?)

		# Create the canvas
		self.canvas = glcanvas.GLCanvas(self, attribList=attribList)

		# Set the event handlers
		self.canvas.Bind(wx.EVT_ERASE_BACKGROUND, self.processEraseBackgroundEvent)
		self.canvas.Bind(wx.EVT_SIZE, self.processSizeEvent)
		self.canvas.Bind(wx.EVT_PAINT, self.processPaintEvent)


		# Move the window
		self.Centre()

		# Show the window
		#self.Show()

	# 
	# Canvas Proxy Methods
	def GetGLExtents(self):
		'''Get the extents of the OpenGL canvas.'''
		return self.canvas.GetClientSize()

	def SwapBuffers(self):
		'''Swap the OpenGL buffers.'''
		self.canvas.SwapBuffers


	#
	# wxPython Window Handlers

	def processEraseBackgroundEvent(self, event):
		'''Process the erase background event'''
		pass # Do nothing, to avoid flashing on MSWin

	def processSizeEvent(self, event):
		'''Process the resize event.'''
		if self.canvas.GetContext():
			# Make sure the frame is shown before calling SetCurrent
			self.Show()
			self.canvas.SetCurrent()

			size = self.GetGLExtents()
			self.OnReshape(size.width, size.height)
			self.canvas.Refresh(False)
		event.Skip()

	def processPaintEvent(self, event):
		'''Process the drawing event.'''
		self.canvas.SetCurrent()

		# This is a 'perfect' time to initialize OpenGL ... only if we need to 
		if not self.GLinitialized:
			self.OnInitGL()
			self.GLinitialized = True

		self.OnDraw()
		event.Skip()

	#
	# GLFrame OpenGL Event Handlers

	def OnInitGL(self):
		'''Initialize OpenGL for use in the window.'''
		glClearColor(0.2, 0.2, 1, 0.5)

	def OnReshape(self, width, height):
		'''Reshape the OpenGL viewport based on the dimensions of the window.'''
		glViewport(0, 0, width, height)

		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(-0.5, 0.5, -0.5, 0.5, -1, 1)

		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()

	def OnDraw(self, *args, **kwargs):
		"""Draw the window"""
		#glClear(GL_COLOR_BUFFER_BIT)

		# Drawing an example triangle in the middle of the screen
		glBegin(GL_TRIANGLES)
		glColor(0, 0, 0)
		glVertex(-.25, -.25)
		glVertex(.25, -.25)
		glVertex(0, .25)
		glEnd()

		self.SwapBuffers()








# create a class (itself an object) "Primary", which inherits from the base class wx.Frame:
class Primary(wx.Frame):

	def __init__(self, *args, **kwargs):
		super(Primary, self).__init__(*args, **kwargs)

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

		self.SetSize((600,400))
		self.SetTitle('Data Analysis')
		self.Centre()
		self.Show(True)

	def OnQuit(self, e):
		self.Close()



def main():
	#app = wx.App()
	#Primary(None)
	#app.MainLoop()

	#app = wx.PySimpleApp()
	app = wx.App()
	frame = GLFrame(None, -1, 'GL Window')
	frame.Show()

	app.MainLoop()
	app.Destroy()


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


