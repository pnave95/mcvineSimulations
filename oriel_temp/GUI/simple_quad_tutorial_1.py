'''
This program is a first experiment in using pyopengl in python 2
'''

#from OpenGL.GL import *
#from OpenGL.GLU import *
from OpenGLContext import testingcontext	# this depends on pyopengl, but does reportedly have bugs and is only usable from python2;  it is only designed for testing and learning


BaseContext = testingcontext.getInteractive()
from OpenGL.GL import *

class TestContext( BaseContext ):
	'''
	Rendering Context with custom viewpoint and Render 
	Note:  will have slightly different results as OpenGLContext automatically enables lighting
	'''
	initialPosition = (0,0,0)  # set initial camera position; tutorial does the re-positioning

	def Render( self, mode):
		'''Render the geometry for the scene.'''

		# prevent OpenGL from removing faces which face backward
		glDisable( GL_CULL_FACE )

		# move the drawing origin 6 units into the screen and 1.5 units to the left
		glTranslate(-1.5, 0.0, -6.0)

		# start the (legacy) geometry generation mode
		glBegin(GL_TRIANGLES)
		glVertex3f(0.0, 1.0, 0.0)
		glVertex3f(-1.0, -1.0, 0.0)
		glVertex3f(1.0, -1.0, 0.0)
		glEnd()

		# move the drawing origin again:  cumulative change is now (1.5, 0.0, 6.0)
		glTranslate(3.0, 0.0, 0.0)

		# start a different geometry generation mode
		glBegin(GL_QUADS)
		glVertex3f(-1.0,-1.0, 0.0)
		glVertex3f( 1.0,-1.0, 0.0)
		glVertex3f( 1.0, 1.0, 0.0)
		glVertex3f(-1.0, 1.0, 0.0)
		glEnd()


if __name__ == "__main__":
	TestContext.ContextMainLoop()



