# experiment with mayavi2 and traits ui for gui building and visualization

from traits.api import HasTraits, Instance
from traitsui.api import View, Item
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor

# option 2
from mayavi.core.ui.mayavi_scene import MayaviScene 

# option 3 -- doesn't seem to work correctly
#from tvtk.pyface.api import Scene

class MyModel(HasTraits):
	scene = Instance(SceneModel, ())

	#view = View(Item('scene', height=400, show_label=False, editor=SceneEditor() ) )
	view = View(Item('scene', height=400, show_label=False, editor=SceneEditor(scene_class=MayaviScene) ) )
	#view = View(Item('scene', height=400, show_label=False, editor=SceneEditor(scene_class=Scene) ) )


MyModel().configure_traits()