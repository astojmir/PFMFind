#! /usr/bin/python
from ShortFrags.GUI.main_gui import MainGui
import Pmw
import sys

root = Pmw.initialise()
root.title("Short Fragment Toolbox")
mg = MainGui(root, dict=globals())
root.mainloop()
