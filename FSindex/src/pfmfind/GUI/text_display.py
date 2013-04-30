#
# Copyright (C) 2004-2006 Victoria University of Wellington
#
# This file is part of the PFMFind module.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#


import Tkinter
import Pmw
from pfmfind.GUI.view import View

class TextDisplay(Pmw.ScrolledText):
    def __init__(self, parent=None, label=None, text=None):
        Pmw.ScrolledText.__init__(self, parent,
                # borderframe = 1,
                labelpos = 'n',
                label_text=label,
                usehullsize = 1,
                text_wrap='none',
                text_font = View.ffont,
                text_padx = 4,
                text_pady = 4)

        if text:
            self.insert('end', text)
        self.configure(text_state = 'disabled')

    def set_size(self, height=400, width=300):
        self.configure(hull_height=height, hull_width=width)

    def set_text(self, label=None, text=None):
        if label:
            self.configure(label_text=label)
        if text:
            self.configure(text_state = 'normal')
            self.setvalue(text)
            self.configure(text_state = 'disabled')

    def insert_text(self, text):
        self.configure(text_state = 'normal')
        self.insert('end', text)
        self.configure(text_state = 'disabled')
