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


from pfmfind.GUI.view import View
from pfmfind.GUI.text_display import TextDisplay

class MatrixView(TextDisplay, View):
    def __init__(self, parent=None):
        TextDisplay.__init__(self, parent, label='Matrix')
        View.__init__(self)

    def reset(self, state):
        FE = state['FE']
        text = FE.print_matrix(state['matrix'],
                               state['conv_type'],
                               state['scale'],
                               state['weight'],
                               state['reg'],
                               list(state['coords']) + [state['iter']],
                               )
        self.set_text('Matrix', text)
        # if PSSM show distance as well
        if state['matrix'] == 'PSSM':
            text = FE.print_matrix(state['matrix'],
                                   "Quasi",
                                   state['scale'],
                                   state['weight'],
                                   state['reg'],
                                   list(state['coords']) + [state['iter']],
                                )
            self.insert_text(text)
            # Additional data
            text = FE.print_align(state['weight'],
                                  state['reg'],
                                  list(state['coords']) + [state['iter']],
                                  )
            self.insert_text(text)
        

