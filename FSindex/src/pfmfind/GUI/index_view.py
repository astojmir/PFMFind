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

class IndexView(TextDisplay, View):
    def __init__(self, parent, PFMF_client):
        TextDisplay.__init__(self, parent, label='Index')
        View.__init__(self)
        self.PFMF_client = PFMF_client
        
    def reset(self, state):
        text = self.PFMF_client.index_data_str()
        self.set_text('Index', text)
