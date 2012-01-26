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


class View:
    """
    Meant to be used as mixin class with something
    derived from Tkinter.Frame (or set_size needs
    to be overwritten).
    """

    ffont = 'TkFixedFont'

    def __init__(self):
        pass

    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)

    def reset(self, state=(None, None, None)):
        pass
