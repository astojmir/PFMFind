#
# Copyright (C) 2005-2006 Victoria University of Wellington
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


import Tkinter, Pmw, string

_ffont = Pmw.logicalfont('Fixed')
BSIZE = 10
def col_f(x0):
    return x0 + x0 // BSIZE

def col_f_inv(x):
    z = x % (BSIZE+1)
    if z == BSIZE:
        z -= 1
    return BSIZE * (x // (BSIZE+1)) + z

def default_set_func(l, f, i):
    pass

class ScrolledSeq(Pmw.ScrolledText):
    def __init__(self, parent, query_seq,
                 set_func=default_set_func):
        Pmw.ScrolledText.__init__(self, parent,
                                  text_height = 3,
                                  text_wrap='none',
                                  text_font = _ffont,
                                  text_padx = 4,
                                  text_pady = 4,
                                  text_state = 'disabled',
                                  text_cursor = 'left_ptr',
                                  hscrollmode = 'static')
        self.state = (None, None, None)
        self.qseq = query_seq
        self.set_func = set_func
        self.tag_config('yellow', background = 'yellow')
        self.tag_bind('query_click', "<Button-1>", self._click_query)
        if self.qseq:
            self._set_query_sequence()

    def _click_query(self, event):
        pos = string.split(self.index( \
            "@%d,%d" % (event.x, event.y)),'.')
        x1 = col_f_inv(int(pos[1]))
        self.set_func(self.state[0], x1, self.state[2])

    def _set_query_sequence(self):
        self.configure(text_state = 'normal')
        self.delete(1.0, 'end') 
        line_len = col_f(len(self.qseq)-1) + 2
        self.insert('end', ' '*(line_len-1))
        self.insert('end', '\n')
        for i in xrange(0, len(self.qseq), BSIZE):
            self.insert('end', self.qseq[i:i+BSIZE])
            self.insert('end', ' ')
        self.insert('end', '\n')
        for i in xrange(0, len(self.qseq), BSIZE):
            self.insert('end', '%-*d' % (BSIZE+1, i))
        self.configure(text_state = 'disabled')

    def _highlight_current(self):
        l, x = self.state[0:2]
        pos0 = "3.%d" % col_f(x)
        pos1 = "3.%d" % col_f(x+l)
        self.tag_remove('yellow', '3.0', 'end')
        self.tag_add('yellow', pos0, pos1)

    def reset(self, state, qseq=None):

        if state == self.state: return
        self.state = state

        if qseq:
            self.qseq = qseq
            self._set_query_sequence()
        self._highlight_current()

        # Set clicking tags
        self.tag_remove('query_click', '2.0', 'end')
        for x in xrange(len(self.qseq) - self.state[0] + 1):
            pos = "2.%d" % col_f(x)
            self.tag_add('query_click', pos)

        self.see("2.%d" % col_f(self.state[1]))

