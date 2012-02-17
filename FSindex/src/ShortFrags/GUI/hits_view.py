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


from ShortFrags.GUI.view import View
from ShortFrags.Expt.hit_list import HitList
import Pmw, Tkinter, string, threading
from ShortFrags.Expt.matrix import ScoreMatrix, ProfileMatrix
from ShortFrags.GUI.ScrolledSeq import ScrolledSeq, col_f
from ShortFrags.GUI.optionmenu import OptionMenu


def _row(ind):
    pos = string.split(ind, '.')
    return int(pos[0])


def _default_set_func():
    pass


class HitsView(Tkinter.Frame, View):

    _scroll_buttons = ['Header', 'Summary', 'Details', 'Matrix']
    _scroll_marks = ['***** Query Parameters *****',
                     '***** Summary *****',
                     '***** Full Details *****',
                     '***** Score Matrix *****',
                     ]

    _sorts = {'Distance': HitList.sort_by_distance,
             'Similarity': HitList.sort_by_similarity,
             'Sequence': HitList.sort_by_seq,
             'Seqid': HitList.sort_by_seqid,
             }

    def __init__(self, parent, PFMF_client,
                 set_func = _default_set_func):

        Tkinter.Frame.__init__(self, parent)
        View.__init__(self)
        self.pack_propagate(0)


        self.HL = None
        self.thread_lock = threading.Lock()
        self.db_access_lock = threading.Lock()
        self.current_thread = 0
        self.HL_done = False

        self.state=(None, None, None)
        self.PFMF_client = PFMF_client
        self.set_func = set_func

        # ************************************************************
        # ******* Options on the left ********************************
        # ************************************************************

        self.wMenu = Tkinter.Frame(self, width=130)
        self.wMenu.pack(side='left', anchor='nw', fill='y')
        self.wMenu.pack_propagate(0)


        # ***** Sorting Options **************************************

        self.wSortOpts = Pmw.Group(self.wMenu,
                                   tag_text='Sort by')
        self.wSortOpts.pack(anchor='nw', fill='x', padx=5, pady=7)
        self.vSortVar = Tkinter.StringVar()
        self.wSortMenu = OptionMenu(self.wSortOpts.interior(),
                                    menubutton_textvariable = self.vSortVar,
                                    items = self._sorts.keys(),
                                    menubutton_width = 10,
                                    command = lambda a: self._show_hits())
        self.wSortMenu.pack(anchor='w', padx=5, pady=3)
        self.vSortVar.set('Distance')
        self.vSortIncr = Tkinter.IntVar()
        self.wSortIncr = Tkinter.Checkbutton(self.wSortOpts.interior(),
                                             text="Increasing",
                                             variable=self.vSortIncr,
                                             command = self._show_hits)
        self.wSortIncr.pack(anchor='w', padx=5, pady=3)
        self.vSortIncr.set(1)


        # ***** Scroll Buttons ***************************************

        self.wScrollOpts = Pmw.Group(self.wMenu,
                                   tag_text='Scroll to')
        self.wScrollOpts.pack(anchor='nw', fill='x', padx=5, pady=7)
        self.wButtonBar = Pmw.ButtonBox(self.wScrollOpts.interior(),
                                        orient='vertical')
        self.wButtonBar.pack(anchor='nw', padx=5, pady=5, fill='x')


        # ************************************************************
        # ******* Sequence and Hits on the right *********************
        # ************************************************************

        self.wSeqHit =  Tkinter.Frame(self)
        self.wSeqHit.pack(side='left', anchor='nw',
                          fill='both', expand = 1)

        # ***** Full sequence viewer *********************************

        self.wSeqText = ScrolledSeq(self.wSeqHit,
                                    self.PFMF_client.query_sequence,
                                    self.set_func)
        self.wSeqText.pack(anchor='nw', fill='x')

        self.wSeqText.tag_config('orchid',
                                 background = 'orchid')
        self.wSeqText.tag_config('royal blue',
                                 background = 'royal blue')
        self.wSeqText.tag_config('medium turquoise',
                                 background = 'medium turquoise')
        self.wSeqText.tag_config('medium aquamarine',
                                 background = 'medium aquamarine')


        # ***** Hit list ***********************************************

        self.wScText = Pmw.ScrolledText(self.wSeqHit,
                                        vscrollmode='static',
                                        hscrollmode='static',
                                        text_wrap='none',
                                        text_font = self.ffont,
                                        text_padx = 4,
                                        text_pady = 4,
                                        text_state = 'disabled',
                                        text_cursor = 'left_ptr',
                                        )
        self.wScText.pack(anchor='nw', fill='both', expand = 1)
        self.wScText.tag_bind('click', "<Button-1>", self._click_event)

        for text in self._scroll_buttons:
            self.wButtonBar.add(text, command=lambda txt=text: self._scroll_to_details(txt))

        self.ind = dict.fromkeys(self._scroll_buttons, '1.0')

    def reset(self, state=(None, None, None)):

        if state == self.state: return

        self.state = length, fragment, iteration = state

        self.wSeqText.reset(state, self.PFMF_client.query_sequence)
        self._highlight_searches()

        self.wScText.tag_remove('click', '1.0', 'end')
        self.wScText.configure(text_state = 'normal')
        self.wScText.setvalue("Retrieving hits...")
        self.wScText.configure(text_state = 'disabled')
        self.ind = dict.fromkeys(self._scroll_buttons, '1.0')


        # Start a new thread
        self.thread_lock.acquire()
        self.current_thread += 1
        self.HL_done = False
        thr = threading.Thread(target=self._db_select_hits_thread,
                               args=(self.current_thread,
                                     length, fragment, iteration))
        thr.start()
        self.thread_lock.release()
        self._process_hits(self.current_thread)

    def _process_hits(self, thread_id):

        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()

        if not is_current_thread: return

        if self.HL_done:
            if self.HL == None:
                self.wScText.configure(text_state = 'normal')
                self.wScText.setvalue("No search hits.")
                self.wScText.configure(text_state = 'disabled')
            else:
                self._show_hits()
        else:
            self.after(100, self._process_hits, thread_id)

    def _db_select_hits_thread(self, thread_id, length,
                               fragment, iteration):

        self.db_access_lock.acquire()

        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()

        if is_current_thread:
            self.HL = self.PFMF_client.select_lif_search(length,
                                                         iteration,
                                                         fragment,
                                                         True)

        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()
        if is_current_thread:
            self.HL_done = True

        self.db_access_lock.release()


    def _show_hits(self):
        if self.HL == None: return

        HL = self.HL
        self.wScText.tag_remove('click', '1.0', 'end')

        # Sort list
        self._sorts[self.vSortVar.get()](HL, self.vSortIncr.get())

        self.wScText.configure(text_state = 'normal')
        text, self._offsets = HL.print_str(get_offsets=True, qs=self.state[1])
        self.wScText.setvalue(text)

        self.wScText.appendtext('***** Score Matrix *****\n')
        if hasattr(HL, 'matrix') and HL.matrix != None:
            M = ProfileMatrix(HL.matrix.pssm)
            self.wScText.appendtext(M.__str__())
        elif hasattr(HL, 'matrix_name') and isinstance(HL.matrix_name, str):
            self.wScText.appendtext(HL.matrix_name)
        self.wScText.configure(text_state = 'disabled')

        for sb, sm in zip(self._scroll_buttons, self._scroll_marks):
            self.ind[sb] = str(self.wScText.search(sm, '1.0'))

        s = _row(self.ind['Summary'])
        d = _row(self.ind['Details'])
        self.wScText.tag_add('click','%d.0' % (s+2), '%d.0' % (d-1))

    def _get_pos(self, x, y):
        pos = string.split(self.wScText.index("@%d,%d" % (x,y)),'.')
        return (int(pos[1]), int(pos[0]))

    def _click_event(self, event):
        if self.HL == None:
            return
        x1,y1 = self._get_pos(event.x, event.y)

        s = _row(self.ind['Summary'])
        d = _row(self.ind['Details'])

##         if y1 < s+2 or y1 > d-2:
##             return
        r = y1 - s - 2
        y2 = d + 1 + self._offsets[r]
        self.wScText.yview('%d.0' % y2)

    def _scroll_to_details(self, text):
        y = _row(self.ind[text])
        self.wScText.yview('%d.0' % y)

    def _highlight_searches(self):
        self.wSeqText.tag_remove('orchid', '1.0', '2.0')
        self.wSeqText.tag_remove('royal blue', '1.0', '2.0')
        self.wSeqText.tag_remove('medium turquoise', '1.0', '2.0')
        self.wSeqText.tag_remove('medium aquamarine', '1.0', '2.0')

        for f in xrange(len(self.PFMF_client.query_sequence)):
            if f in self.PFMF_client.max_iters:
                pos1 = "1.%d" % col_f(f)
                pos2 = "2.%d" % col_f(f)
                l = self.state[0]
                if l in self.PFMF_client.max_iters[f]:
                    i = self.PFMF_client.max_iters[f][l]
                    if i > self.state[2]:
                        self.wSeqText.tag_add('orchid', pos2)
                    elif i == self.state[2]:
                        self.wSeqText.tag_add('royal blue', pos2)
                    else:
                        self.wSeqText.tag_add('medium turquoise', pos2)
                else:
                    self.wSeqText.tag_add('medium aquamarine', pos2)
