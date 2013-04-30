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
from pfmfind.GUI.optionmenu import OptionMenu
from pfmfind.Expt.hit_list import HitList
from pfmfind.Expt.hit_list import description
from array import array
from cStringIO import StringIO
import Pmw, Tkinter, string, bisect, threading


FEATURE = 'Uniprot Feature Keys'
FEATURE = 'Uniprot Keywords'
# Column conversion function (maps actual
# column to the column that takes into account
# the gaps we introduce at every tenth residue
BSIZE = 10
def col_f(x0):
    return x0 + x0 // BSIZE

def col_f_inv(x):
    z = x % (BSIZE+1)
    if z == BSIZE:
        z -= 1
    return BSIZE * (x // (BSIZE+1)) + z

class KeywordView(Tkinter.Frame, View):

    _sorts = ['Significance', 'Alphabetical']

    def __init__(self, parent, PFMF_client,
                 msg_func=None, set_func=lambda l,f,i:0):
        Tkinter.Frame.__init__(self, parent)
        View.__init__(self)
        self.pack_propagate(0)

        #Menu
        self.wMenu = Tkinter.Frame(self)

        self.wRng = Tkinter.Frame(self.wMenu, bd=2, relief='ridge')
        self.wRng.pack(side='left', anchor='w', padx=10, pady=5)
        Tkinter.Label(self.wRng, text='Fragment Range:').pack(side='left', padx=3)
        self.wRng1 = Pmw.EntryField(self.wRng, labelpos = 'w', label_text = 'From',
                                    entry_width = 4)
        self.wRng1.pack(side='left', pady=5)
        self.wRng2 = Pmw.EntryField(self.wRng, labelpos = 'w', label_text = 'To',
                                    entry_width = 4, hull_padx=2)
        self.wRng2.pack(side='left', pady=5)


        self.ontologies = []
        self.vOnt = Tkinter.StringVar()
        self.wOnt = OptionMenu(self.wMenu,
                               labelpos = 'w',
                               label_text = 'Ontology:',
                               menubutton_textvariable = self.vOnt,
                               items = self.ontologies,
                               menubutton_width = 20,
                               command = None,
                               hull_bd=2, hull_relief='ridge',
                               )
        self.wOnt.pack(side='left', anchor='w', padx=5, pady=5)


        self.vSort = Tkinter.StringVar()
        self.wSort = OptionMenu(self.wMenu,
                                labelpos = 'w',
                                label_text = 'Sort:',
                                menubutton_textvariable = self.vSort,
                                items = self._sorts,
                                menubutton_width = 15,
                                hull_bd=2, hull_relief='ridge',
                                )
        self.wSort.pack(side='left',  padx=10, pady=5)


        self.wBut = Tkinter.Button(self.wMenu, text='Compute',
                                   command=self.compute_view)
        self.wBut.pack(side='left', padx=15)

        self.wMenu.pack(anchor='w')

        # Hit matrix
        self.wScText = Pmw.ScrolledText(self,
                                        vscrollmode='static',
                                        hscrollmode='static',
                                        columnheader = 1,
                                        rowheader = 1,
                                        usehullsize = 1,
                                        text_wrap='none',
                                        text_font = self.ffont,
                                        Header_font = self.ffont,
                                        text_padx = 4,
                                        text_pady = 4,
                                        text_state='disabled',
                                        Header_state='disabled',
                                        Header_padx = 4,
                                        rowheader_pady = 4,
                                        text_cursor = 'left_ptr',
                                        Header_cursor = 'left_ptr',
                                        rowheader_width = 35,
                                        )
        self.wScText.pack(fill = 'both', expand = 1)

        # Details of hits - scrolled text
        self.wHitsText = Pmw.ScrolledText(self,
                                          usehullsize = 1,
                                          text_wrap='none',
                                          text_font = self.ffont,
                                          vscrollmode='static',
                                          text_state = 'disabled',
                                          text_cursor = 'left_ptr',
                                          )
        self.wHitsText.pack()

        # State
        self.PFMF_client = PFMF_client
        self.computed_len_iter = (None, None)
        self.VL = None
        self.hits_text = ""
        self.thread_lock = threading.Lock()
        self.db_access_lock = threading.Lock()
        self.current_thread = 0

        self.msg_func = msg_func
        self.set_func = set_func
        self.state = (None, None, None)

        # Tags
        self.wScText.tag_bind('enter', "<Enter>", self._enter_event)
        self.wScText.tag_bind('enter', "<Motion>", self._enter_event,"+")
        self.wScText.tag_bind('leave', "<Leave>", self._leave_event)
        self.wScText.tag_bind('click', "<Button-1>", self._click_event)
        self.wScText.tag_config('red', background = 'red')
        self.wScText.tag_config('orange', background = 'orange')
        self.wScText.tag_config('yellow', background = 'yellow')
        self.wScText.tag_config('green yellow', background = 'green yellow')
        self.wScText.tag_config('green', background = 'green')
        self.wScText.tag_config('cyan', background = 'cyan')
        self.wScText.tag_config('blue', background = 'blue')
        self.wScText.tag_config('purple', background = 'purple')

        colhead = self.wScText.component('columnheader')
        colhead.tag_config('yellow', background = 'yellow')
        colhead.tag_bind('query_click', "<Button-1>", self._click_query)

    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)
        self.update_idletasks()
        self.wHitsText.configure(hull_height=120, hull_width = width)
        self.wHitsText.update_idletasks()
        hght = self.winfo_height() - self.wMenu.winfo_height() -\
               self.wHitsText.component('hull').winfo_height()
        self.wScText.configure(hull_height=hght, hull_width=width)

    def compute_view(self):

        r1 = self.wRng1.getvalue()
        r2 = self.wRng2.getvalue()

        try:
            r1 = int(r1)
            r2 = int(r2)
        except ValueError:
            return

        length, fragment, iteration = self.state

        if r1 < 0 or r2 > len(self.rqseq)-length+1:
            return

        if (length, iteration) != self.computed_len_iter or \
           (r1,r2) != self.frag_rng or self.ont != self.vOnt.get():
            self.VL = None
            self.msg_func("Computing view.")
            self.wBut.configure(state='disabled')

            thr = threading.Thread(target=self._db_select_view_thread,
                                   args=(length, iteration,
                                         self.vOnt.get(), (r1, r2)))
            thr.start()
        else:
            self.VL = self.view_list

        self._process_view(length, iteration, (r1,r2))


    def _db_select_view_thread(self, *args, **kwargs):
        self.VL = self.PFMF_client.select_li_features(*args, **kwargs)

    def _process_view(self, length, iteration, frag_rng):
        # Check if ready
        if self.VL == None:
            self.after(100, self._process_view,
                       length, iteration, frag_rng)
            return


        self.wBut.configure(state='normal')
        self.msg_func("")
        if self.VL == []: return


        self.view_list = self.VL
        self.computed_len_iter = (length, iteration)
        self.frag_rng = frag_rng
        self.ont = self.vOnt.get()


        self.sort_type = self.vSort.get()

        HS = {}
        kw_accessions = []
        total_accessions = 0
        for i, (keyword, items) in enumerate(self.view_list):
            kw_accessions.append(0)
            for j, accessions in items:
                k = len(accessions)
                if j in HS:
                    HS[j] += k
                else:
                    HS[j] = k
                kw_accessions[-1] += k

        self.HL_sizes = HS

        # Sort view_list
        if self.sort_type == 'Alphabetical':
            self.view_list.sort()
        elif self.sort_type == 'Significance':
            tmp = zip(kw_accessions, self.view_list)
            tmp.sort()
            tmp.reverse()
            self.view_list = [t[1] for t in tmp]

        line_len = col_f(len(self.rqseq)-1) + 2
        ts = StringIO()
        ks = StringIO()

        for i, (keyword, items) in enumerate(self.view_list):
            a = array('c', ' '*line_len)

            for j, accessions in items:
                col_new = col_f(j)
                a[col_new] = '*'

            if i < len(self.view_list) - 1:
                a[line_len - 1] = '\n'
                ts.write(a.tostring())
                ks.write('%s\n' % keyword)
            else:
                ts.write(a.tostring())
                ks.write('%s' % keyword)

        self.kw_text = ks.getvalue()
        self.view_text = ts.getvalue()


        self._show_view(self.kw_text, self.view_text)



    def _reset_params(self, default=True):
        fragment = self.state[1]
        self.wRng1.setentry(str(fragment))
        self.wRng2.setentry(str(fragment+1))
        if default:
            self.vOnt.set(self.ontologies[0])
            self.vSort.set(self._sorts[0])
        else:
            self.vOnt.set(self.ont)
            self.vSort.set(self.sort_type)

    def _show_view(self, kw_text, view_text):

        colhead = self.wScText.component('columnheader')
        colhead.configure(state = 'normal')
        colhead.delete(1.0, 'end')

        rowhead = self.wScText.component('rowheader')
        rowhead.configure(state = 'normal')
        rowhead.delete(1.0, 'end')

        self.wScText.configure(text_state = 'normal')
        self.wScText.delete(1.0, 'end')

        line_len = col_f(len(self.rqseq)-1) + 2

        # Query sequence - introduce gaps
        # Note: I wasn't able to add another line and
        # have a good alignment with the main scrolled text
        a = array('c', ' '*(line_len-1))
        for i in range(0, len(self.rqseq), BSIZE):
            col_new = col_f(i)
            a[col_new:col_new+BSIZE] = array('c', self.rqseq[i:i+BSIZE])
        colhead.insert('end', a.tostring())
        colhead.configure(state = 'disabled')

        # Main string

        if view_text == None:
            empty = True
            kw_text = ""
            view_text = " " * line_len
        else:
            empty = False

        self.wScText.setvalue(view_text)
        rowhead.insert('end', kw_text)
        self.wScText.configure(text_state = 'disabled')
        rowhead.configure(state = 'disabled')

        self.wHitsText.configure(text_state = 'normal')
        self.wHitsText.setvalue("")
        self.wHitsText.configure(text_state = 'disabled')

        self._construct_tags(empty)
        self._highlight_query(self.state[1])


    def reset(self, state=(None, None, None)):

        if state == self.state: return
        self.state = state
        (length, fragment, iteration) = self.state

        self.rqseq = self.PFMF_client.query_sequence
        self.ontologies = self.PFMF_client.feature_types
        self.wOnt.setitems(self.ontologies)

        if (length, iteration) == self.computed_len_iter:
            self._show_view(self.kw_text, self.view_text)
            self._reset_params(default=False)
        else:
            self._show_view(None, None)
            self._reset_params(default=True)

        # Highlight query
        self._highlight_query(fragment)

        # Scroll to the selected fragment
        self.wScText.see("1.%d" % col_f(self.state[1]))

    def clear(self):
        self.state = (None, None, None)

        colhead = self.wScText.component('columnheader')
        colhead.configure(state='normal')
        colhead.delete(1.0, 'end')
        colhead.configure(state='disabled')

        rowhead = self.wScText.component('rowheader')
        rowhead.configure(state='normal')
        rowhead.delete(1.0, 'end')
        rowhead.configure(state='disabled')

        self.wScText.configure(text_state='normal')
        self.wScText.setvalue("")
        self.wScText.configure(text_state='disabled')

        self.wHitsText.configure(text_state='normal')
        self.wHitsText.setvalue("")
        self.wHitsText.configure(text_state='disabled')

    def _construct_tags(self, empty=False):
        # Tags

        # Remove old tags
        for tag in self.wScText.tag_names():
            self.wScText.tag_remove(tag, '1.0', 'end')

        colhead = self.wScText.component('columnheader')
        colhead.tag_add('query_click', '1.0', '2.0')
        if empty: return

        # Set new tags
        self.wScText.tag_add('enter', '1.0', 'end')
        self.wScText.tag_add('leave', '1.0', 'end')
        self.wScText.tag_add('click', '1.0', 'end')

        for i, (keyword, items) in enumerate(self.view_list):
            for j, accessions in items:
                col_new = col_f(j)
                pos = "%d.%d" % (1+i, col_new)

                # This is 'significance' of the dot for now
                kw_p = float(len(accessions)) / self.HL_sizes[j]

                if kw_p >= 0.5:
                    self.wScText.tag_add('purple', pos)
                elif kw_p >= 0.4:
                    self.wScText.tag_add('blue', pos)
                elif kw_p >= 0.3:
                    self.wScText.tag_add('cyan', pos)
                elif kw_p >= 0.2:
                    self.wScText.tag_add('green', pos)
                elif kw_p >= 0.1:
                    self.wScText.tag_add('green yellow', pos)
                elif kw_p >= 0.05:
                    self.wScText.tag_add('yellow', pos)
                elif kw_p >= 0.01:
                    self.wScText.tag_add('orange', pos)
                else:
                    self.wScText.tag_add('red', pos)


    def _click_query(self, event):
        colhead = self.wScText.component('columnheader')
        pos = string.split(colhead.index("@%d,%d" % (event.x, event.y)),'.')
        x1 = col_f_inv(int(pos[1]))
        self.wHitsText.configure(text_state='normal')
        self.wHitsText.setvalue("")
        self.wHitsText.configure(text_state='disabled')
        self.set_func(self.state[0], x1, self.state[2])

    def _get_pos(self, x, y):
        pos = string.split(self.wScText.index("@%d,%d" % (x,y)),'.')
        return (int(pos[1]), int(pos[0]) - 1)

    def _enter_event(self, event):
        (x1,y1) = self._get_pos(event.x, event.y)
        if y1 < len(self.view_list):
            kw = self.view_list[y1][0]
            self.msg_func(kw)

    def _leave_event(self, event):
        self.msg_func("")

    def _click_event(self, event):
        if self.hits_text == None: return

        x1,y1 = self._get_pos(event.x, event.y)
        x1 = col_f_inv(x1)
        if y1 >= len(self.view_list): return

        keyword, items = self.view_list[y1]
        length, fragment, iteration = self.state

        i = bisect.bisect(items, (x1, []))
        if i < len(items) and x1 == items[i][0]:
            self.wHitsText.configure(text_state='normal')
            self.wHitsText.setvalue("Retrieving hits...")
            self.wHitsText.configure(text_state='disabled')

            self.thread_lock.acquire()
            self.current_thread += 1
            self.hits_text = None
            thr = threading.Thread(target=self._db_select_hits_thread,
                                   args=(self.current_thread,
                                         items[i][1], length,
                                         x1, iteration),
                                   kwargs={'features': [self.vOnt.get()]})
            thr.start()
            self.thread_lock.release()
        else:
            self.hits_text = ""

        self._process_hits(self.current_thread)

    def _process_hits(self, thread_id):

        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()

        if not is_current_thread: return

        if self.hits_text == None:
            self.after(100, self._process_hits, thread_id)
        else:
            self.wHitsText.configure(text_state='normal')
            self.wHitsText.setvalue(self.hits_text)
            self.wHitsText.configure(text_state='disabled')

    def _db_select_hits_thread(self, thread_id, accessions, length,
                               fragment, iteration, features):

        self.db_access_lock.acquire()

        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()

        if is_current_thread:
            HL = self.PFMF_client.select_lif_search(length,
                                                    iteration,
                                                    fragment,
                                                    features=features)
        self.thread_lock.acquire()
        is_current_thread = (thread_id == self.current_thread)
        self.thread_lock.release()
        if is_current_thread:
            HL1 = [h for h in HL if h._bioentry_id in accessions]
            self.hits_text = description(HL1,
                                         self.rqseq[fragment:fragment+length],
                                         fragment)
        self.db_access_lock.release()


    def _highlight_query(self, x):
        l = self.state[0]
        colhead = self.wScText.component('columnheader')
        pos0 = "1.%d" % col_f(x)
        pos1 = "1.%d" % col_f(x+l)
        colhead.tag_remove('yellow', '1.0', 'end')
        colhead.tag_add('yellow', pos0, pos1)
