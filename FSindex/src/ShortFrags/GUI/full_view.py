from ShortFrags.GUI.view import View
from ShortFrags.Expt.hit_list import HitList
from ShortFrags.Expt.hit_list import description
from array import array
import Pmw
import Tkinter
import string

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

class FullView(Tkinter.Frame, View):
    def __init__(self, parent=None, size=500, msg_func=None, set_func=lambda l,f,i:0):
        Tkinter.Frame.__init__(self, parent)
        View.__init__(self)
        self.pack_propagate(0)
        self.ffont = Pmw.logicalfont('Fixed')

        #Menu
        self.wMenu = Tkinter.Frame(self)
        self.vGroup = Tkinter.IntVar()
        self.wGroup = Tkinter.Checkbutton(self.wMenu, text="Group by cluster",
                                          variable=self.vGroup,
                                          command = self.reset)
        self.wGroup.pack(anchor='w',pady=5)
        self.vGroup.set(1)
        self.wMenu.pack(anchor='w')

        # Hit matrix
        self.wScText = Pmw.ScrolledText(self,
                                        vscrollmode='static',
                                        hscrollmode='static',
                                        columnheader = 1,
                                        usehullsize = 1,
                                        text_wrap='none',
                                        text_font = self.ffont,
                                        Header_font = self.ffont,
                                        text_padx = 4,
                                        text_pady = 4,
                                        text_state = 'disabled',
                                        Header_state = 'disabled',
                                        Header_padx = 4,
                                        text_cursor = 'left_ptr',
                                        Header_cursor = 'left_ptr',
                                        )
        #bdf = self.wScText.component('borderframe')
        #bdf.configure(relief='flat')
        self.wScText.pack()
        
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
        self.l = None
        self.iter = None
        self.qlen = None
        self.msg_func = msg_func
        self.set_func = set_func
        self.state = {}

        # Tags
        self.wScText.tag_bind('enter', "<Enter>", self._enter_event)
        self.wScText.tag_bind('enter', "<Motion>", self._enter_event,"+")
        self.wScText.tag_bind('leave', "<Leave>", self._leave_event)
        self.wScText.tag_bind('click', "<Button-1>", self._click_event)
        self.wScText.tag_config('blue', background = 'blue')
        self.wScText.tag_config('red', background = 'red')
        self.wScText.tag_config('green', background = 'green')
        self.wScText.tag_config('yellow', background = 'yellow')

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
                                 

    def reset(self, state={}):
        if state == {}:
            state = self.state
        else:
            self.state = state

        self.FE = state['FE']
        self.defline_func = self.FE.Idata.get_def
        self.rqseq = self.FE.rqseq
        l, f = state['coords']
        iter = state['iter']
        group = self.vGroup.get()

        if l != self.l or iter != self.iter or group != self.group:
            self.l = l
            self.iter = iter
            self.group = group
            # Get the dictionary
            if self.group:
                grp, dict = self.FE.full_view_by_cluster(l, iter)
            else:
                grp, dict = self.FE.full_view_by_seqid(l, iter)
            self.grp = grp
            self.dict = dict
            query_str, view_str = self._construct_strings(self.rqseq)

            # Put strings
            colhead = self.wScText.component('columnheader')
            colhead.configure(state = 'normal')
            colhead.delete(1.0, 'end') 
            colhead.insert('end', query_str)
            colhead.configure(state = 'disabled')

            self.wScText.configure(text_state = 'normal')
            self.wScText.setvalue(view_str)
            self.wScText.configure(text_state = 'disabled')

            self.wHitsText.configure(text_state = 'normal')
            self.wHitsText.setvalue("")
            self.wHitsText.configure(text_state = 'disabled')

            self._construct_tags()

        # Highlight query
        self._highlight_query(f)

    def clear(self):
        self.state = {}
        self.l = None
        self.iter = None
        self.vGroup.set(1)

        colhead = self.wScText.component('columnheader')
        colhead.configure(state = 'normal')
        colhead.delete(1.0, 'end') 
        colhead.configure(state = 'disabled')
        
        self.wScText.configure(text_state = 'normal')
        self.wScText.setvalue("")
        self.wScText.configure(text_state = 'disabled')
        
        self.wHitsText.configure(text_state = 'normal')
        self.wHitsText.setvalue("")
        self.wHitsText.configure(text_state = 'disabled')
        
    def _construct_strings(self, qseq):
        # Calculate total length of the main string
        line_len = col_f(len(qseq)-1) + 2
        tot_len = line_len * len(self.grp)

        # Query sequence - introduce gaps
        a = array('c', ' '*(line_len-1))
        for i in range(0, len(qseq), BSIZE):
            col_new = col_f(i)
            a[col_new:col_new+BSIZE] = array('c', qseq[i:i+BSIZE])
        query_str = a.tostring()

        # Main string
        if len(self.grp) == 0:
            return (query_str, " " * len(query_str))

        a = array('c', ' '*tot_len)
        for i, (cluster, seqid) in enumerate(self.grp):
            offset = i*line_len
            for j in self.dict[seqid].keys():
                col_new = col_f(j)
                if self.dict[seqid][j][1]:
                    a[offset + col_new] = '?'
                else:
                    a[offset + col_new] = '*'
            a[offset + line_len - 1] = '\n'

        # No '\n' on the last line
        a[tot_len-1] = ' '
        return (query_str, a.tostring())

    def _construct_tags(self):
        # Tags

        # Remove old tags
        for tag in self.wScText.tag_names():
            self.wScText.tag_remove(tag, '1.0', 'end')
            
        colhead = self.wScText.component('columnheader')
        colhead.tag_add('query_click', '1.0', 'end')
        if len(self.grp) == 0: return
        
        # Set new tags
        self.wScText.tag_add('enter', '1.0', 'end')
        self.wScText.tag_add('leave', '1.0', 'end')
        self.wScText.tag_add('click', '1.0', 'end')

        for i, (cluster, seqid) in enumerate(self.grp):
            for j in self.dict[seqid].keys():
                col_new = col_f(j)
                pos = "%d.%d" % (1+i, col_new)
                if self.dict[seqid][j][1]:
                    self.wScText.tag_add('blue', pos)
                else:
                    self.wScText.tag_add('green', pos)

    def _click_query(self, event):
        colhead = self.wScText.component('columnheader')
        pos = string.split(colhead.index("@%d,%d" % (event.x, event.y)),'.')
        x1 = col_f_inv(int(pos[1]))
        self.wHitsText.configure(text_state = 'normal')
        self.wHitsText.setvalue("")
        self.wHitsText.configure(text_state = 'disabled')
        self.set_func(self.l, x1, self.iter)

    def _get_pos(self, x, y):
        pos = string.split(self.wScText.index("@%d,%d" % (x,y)),'.')
        return (int(pos[1]), int(pos[0]) - 1)

    def _enter_event(self, event):
        (x1,y1) = self._get_pos(event.x, event.y)
        if y1 < len(self.grp):
            seqid = self.grp[y1][1]
            self.msg_func(self.defline_func(seqid))        

    def _leave_event(self, event):
        self.msg_func("")        
        
    def _click_event(self, event):
        x1,y1 = self._get_pos(event.x, event.y)
        x1 = col_f_inv(x1)
        if y1 >= len(self.grp): return
        seqid = self.grp[y1][1]
        if x1 in self.dict[seqid]:
            txt = description(self.dict[seqid][x1][0],
                              self.rqseq[x1:x1+self.l],
                              x1)
        else:
            txt = ""
        self.wHitsText.configure(text_state = 'normal')
        self.wHitsText.setvalue(txt)
        self.wHitsText.configure(text_state = 'disabled')

        # Change counter
        self.set_func(self.l, x1, self.iter)

    def _highlight_query(self, x):
        colhead = self.wScText.component('columnheader')
        pos0 = "1.%d" % col_f(x)
        pos1 = "1.%d" % col_f(x+self.l)
        colhead.tag_remove('yellow', '1.0', 'end')
        colhead.tag_add('yellow', pos0, pos1)
