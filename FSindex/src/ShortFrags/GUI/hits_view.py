from ShortFrags.GUI.view import View
from ShortFrags.Expt.hit_list import *
import Pmw
import Tkinter
import string

def _row(ind):
    pos = string.split(ind, '.')
    return int(pos[0])
   

class HitsView(Tkinter.Frame, View):
    def __init__(self, parent=None, size=500):
        Tkinter.Frame.__init__(self, parent)
        View.__init__(self)
        self.pack_propagate(0)
        self.ffont = Pmw.logicalfont('Fixed')
        self.wMenu = Tkinter.Frame(self)
        self.wMenu.pack(anchor='w')


        self.HL = None
        self.sorts = {'Distance': HitList.sort_by_distance,
                      'Similarity': HitList.sort_by_similarity,
                      'Sequence': HitList.sort_by_seq,
                      'Seqid': HitList.sort_by_seqid,
                      }

        self.vSortVar = Tkinter.StringVar()
        self.wSortMenu = Pmw.OptionMenu(self.wMenu,
                                        labelpos = 'w',
                                        label_text = 'Sort by:',
                                        menubutton_textvariable = self.vSortVar,
                                        items = self.sorts.keys(),
                                        menubutton_width = 10,
                                        command = lambda a: self._show_hits())
        self.wSortMenu.pack(side='left')
        self.vSortVar.set('Distance')
        self.vSortIncr = Tkinter.IntVar()
        self.wSortIncr = Tkinter.Checkbutton(self.wMenu, text="Increasing",
                                             variable=self.vSortIncr,
                                             command = self._show_hits)
        self.wSortIncr.pack(side='left', padx=5)
        self.vSortIncr.set(1)
        self.wButtonBar = Pmw.ButtonBox(self.wMenu)
        self.wButtonBar.pack(side='left', expand=1, fill='x')


        self.wScText = Pmw.ScrolledText(self,
                                        vscrollmode='static',
                                        hscrollmode='static',
                                        usehullsize = 1,
                                        text_wrap='none',
                                        text_font = self.ffont,
                                        text_padx = 4,
                                        text_pady = 4,
                                        text_state = 'disabled',
                                        text_cursor = 'left_ptr',
                                        )
        self.wScText.pack()
        self.wScText.tag_bind('click', "<Button-1>", self._click_event)

        for text in ['Header', 'Summary', 'Details', 'Performance']:
            self.wButtonBar.add(text, command=lambda txt=text: self._scroll_to_details(txt))

        self.ind = {'Header':'1.0', 'Summary':'1.0', 'Details':'1.0', 'Performance':'1.0'}

    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)
        self.update_idletasks()
        hght = self.winfo_height() - self.wMenu.winfo_height()
        self.wScText.configure(hull_height=hght, hull_width=width)

    def reset(self, state={}):
        self.FE = state['FE']
        l, f = state['coords']
        i = state['iter']
        iters = self.FE.data.get_frag(l,f)

        if i < len(iters):
            self.HL = iters[i]
        else:
            self.HL = None
            self.wScText.tag_remove('click', '1.0', 'end')
            self.wScText.configure(text_state = 'normal')
            self.wScText.setvalue("No search hits.")
            self.wScText.configure(text_state = 'disabled')
            self.ind = {'Header':'1.0', 'Summary':'1.0', 'Details':'1.0', 'Performance':'1.0'}
            return

        self._show_hits()

    def _show_hits(self):
        if self.HL == None:
            return
        HL = self.HL
        self.wScText.tag_remove('click', '1.0', 'end')

        # Sort list
        self.sorts[self.vSortVar.get()](HL, self.vSortIncr.get())
        
        self.wScText.configure(text_state = 'normal')
        self.wScText.setvalue(HL.print_str())
        self.wScText.configure(text_state = 'disabled')

        self.ind['Header'] = self.wScText.search('***** Query Parameters *****', '1.0')
        self.ind['Summary'] = self.wScText.search('***** Summary *****', '1.0')
        self.ind['Details'] = self.wScText.search('***** Full Details *****', '1.0')
        self.ind['Performance'] = self.wScText.search('***** Index Performance *****', '1.0')

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
        y2 = d + 1 + 8*r
        self.wScText.yview('%d.0' % y2)

    def _scroll_to_details(self, text):
        y = _row(self.ind[text])
        self.wScText.yview('%d.0' % y)
