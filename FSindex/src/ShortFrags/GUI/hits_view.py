from ShortFrags.GUI.view import View
from ShortFrags.Expt.hit_list import *
import Pmw
import Tkinter


class HitsView(Tkinter.Frame, View):
    def __init__(self, parent=None, size=500):
        Tkinter.Frame.__init__(self, parent)
        View.__init__(self)
        self.grid_propagate(0)
        self.ffont = Pmw.logicalfont('Fixed')
        self.wMenu = Tkinter.Frame(self)
        self.wScFrame = Pmw.ScrolledFrame(self)
        bdf = self.wScFrame.component('borderframe')
        bdf.configure(relief='flat')

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

        #self.wMenu.grid(row=0, sticky='ew')
        
        self.wScFrame.grid(row=1, sticky='ew')
        self.wFrame = self.wScFrame.interior()
        self.wHeader = Tkinter.Label(self.wFrame,
                                     font=self.ffont,
                                     anchor='w',
                                     justify='left')
        stxt = "***** Summary *****\n" + \
               '  %4.4s  %-*.*s %4.4s %4.4s' % ('Rank', 57, 57,
                                               'Description', 'Dist', 'Sim')
        self.wSummaryTitle = Tkinter.Label(self.wFrame,
                                           font=self.ffont,
                                           anchor='w',
                                           justify='left',
                                           text=stxt)
        self.wSummaryButtons = []
        self.wFullTitle = Tkinter.Label(self.wFrame,
                                        font=self.ffont,
                                        anchor='w',
                                        justify='left',
                                        text='\n\n***** Full Details *****')
        self.wFullLines = []
        self.wPerf = Tkinter.Label(self.wFrame,
                                   font=self.ffont,
                                   anchor='w',
                                   justify='left')
        
        for i in range(size):
            self.wSummaryButtons.append(Tkinter.Button(
                self.wFrame,
                command = lambda i=i:self._scroll_to_details(i),
                font=self.ffont,
                anchor='w',
                justify='left',
                relief='groove'))
            self.wFullLines.append(Tkinter.Label(self.wFrame,
                                                 font=self.ffont,
                                                 anchor='w',
                                                 justify='left',
                                                 relief='ridge'))
            self.glist = []

        self.sections = {'Header': self.wHeader,
                         'Summary': self.wSummaryTitle,
                         'Details': self.wFullTitle,
                         'Performance': self.wPerf,
                         }
                         
        for text in ['Header', 'Summary', 'Details', 'Performance']:
            self.wButtonBar.add(text,
                                command = lambda w=self.sections[text]:
                                self._scroll_to(w))

    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)
        clp = self.wScFrame.component('clipper')
        clp.configure(height=height-40, width=width-40)
        
    def reset(self, state={}):
        self.FE = state['FE']
        l, f = state['coords']
        i = state['iter']
        iters = self.FE.get_iters(l,f)

        for w in self.glist:
            w.grid_forget()
              
        if i < len(iters):
            self.HL = iters[i]
        else:
            self.HL = None
            self.wHeader.configure(text='No search hits.')
            self.glist = [self.wHeader]
            self.wHeader.grid(row=0, column=0)
            return

        self._show_hits()

    def _show_hits(self):
        HL = self.HL

        # Sort list
        self.sorts[self.vSortVar.get()](HL, self.vSortIncr.get())

        if len(HL.hits) > self.size:
            old_size = size
            self.size = len(HL.hits)
            for j in range(old_size, self.size):
                self.wSummaryButtons.append(Tkinter.Button(self.wFrame,
                                                           font=self.ffont,
                                                           anchor='w',
                                                           justify='left',
                                                           relief='groove'))
                self.wFullLines.append(Tkinter.Label(self.wFrame,
                                                     font=self.ffont,
                                                     anchor='w',
                                                     justify='left',
                                                     relief='ridge'))

        
        self.wHeader.configure(text=HL.header_str())
        self.glist = [self.wMenu, self.wHeader, self.wSummaryTitle]

        for j, ht in enumerate(HL.hits):
            self.wSummaryButtons[j].configure(text=summary([ht], 1,
                                                           rank_offset=j)) 
            self.glist.append(self.wSummaryButtons[j])

        self.glist.append(self.wFullTitle)
        for j, ht in enumerate(HL.hits):
            self.wFullLines[j].configure(text=description([ht], HL.query_seq)) 
            self.glist.append(self.wFullLines[j])

        self.wPerf.configure(text='\n\n'+HL.perf_str(self.FE.Idata))
        self.glist.append(self.wPerf)

        for j, w in enumerate(self.glist):
            w.grid(row=j, sticky='ew')
        
    def _scroll_to(self, widget):
        y = float(widget.winfo_y()) / self.wFrame.winfo_height()
        self.wScFrame.yview('moveto', y)

    def _scroll_to_details(self, i):
        widget = self.wFullLines[i]
        self._scroll_to(widget)
