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
        self.wScText.grid(row=1, sticky='ew')

        for text in ['Header', 'Summary', 'Details', 'Performance']:
            self.wButtonBar.add(text)

    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)
        self.update_idletasks()
        hght = self.winfo_height() - self.wMenu.winfo_height()
        self.wScText.configure(hull_height=hght, hull_width=width-40)

    def reset(self, state={}):
        self.FE = state['FE']
        l, f = state['coords']
        i = state['iter']
        iters = self.FE.get_iters(l,f)

        self.wMenu.grid_forget()
              
        if i < len(iters):
            self.HL = iters[i]
        else:
            self.HL = None
            self.wScText.configure(text_state = 'normal')
            self.wScText.setvalue("No search hits.")
            self.wScText.configure(text_state = 'disabled')
            return

        self._show_hits()

    def _show_hits(self):
         HL = self.HL
         self.wMenu.grid(row=0, sticky='ew')

         # Sort list
         self.sorts[self.vSortVar.get()](HL, self.vSortIncr.get())

         self.wScText.configure(text_state = 'normal')
         self.wScText.setvalue(HL.print_str(self.FE.Idata))
         self.wScText.configure(text_state = 'disabled')
         
        
    def _scroll_to(self, widget):
        pass

    def _scroll_to_details(self, i):
        pass
