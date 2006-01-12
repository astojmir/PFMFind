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


import Tkinter, Pmw

def _null_func(*args, **kwargs):
    pass
                    
class StatusShow(Tkinter.Frame):
    def __init__(self, parent, PFMF_client, ufunc=_null_func,
                 set_func=_null_func, unset_func=_null_func):

        Tkinter.Frame.__init__(self, parent, bd=2, relief='groove')

        self.PC = PFMF_client
        self.update_func = ufunc
        self.set_func = set_func
        self.unset_func = unset_func

        self.wFrame = Tkinter.Frame(self)
        self.wFrame.pack(side='left', anchor='w', pady=2)
        self.wLabel = Tkinter.Label(self.wFrame, text="Query Fragment")
        self.wLabel.grid(row=0, column=0, sticky='w')
        self.wSeqLabel = Tkinter.Label(self.wFrame, text="",
                                       font=Pmw.logicalfont('Fixed'))
        self.wSeqLabel.grid(row=0, column=1, sticky='w')
                                            
        self.wFragLen = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Fragment\nLength',
                                    entry_width = 3,
                                    entryfield_value = ' ',
                                    )
        self.wCurFrag = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Current\nFragment',
                                    entry_width = 5,
                                    entryfield_value = ' ',
                                    )
        self.wCurIter = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Current\nIteration',
                                    entry_width = 3,
                                    entryfield_value = ' ',
                                    )
        for w in [self.wFragLen, self.wCurIter, self.wCurFrag]: 
            w.pack(side='right', anchor='w', padx=5)
        self.reset(init=True)

    def set_values(self, length, fragment, iteration=None):
        self.wFragLen.setvalue(str(length))
        self.wCurFrag.setvalue(str(fragment))
        if iteration == None:
            iteration = self._max_iters(length, fragment)
        self.wCurIter.setvalue(str(iteration))
        
    def get_values(self):
        l = int(self.wFragLen.getvalue())
        f = int(self.wCurFrag.getvalue())
        i = int(self.wCurIter.getvalue())
        return (l,f,i)


    def _max_iters(self, l, f):
        if f in self.PC.max_iters:
            return self.PC.max_iters[f].get(l, -1) + 1
        else:
            return 0


    def _set_qseq(self):
        l = int(self.wFragLen.getvalue())
        f = int(self.wCurFrag.getvalue())
        s = self.PC.query_sequence[f:f+l]
        self.wSeqLabel.configure(text=s)  

    def _validate_none(self, text):
        print "validate_none", text
        if text == '':
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _validate_len(self, text):
        try:
            v = int(text)
        except ValueError:
            return Pmw.ERROR
        
        if v >= self.PC.min_len and v < self.PC.max_len:
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _validate_frag(self, text):
        try:
            v = int(text)
        except ValueError:
            return Pmw.ERROR

        l = int(self.wFragLen.getvalue())
        if v >= 0 and v <= len(self.PC.query_sequence) - l:
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _validate_iter(self, text):
        try:
            v = int(text)
        except ValueError:
            return Pmw.ERROR
        
        l = int(self.wFragLen.getvalue())
        f = int(self.wCurFrag.getvalue())
        if v >= 0 and v <= self._max_iters(l, f):
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _len_changed(self):
        l = int(self.wFragLen.getvalue())
        f = int(self.wCurFrag.getvalue())
        if f > len(self.PC.query_sequence) - l:
            f = len(self.PC.query_sequence) - l
            self.wCurFrag.setvalue(str(f))
        else:
            self._frag_changed()

    def _frag_changed(self):
        l = int(self.wFragLen.getvalue())
        f = int(self.wCurFrag.getvalue())
        i = int(self.wCurIter.getvalue())
        if i > self._max_iters(l, f):
            i = self._max_iters(l, f)
            self.wCurIter.setvalue(str(i))
        else:
            self._set_qseq()
            self.update_func(l, f, i)
            
    def _iter_changed(self):
        self._set_qseq()
        self.update_func(*self.get_values())

    def reset(self, init=False):
        if self.PC.cur_expt:
            self.set_values(self.PC.min_len, 0, 0)
            self.wFragLen.configure(entryfield_validate =
                                    self._validate_len, 
                                    entryfield_entry_state = 'normal',
                                    entryfield_modifiedcommand =
                                    self._len_changed,
                                    downarrow_state = 'normal',
                                    uparrow_state = 'normal')
            self.wCurFrag.configure(entryfield_validate =
                                    self._validate_frag, 
                                    entryfield_entry_state = 'normal',
                                    entryfield_modifiedcommand =
                                    self._frag_changed,
                                    downarrow_state = 'normal',
                                    uparrow_state = 'normal')
            self.wCurIter.configure(entryfield_validate =
                                    self._validate_iter, 
                                    entryfield_entry_state = 'normal',
                                    entryfield_modifiedcommand =
                                    self._iter_changed,
                                    downarrow_state = 'normal',
                                    uparrow_state = 'normal')
            self._set_qseq()
            if not init:
                self.update_func(*self.get_values())
                self.set_func()
        else:
            for w in [self.wFragLen, self.wCurFrag, self.wCurIter]:  
                w.configure(entryfield_modifiedcommand = None,
                            entryfield_validate = None,
                            entryfield_entry_state = 'disabled',
                            downarrow_state = 'disabled',
                            uparrow_state = 'disabled')
                w.setvalue(' ')
            self.wSeqLabel.configure(text='')  
            if not init:
                self.unset_func()
