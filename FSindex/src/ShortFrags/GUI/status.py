import Tkinter
import Pmw

                    
class StatusShow(Tkinter.Frame):
    def __init__(self, parent=None, lrange=(5,15),
                 seqlen = 99, iter_func = lambda l,f: 0,
                 ufunc=lambda : 0, rqseq=""):
        Tkinter.Frame.__init__(self,parent)

        self.seqlen = seqlen
        self.max_iter_func = iter_func
        self.update_func = ufunc
        self.rqseq = rqseq
        
        self.ranges = [lrange,
                       (0, seqlen - lrange[0]+1),
                       (0, iter_func(lrange[0],0)+1)]
        
        self.wFrame = Tkinter.Frame(self)
        self.wFrame.pack(anchor='w', pady=2)
        self.wLabel = Tkinter.Label(self.wFrame, text="Query Frag:")
        self.wLabel.grid(row=0, column=0, sticky='w')
        self.wSeqLabel = Tkinter.Label(self.wFrame, text="",
                                       font=Pmw.logicalfont('Fixed'))
        self.wSeqLabel.grid(row=0, column=1, sticky='w')
                                            
        self.frag_len = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Fragment Length',
                                    entry_width = 6,
                                    entryfield_value = lrange[0],
                                    entryfield_validate =
                                    lambda t: self._validate(t,0),
                                    entryfield_modifiedcommand =
                                    self._len_changed)
        self.frag_len.pack(anchor='w', pady=2)

        self.cur_frag = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Current Fragment',
                                    entry_width = 6,
                                    entryfield_value = 0,
                                    entryfield_modifiedcommand =
                                    self._frag_changed,
                                    entryfield_validate =
                                    lambda t: self._validate(t,1),
                                    )
        self.cur_frag.pack(anchor='w', pady=2)
        self.cur_iter = Pmw.Counter(self,
                                    labelpos = 'w',
                                    label_text = 'Current Iteration',
                                    entry_width = 6,
                                    entryfield_value = 0,
                                    entryfield_modifiedcommand =
                                    self._iter_changed,
                                    entryfield_validate =
                                    self._validate_iter,
                                    )
        self.cur_iter.pack(anchor='w', pady=2)

    def reset_params(self, lrange, seqlen, iter_func, ufunc, rqseq):
        self.seqlen = seqlen
        self.max_iter_func = iter_func
        self.update_func = ufunc
        
        self.ranges = [lrange,
                       (0, seqlen - lrange[0]+1),
                       (0, iter_func(lrange[0],0)+1)]
        self.rqseq = rqseq
        self.reset_values()


    def reset_values(self):
        self.set_values(self.ranges[0][0],
                        0, self.ranges[2][1]-1)

    def set_values(self, len, frag, itr=None):
        self.frag_len.setvalue(str(len))
        self.cur_frag.setvalue(str(frag))
        if itr == None:
            itr = self.max_iter_func(len,frag)
        self.cur_iter.setvalue(str(itr))
        
    def get_values(self):
        l = int(self.frag_len.getvalue())
        f = int(self.cur_frag.getvalue())
        i = int(self.cur_iter.getvalue())
        return (l,f,i)

    def _validate(self, text, which):
        try:
            v = int(text)
        except ValueError:
            return Pmw.ERROR

        a = self.ranges[which][0]
        b = self.ranges[which][1]
        if v >= a and v < b:
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _set_qseq(self):
        l = int(self.frag_len.getvalue())
        f = int(self.cur_frag.getvalue())
        s = self.rqseq[f:f+l]
        self.wSeqLabel.configure(text=s)  

    def _validate_iter(self, text):
        try:
            v = int(text)
        except ValueError:
            return Pmw.ERROR
        
        l = int(self.frag_len.getvalue())
        f = int(self.cur_frag.getvalue())
        if v >= 0 and v <= self.max_iter_func(l,f):
            return Pmw.OK
        else:
            return Pmw.ERROR

    def _len_changed(self):
        l = int(self.frag_len.getvalue())
        self.ranges[1] = (0, self.seqlen - l+1)
        f = int(self.cur_frag.getvalue())
        if f >= self.ranges[1][1]:
            f = self.ranges[1][1] - 1
            self.cur_frag.setvalue(str(f))
        else:
            self._set_qseq()
            self.update_func()

    def _frag_changed(self):
        l = int(self.frag_len.getvalue())
        f = int(self.cur_frag.getvalue())
        self.ranges[2] = (0, self.max_iter_func(l,f)+1) 
        i = int(self.cur_iter.getvalue())
        if i >= self.ranges[2][1]:
            i = self.ranges[2][1] - 1
            self.cur_iter.setvalue(str(i))
        else:
            self._set_qseq()
            self.update_func()
            
    def _iter_changed(self):
        self._set_qseq()
        self.update_func()
