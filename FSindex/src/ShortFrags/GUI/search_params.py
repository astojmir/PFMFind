import Tkinter
import Pmw

class SrchCutoff(Tkinter.Frame):
    def __init__(self, parent=None, rng_default=25, kNN_default=100,
                 rel_default=0.1):
        Tkinter.Frame.__init__(self,parent)
        Tkinter.Label(self, text='Search Cutoff:').grid(sticky='w', pady=5)

        vld1 = {'validator':'numeric', 'min':0}
        vld2 = {'validator':'real', 'min':0.0, 'max':0.1}

        self.srch_types = [
            ['range','Range',rng_default,int, None, None, vld1],
            ['kNN', 'Nearest\nNeigbours', kNN_default, int, None, None, vld1],
            ['rel', 'Relative\nRange', rel_default, float, None, None, vld2],
            ]
        self.srch_type=Tkinter.IntVar()

        for i in range(len(self.srch_types)):
            t = self.srch_types[i]
            t[4] = Tkinter.Radiobutton(self, text=t[1],
                                    variable=self.srch_type,
                                    value=i)
            t[4].grid(row=i+1,column=0, sticky='w')
            t[5] = Pmw.EntryField(self, entry_width = 8,
                                  validate = t[6],
                                  value=t[2])
            t[5].grid(row=i+1, column=1, sticky='w')
            
        self.srch_type.set(0)
        
    def restore_defaults(self):
        self.srch_type.set(0)
        for t in self.srch_types:
            t[5].setvalue(str(t[2]))
            
    def get_value(self):
        i = self.srch_type.get()
        conv = self.srch_types[i][3]
        v = conv(self.srch_types[i][5].getvalue())
        return (i, v)
    
    def set_value(self, val):
        self.restore_defaults()
        self.srch_type.set(val[0])
        self.srch_types[val[0]][5].setvalue(str(val[1]))

    def disable(self):
        for t in self.srch_types:
            t[4].config(state='disabled')
            t[4].update()
            t[5].entry.config(state='disabled')
            t[5].update()

    def enable(self):
        for t in self.srch_types:
            t[4].config(state='normal')
            t[4].update()
            t[5].entry.config(state='normal')
            t[5].update()
            
