import Tkinter
import Pmw
from tkMessageBox import showerror, showinfo, Message
from tkFileDialog import *
from progress_dialog import ProgressDialog

class CreateIndexDialog(Pmw.Dialog):
    def __init__(self, parent=None):

        Pmw.Dialog.__init__(self, parent, 
                            buttons = ('Accept', 'Reset', 'Cancel'),
                            defaultbutton = 'Accept',
                            title = 'Create Index',
                            command = self.command)
        wInt = self.interior()
        self.wParent = parent
       
        # Label with instructions
        text = "Make sure you know what you are doing. Use # as the separation " + \
        "character.\nKeep to the standard alphabet and make sure alphabets " + \
        "at each position are consistent."
        
        self.wInstr = Tkinter.Label(wInt, text=text, justify='left')
        self.wInstr.grid(row=0, column=0, columnspan=2, sticky='wnes',pady=10, padx=5)

        # Length and use suffix array flag
        self.index_length = 0
        self.old_length = 0
        self.wIlen = Pmw.Counter(wInt,
                                 labelpos = 'w',
                                 label_text = 'Index Length',
                                 entry_width = 6,
                                 entryfield_validate = {'validator' : 'integer',
                                                        'min' : 4, 'max' : 20},
                                 entryfield_value = 4,
                                 entryfield_modifiedcommand =
                                 self._Ilen_changed)
        self.wIlen.grid(row=1, column=0, padx=5, pady=5, sticky='w')

        self.use_sa = Tkinter.IntVar()
        self.wSarr = Tkinter.Checkbutton(wInt, text="Generate suffix array",
                                         variable=self.use_sa)
        self.wSarr.grid(row=1, column=1, padx=5, pady=5, sticky='w')
        self.use_sa.set(0)

        # Fasta db
        self.fasta_name = None
        self.wF = Tkinter.Frame(wInt)

        Tkinter.Label(self.wF, text='FASTA Filename:').grid(column=0, sticky='w')
        self.wFtxt = Tkinter.Text(self.wF, width=82, height=1, wrap='none', state='disabled')
        self.wFtxt.grid(row=0, column=1, sticky='w') 
        self.wFbut = Tkinter.Button(self.wF, text='Choose...', width=5,
                                   command = self.set_fasta_file)
        self.wFbut.grid(row=0, column=2, sticky='w', padx=5, pady=5)

        # Save as
        self.index_name = None
        Tkinter.Label(self.wF, text='Save Index As:').grid(column=0, sticky='w')
        self.wItxt = Tkinter.Text(self.wF, width=82, height=1, wrap='none', state='disabled')
        self.wItxt.grid(row=1, column=1, sticky='w')
        self.wIbut = Tkinter.Button(self.wF, text='Choose...', width=5,
                                    command = self.set_index_file)
        self.wIbut.grid(row=1, column=2, sticky='w', padx=5, pady=5)

        self.wF.grid(row=2, column=0, columnspan=2, sticky='w')
        
        # Entry fields for partitions - no validation
        self.wPartitions = []
        self.default_partition = "STA#N#ILVM#KR#DEQ#W#FY#H#G#P#C"
        for i in range(20):
            w = Pmw.EntryField(wInt, labelpos = 'w',
                               label_text = 'Position #%2d' % (i+1),
                               value = self.default_partition,
                               entry_width = 40)
            self.wPartitions.append(w)  

        # Set the length and put 12 positions
        self.wIlen.setvalue(12)


    def _Ilen_changed(self):
        self.old_length = self.index_length
        self.index_length = int(self.wIlen.getvalue())
        self.put_entry_fields()

    def put_entry_fields(self):
        if (self.old_length < self.index_length):
            for i in range(self.old_length,
                           self.index_length):
                j = 4 + i // 2
                k = i % 2
                self.wPartitions[i].grid(row=j, column=k, padx=5, pady=5,sticky='w')
                self.wPartitions[i].setvalue(self.default_partition)
        elif (self.old_length > self.index_length):
            for i in range(self.index_length,
                           self.old_length):
                self.wPartitions[i].grid_forget()
            
    def set_fasta_file(self):
       path = askopenfilename(defaultextension='.fas',
                                 filetypes=[('FASTA database','.fas'),
                                            ('FASTA database','.aa')],
                                 parent = self.interior(),
                                 title = 'Choose FASTA Filename',
                                 )
       if path == ():
           return
       self.wFtxt.configure(state='normal')
       self.wFtxt.delete(1.0, 'end')
       self.wFtxt.insert(1.0, path)
       self.wFtxt.configure(state='disabled')
       self.fasta_name = path

    def set_index_file(self):
       path = asksaveasfilename(defaultextension='.ix',
                                 filetypes=[('Fragment Index','.ix')],
                                 parent = self.interior(),
                                 title = 'Choose Save Filename',
                                 )
       if path == ():
           return
       self.wItxt.configure(state='normal')
       self.wItxt.delete(1.0, 'end')
       self.wItxt.insert(1.0, path)
       self.wItxt.configure(state='disabled')
       self.index_name = path

    def reset(self):
        for wP in self.wPartitions:
            wP.setvalue(self.default_partition)

    def command(self, result):
        if result == 'Reset':
            self.reset()
            return
        elif result == 'Cancel':
            self.deactivate({})
            return        
        # result == 'Accept'

        # Check that filenames are entered
        hull = self.component('hull')
        if self.fasta_name == None:
            showerror('Input Error', 'Missing FASTA Filename.', parent=hull)
            return
        if self.index_name == None:
            showerror('Input Error', 'Missing Index Filename.', parent=hull)
            return

        pttn = []
        for i in range(self.index_length):
            pttn.append(self.wPartitions[i].getvalue())

        dict = {'pttn': pttn,
                'use_sa': self.use_sa.get(),
                'fasta_name': self.fasta_name,
                'index_name': self.index_name,
                }
        self.deactivate(dict)
                        
if __name__ == '__main__':
    root = Pmw.initialise()
    root.title("Short Fragment Toolbox")
    CreateIndexDialog()
    root.mainloop()
