import Tkinter
import Pmw
from string import upper, translate, strip, replace
from tkMessageBox import showerror, showinfo, Message
from string import upper, translate, strip, replace
from tkFileDialog import *
import os

class QuerySeqInput(Pmw.Dialog):
    def __init__(self, parent=None, ind_lst=[], 
                 ind_default=0):

        Pmw.Dialog.__init__(self, parent, 
                            buttons = ('Accept', 'Cancel', 'Reset'),
                            defaultbutton = 'Accept',
                            title = 'Setup index and query sequence',
                            command = self.command)
        w = self.interior()
        self.parent = parent
       
##         # Index
##         self.Ilst = ind_lst
##         self.Idef = ind_default
##         if len(ind_lst) > 0:
##             Tkinter.Label(w, text='Index:').grid(row=0, column=0)
##             self.Ivar = Tkinter.StringVar()
##             self.Ivar.set(ind_lst[ind_default])
##             self.Imnu = Pmw.OptionMenu(w,
##                                        menubutton_textvariable = self.Ivar,
##                                        items = self.Ilst,
##                                        menubutton_width = 25)
##             self.Imnu.grid(row=0, column=1, columnspan=2, sticky='w')

        # Description
        Tkinter.Label(w, text='Query\nDescription:').grid(row=1, column=0)
        self.Dtxt = Pmw.ScrolledText(w, text_width=65, text_height=2,
                                     hscrollmode='none', vscrollmode='static',
                                     text_wrap='word')
        self.Dtxt.grid(row=1, column=1, columnspan=3, padx=5, pady=5)

        # Sequence
        Tkinter.Label(w, text='Query\nSequence:').grid(row=2, column=0)
        self.Lbut = Tkinter.Button(w, text='Load', width=5,
                                   command = self.load, state='disabled')
        self.Lbut.grid(sticky='n', row=3, column=0)
        self.Sbut = Tkinter.Button(w, text='Format', width=5,
                                    command = self.format_seq)
        self.Sbut.grid(sticky='n', row=4, column=0)
        self.Stxt = Pmw.ScrolledText(w, text_width=65, text_height=8,
                                     hscrollmode='none', vscrollmode='static',
                                     text_wrap='none')
        self.Stxt.grid(row=2, column=1, columnspan=3,
                       rowspan=3,padx=5, pady=5)

        # Configure the tag to mark the unused sequence
        self.Stxt.tag_config("n", foreground="grey")

        # Range
        self.V0 = {'validator':'numeric', 'min':0}
        self.V1 = {'validator':'numeric', 'min':0}
        RF = Tkinter.Frame(w)
        RF.grid(sticky='w', row=5, column=1, columnspan=3)
        
        Tkinter.Label(w, text='Range:').grid(row=5, column=0)
        self.R0 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'From',
                                 entry_width = 5, validate = self.V0, value=0)
        self.R0.pack(side='left', pady=5)
        self.R1 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'To',
                                 entry_width = 5, validate = self.V1, value=0)
        self.R1.pack(side='left', pady=10)

        self.Rbut = Tkinter.Button(RF, text='Show', command = self.format_range)
        self.Rbut.pack(side='left', padx=10)

        # Length Range
        self.V2 = {'validator':'numeric', 'min':5, 'max':21}
##         LF = Tkinter.Frame(w)
##         LF.grid(sticky='w', row=5, column=3, padx=5)
        
        Tkinter.Label(RF, text='Fragment\nLengths:').pack(side='left', padx=10)
        self.L0 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'From',
                                 entry_width = 3, validate = self.V2, value=7)
        self.L0.pack(side='left', pady=5)
        self.L1 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'To',
                                 entry_width = 3, validate = self.V2, value=13)
        self.L1.pack(side='left', pady=5)

        # Index name
        self.index_name = None
        Tkinter.Label(w, text='Index File:').grid(row=6, column=0)
        self.Itxt = Tkinter.Text(w, width=55, height=1, wrap='none', state='disabled')
        self.Itxt.grid(row=6, column=1, columnspan=2, padx=5, pady=5)
        self.Ibut = Tkinter.Button(w, text='Choose...', width=5,
                                   command = self.set_index_file)
        self.Ibut.grid(row=6, column=3, padx=5, pady=5)

        # Clusters name
        self.clusters_file = None
        Tkinter.Label(w, text='Clusters File:').grid(row=7, column=0)
        self.Mtxt = Tkinter.Text(w, width=55, height=1, wrap='none', state='disabled')
        self.Mtxt.grid(row=7, column=1, columnspan=2, padx=5, pady=5)
        self.Mbut = Tkinter.Button(w, text='Choose...', width=5,
                                   command = self.set_clusters_file)
        self.Mbut.grid(row=7, column=3, padx=5, pady=5)
        

        # Save name
        self.save_name = None
        Tkinter.Label(w, text='Save Filename:').grid(row=8, column=0)
        self.Ftxt = Tkinter.Text(w, width=55, height=1, wrap='none', state='disabled')
        self.Ftxt.grid(row=8, column=1, columnspan=2, padx=5, pady=5)
        self.Fbut = Tkinter.Button(w, text='Choose...', width=5,
                                   command = self.set_save_file)
        self.Fbut.grid(row=8, column=3, padx=5, pady=5)



    def command(self, result):
        if result == 'Reset':
            self.reset()
            return
        elif result == 'Cancel':
            self.deactivate({})
            return        
        # result == 'Accept'
        res = {}

        # Get sequence and validate
        hull = self.component('hull')
        seq = self.get_seq()
        a = 'ACDEFGHIKLMNPQRSTVWY'
        if len(seq) == 0:
            showerror('Input Error', 'No sequence input.', parent=hull)
            return           
        if len(translate(seq,'#'*256,a)) > 0:
            showerror('Input Error', 'Invalid sequence.', parent=hull)
            return
        res['sequence'] = seq
        
        # Get lengths
        lerr = 0
        l0 = self.L0.getvalue()
        l1 = self.L1.getvalue()
        if l0 == '' or l1 == '':
            lerr = 1
        else:
            l0 = int(l0)
            l1 = int(l1)
            if (l0 < l1):
                res['lengths'] = (l0,l1)
            else:
                lerr = 1
        if lerr:
            showerror('Input Error', 'Invalid lengths.', parent=hull)
            return

##         # Get index
##         if len(self.Ilst) > 0:
##             res['index'] = self.Ivar.get()

        # Get description
        desc = strip(self.Dtxt.get())
        if len(desc) == 0:
            desc = 'Unknown sequence'
            showinfo('Description',
                     "Description set to: '%s'" % desc,
                     parent=hull)
        res['description'] = desc

        # Get range
        self.set_valid_range()
        r0 = int(self.R0.getvalue())
        r1 = int(self.R1.getvalue())
        res['range'] = (r0,r1)

        # Get Index Name
        if self.index_name == None:
            showerror('Input Error', 'Missing Index Name.', parent=hull)
            return
        res['index'] = self.index_name

        # Get Clusters File
        if self.index_name == None:
            self.clusters_file = ""
            #showerror('Input Error', 'Missing Clusters File.', parent=hull)
            #return
        res['clusters_file'] = self.clusters_file
        
        # Get Save Name
        if self.save_name == None:
            showerror('Input Error', 'Missing Save Filename.', parent=hull)
            return
        res['save_name'] = self.save_name

        # Return res
        self.deactivate(res)

    def get_seq(self):
        s = self.Stxt.get()
        s = replace(s,' ', '')
        s = replace(s,'\n', '')
        s = upper(strip(s))
        return s
                
    def set_valid_range(self):
        r0 = self.R0.getvalue()
        r1 = self.R1.getvalue()
        s = self.get_seq()
        l = len(s)
        if r0 == '':
            r0 = 0
            self.R0.setvalue(r0)
        if r1 == '':
            r1 = l
            self.R1.setvalue(r1)

        r0 = int(r0)
        r1 = int(r1)
        if r1 > l:
            r1 = l
            self.R1.setvalue(r1)
        if r0 < r1:
            return
        r0 = 0
        r1 = l
        self.R0.setvalue(r0)
        self.R1.setvalue(r1)

        
    def format_range(self):
        self.Stxt.tag_remove('n', '1.0', 'end')
        self.set_valid_range()
        r0 = int(self.R0.getvalue())
        r1 = int(self.R1.getvalue())
        c0 = r0%60
        c1 = r1%60
        t0 = '1.0';
        t1 = '%d.%d' % (r0/60+1, c0+c0/10)
        t2 = '%d.%d' % (r1/60+1, c1+c1/10)
        t3 = 'end'
        self.Stxt.tag_add('n', t0, t1)
        self.Stxt.tag_add('n', t2, t3)
            
    def format_seq(self):
        s0 = self.get_seq()
        s1=''
        for i in range(10,len(s0)+10,10):
            s1 = s1 + s0[i-10:i]
            if i % 70 == 0:
                s1 += '\n'
            else:
                s1 += ' '
        self.Stxt.setvalue(s1)
        self.format_range()
        
    def load(self):
        pass

    def reset(self):
        self.Stxt.setvalue('')
        self.Dtxt.setvalue('')
        self.R0.setvalue(0)
        self.R1.setvalue(0)
        self.L0.setvalue(7)
        self.L1.setvalue(13)
        if len(self.Ilst) > 0:
            self.Ivar.set(self.Ilst[self.Idef])
        self.index_name = None
        self.clusters_file = None
        self.save_name = None
        self.Itxt.configure(state='normal')
        self.Itxt.delete(1.0, 'end')
        self.Itxt.configure(state='disabled')
        self.Mtxt.configure(state='normal')
        self.Mtxt.delete(1.0, 'end')
        self.Mtxt.configure(state='disabled')
        self.Ftxt.configure(state='normal')
        self.Ftxt.delete(1.0, 'end')
        self.Ftxt.configure(state='disabled')

    def set_index_file(self):
        path = askopenfilename(defaultextension='.ix',
                               filetypes=[('Index','.ix')],
                               parent = self.interior(),
                               title = 'Choose Index File',
                               initialdir=os.getcwd()
                               )
        if path == ():
            return
        self.Itxt.configure(state='normal')
        self.Itxt.delete(1.0, 'end')
        self.Itxt.insert(1.0, path)
        self.Itxt.configure(state='disabled')
        self.index_name = path

    def set_clusters_file(self):
        path = askopenfilename(parent = self.interior(),
                               title = 'Choose Clusters File',
                               defaultextension='.cls',
                               filetypes=[('Clusters','.cls')],
                               initialdir=os.getcwd()
                               )
        if path == ():
            return
        self.Mtxt.configure(state='normal')
        self.Mtxt.delete(1.0, 'end')
        self.Mtxt.insert(1.0, path)
        self.Mtxt.configure(state='disabled')
        self.clusters_file = path

    def set_save_file(self):
        path = asksaveasfilename(defaultextension='.ftb',
                                 filetypes=[('Fragment Experiment','.ftb')],
                                 parent = self.interior(),
                                 title = 'Choose Save Filename',
                                 initialdir=os.getcwd()
                                 )
        if path == ():
            return
        self.Ftxt.configure(state='normal')
        self.Ftxt.delete(1.0, 'end')
        self.Ftxt.insert(1.0, path)
        self.Ftxt.configure(state='disabled')
        self.save_name = path


if __name__ == '__main__':
    root = Pmw.initialise()
    root.title("Short Fragment Toolbox")
    QuerySeqInput()
    root.mainloop()

