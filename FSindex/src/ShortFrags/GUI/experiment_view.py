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


import Tkinter, tkFont, Pmw, threading
from tkMessageBox import askquestion
from tkMessageBox import showerror, showinfo, Message
from ShortFrags.GUI.view import View


_FIXED_RANGE = True # Allow only lengths from 6 to 20


class ExperimentView(Tkinter.Frame, View):

    def __init__(self, parent, PFMF_client, set_exp_func):
        self.parent = parent
        Tkinter.Frame.__init__(self, parent)
        self.pack_propagate(0)

        self.PFMF_client = PFMF_client
        self.thread_lock = threading.Lock()
        self.db_access_lock = threading.Lock()
        self.current_thread = 0

        self.selected_exp_id = None
        self.inserted_exp_id = None
        self.selected_exp_data = None
        self.exp_list = list()

        self.set_exp_func = set_exp_func
        self.new_exp = False

        # ************************************************************
        # ******* Experiment List ************************************
        # ************************************************************

        self.wExpList = Pmw.Group(self, tag_text='Experiments')
        self.wExpList.pack(anchor='nw', fill='x',
                           padx=5, pady=7)

        self.wButtonBar = Pmw.ButtonBox(self.wExpList.interior(),
                                        orient='vertical')
        self.wButtonBar.grid(row=0, column=2, padx=5, pady=5, sticky='w')

        for text, callback in [('New', self._new_experiment),
                               ('Delete', self._delete_experiment),
                               ('Set', self._set_experiment),
                                ]:
            self.wButtonBar.add(text, command=callback)

        self.wExpCombo = Pmw.ComboBox(self.wExpList.interior(),
                                      selectioncommand =
                                      self._change_experiment,
                                      scrolledlist_items = [],
                                      listbox_width = 45,
                                      listbox_height = 6,
                                      scrolledlist_vscrollmode='static',
                                      scrolledlist_hscrollmode='none',
                                      dropdown=0,
                                      history=0,
                                      listbox_font=self.ffont,
                                      entry_font=self.ffont)
        self.wExpCombo.grid(row=0, column=1, padx=5, pady=10, sticky='we')
        self.wExpList.interior().columnconfigure(1, weight=2)
        self.wDescText = Pmw.ScrolledText(self.wExpList.interior(),
                                          labelpos='w',
                                          label_text='Description:',
                                          labelmargin=5,
                                          vscrollmode='static',
                                          hscrollmode='none',
                                          text_height=2,
                                          text_wrap='word',
                                          text_font=self.ffont,
                                          text_padx = 4,
                                          text_pady = 4,
                                          text_state = 'disabled',
                                          text_cursor = 'left_ptr',
                                          )
        self.wDescText.grid(row=1, columnspan=2, sticky='we', padx=5, pady=7)


        # ************************************************************
        # ******* Experiment details *********************************
        # ************************************************************

        self.wExpOpts = Pmw.Group(self, tag_text='Query Sequence')
        self.wExpOpts.pack(anchor='nw', padx=5, pady=7, fill='x', expand=1)
        self.wExpOpts.interior().columnconfigure(1, weight=2)

        w = self.wExpOpts.interior()
        # Description
        Tkinter.Label(w, text='Description:').grid(row=1, column=0)
        self.wDtxt = Pmw.ScrolledText(w, text_width=65, text_height=2,
                                     hscrollmode='none', vscrollmode='static',
                                     text_wrap='word')
        self.wDtxt.grid(row=1, column=1, columnspan=3, padx=5, pady=5, sticky='we')

        # Sequence
        Tkinter.Label(w, text='Sequence:').grid(row=2, column=0)
        self.Sbut = Tkinter.Button(w, text='Format', width=5,
                                    command = self._format_seq)
        self.Sbut.grid(sticky='n', row=4, column=0)
        self.wSeqText = Pmw.ScrolledText(w, text_width=65, text_height=8,
                                         hscrollmode='none',
                                         vscrollmode='static',
                                         text_wrap='word')
        self.wSeqText.grid(row=2, column=1, columnspan=3,
                           rowspan=3,padx=5, pady=5, sticky='we')


        # ************************************************************
        # ******* Update buttons *************************************
        # ************************************************************
        self.wReset = Tkinter.Button(self, text='Reset',
                                     command=self._reset_form)
        self.wReset.pack(side='left', anchor='nw', padx=5, pady=5)

        self.wUpdate = Tkinter.Button(self, text='Update',
                                      command=self._update_database)
        self.wUpdate.pack(side='right', anchor='ne', padx=5, pady=5)

##         if self.PFMF_client.conn:
##             self.reset_experiment_list()

        # Range - fixed to 6--20
        if _FIXED_RANGE: return
        # ********* NOT USED FOR NOW *********************************
        self.vValid = {'validator':'numeric', 'min':0}
        RF = Tkinter.Frame(w)
        RF.grid(sticky='w', row=5, column=1, pady=10)

        Tkinter.Label(w, text='Length Range:').grid(row=5, column=0)
        self.L0 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'From',
                                 labelmargin=5, entry_width = 5,
                                 validate = self.vValid, value=7)
        self.L0.grid(row=5, column=1, padx=12)
        self.L1 = Pmw.EntryField(RF, labelpos = 'w', label_text = 'To',
                                 labelmargin=5, entry_width = 5,
                                 validate = self.vValid, value=13)
        self.L1.grid(row=5, column=2)
        # ************************************************************



    def _change_experiment(self, text=None):

        if len(self.exp_list) == 0:
            return

        index = int(self.wExpCombo.curselection()[0])
        exp_id = self.exp_list[index][0]

        if exp_id == self.selected_exp_id: return

        self.selected_exp_id = exp_id
        self.inserted_exp_id = None
        self.selected_exp_data = None
        self._clear_data()

        thr = threading.Thread(target=self._db_get_data_thread,
                               args=(exp_id,))
        thr.start()
        self._insert_data(exp_id)


    def _db_get_data_thread(self, exp_id):
        self.db_access_lock.acquire()
        if exp_id == self.selected_exp_id:
            exp_data = self.PFMF_client.get_experiment_data(exp_id)
        else:
            return
        if exp_id == self.selected_exp_id:
            self.selected_exp_data = exp_data
        self.db_access_lock.release()

    def _insert_data(self, exp_id):

        if exp_id != self.selected_exp_id or \
           exp_id == self.inserted_exp_id :
            return

        if self.selected_exp_data:
            name, desc, qseq, qdesc, min_len, max_len =\
                  self.selected_exp_data

            self.wDescText.configure(text_state='normal')
            self.wDescText.setvalue(desc)
            self.wDtxt.configure(text_state='normal')
            self.wDtxt.setvalue(qdesc)
            self.wDtxt.configure(text_state='disabled')
            self.wSeqText.configure(text_state='normal')
            self.wSeqText.setvalue(qseq)
            self._format_seq()
            self.wSeqText.configure(text_state='disabled')
            if _FIXED_RANGE: return
            # ********* NOT USED FOR NOW *****************************
            self.L0.configure(entry_state='normal')
            self.L0.setvalue(str(min_len))
            self.L1.configure(entry_state='normal')
            self.L1.setvalue(str(max_len))
            # ********************************************************
        else:
            self.after(100, self._insert_data, exp_id)



    def _clear_data(self):
        self.wDescText.configure(text_state = 'normal')
        self.wDescText.setvalue("")
        self.wDescText.configure(text_state = 'disabled')
        self.wDtxt.configure(text_state = 'normal')
        self.wDtxt.setvalue("")
        self.wDtxt.configure(text_state = 'disabled')
        self.wSeqText.configure(text_state = 'normal')
        self.wSeqText.setvalue("")
        self.wSeqText.configure(text_state = 'disabled')
        if _FIXED_RANGE: return
        # ********* NOT USED FOR NOW *********************************
        self.L0.configure(entry_state='normal')
        self.L0.setvalue("")
        self.L0.configure(entry_state='disabled')
        self.L1.configure(entry_state='normal')
        self.L1.setvalue("")
        self.L1.configure(entry_state='disabled')
        # ************************************************************

    def _validate_input(self, exp_data, new_exp=False):
        name, desc, qseq, qdesc, min_len, max_len = exp_data
        if new_exp:
            # Name must be non-empty
            if len(name) == 0:
                showerror('Input Error',
                          'No experiment name given.',
                          parent=self.parent)
                return False

            # Name may have no more than 40 chars.
            if len(name) > 40:
                showerror('Input Error',
                          'Experiment name too long.',
                          parent=self.parent)
                return False

            # Name must be unique
            names = [s[1] for s in self.exp_list]
            if name in names:
                showerror('Input Error',
                          'Experiment name must be unique.',
                          parent=self.parent)
                return False

            # Sequence non-empty and valid
            a = 'ACDEFGHIKLMNPQRSTVWY'
            if len(qseq) == 0:
                showerror('Input Error',
                          'No sequence input.',
                          parent=self.parent)
                return False
            if len(qseq.translate('#'*256, a)) > 0:
                showerror('Input Error',
                          'Sequence contains invalid letters.',
                          parent=self.parent)
                return False

            # Sequence description non-empty
            if len(qdesc) == 0:
                showerror('Input Error',
                          'Must have sequence description.',
                          parent=self.parent)
                return False

        # Non-empty experiment description
        if len(desc) == 0:
            showerror('Input Error',
                      'Must have experiment description.',
                      parent=self.parent)
            return False

        if _FIXED_RANGE: return True
        # ********* NOT USED FOR NOW *********************************

        # Range: HI <= LOW + reasonable values
        if min_len == None or max_len == None:
            showerror('Input Error',
                      'At least one of length ranges is missing.',
                      parent=self.parent)
            return False

        if min_len > max_len:
            showerror('Input Error',
                      'Minimum length is greater than maximum length.',
                      parent=self.parent)
            return False

        if min_len < 6:
            showerror('Input Error',
                      'Minimum length is smaller than 6.',
                      parent=self.parent)
            return False

        if max_len > 20:
            showerror('Input Error',
                      'Maximum length is greater than 20.',
                      parent=self.parent)
            return False
        # ************************************************************

    def _new_experiment(self):
        self.new_exp=True
        self.selected_exp_id = None
        self.inserted_exp_id = None
        self.selected_exp_data = None

        self.wExpCombo.component('entryfield').clear()
        self.wDescText.configure(text_state = 'normal')
        self.wDescText.setvalue("")
        self.wDtxt.configure(text_state = 'normal')
        self.wDtxt.setvalue("")
        self.wSeqText.configure(text_state = 'normal')
        self.wSeqText.setvalue("")


    def _delete_experiment(self):
        # Cannot delete the current experiment
        if self.selected_exp_id == self.PFMF_client.cur_expt:
            showerror('Delete Error',
                      'Cannot delete the current experiment.',
                      parent=self.parent)
            return
        if askquestion('Delete Experiment',
                       "Are you sure you want to delete the selected" \
                       " experiment and all searches associated with" \
                       " it?", parent=self.parent) == 'no':
            return

        # This will block on database I/O
        self.PFMF_client.delete_experiment(self.selected_exp_id)
        self.reset_experiment_list()
        self._clear_data()
        self.wExpCombo.component('entryfield').clear()

    def _set_experiment(self):
        index = int(self.wExpCombo.curselection()[0])
        exp_id = self.exp_list[index][0]

        if exp_id == self.PFMF_client.cur_expt:
            self.PFMF_client.reset_current_experiment()
            self._root().title("PFMFind")
        else:
            self.PFMF_client.set_current_experiment(exp_id)
            self._root().title("PFMFind - %s" %\
                self.exp_list[index][1])
        self.set_exp_func()

    def _reset_form(self):
        if self.new_exp:
            self._new_experiment()
        else:
            exp_data = self.get_form_data()
            if exp_data != self.selected_exp_data:
                self._insert_data(self.selected_exp_id)


    def _update_database(self):
        exp_data = self.get_form_data()
        if self.new_exp:
            if not self._validate_input(exp_data, new_exp=True): return
            # This needs to block - too complicated to spawn
            # a thread
            eid = self.PFMF_client.create_experiment(*exp_data)
            self.reset_experiment_list()
            # Select the name
            self.wExpCombo.selectitem(exp_data[0])
            self.new_exp=False
            self._change_experiment()
        else:
            if exp_data[1] == self.selected_exp_data[1]: return
            if not self._validate_input(exp_data): return
            thr = threading.Thread(\
                target=self.PFMF_client.update_experiment,
                args=(self.selected_exp_id, exp_data[1]))
            thr.start()
            self.selected_exp_data = exp_data

    def _format_seq(self):
        s = self.get_form_data()[2]
        self.wSeqText.setvalue("")

        for i in range(10, len(s)+10,10):
            self.wSeqText.insert('end', s[i-10:i])
            self.wSeqText.insert('end', ' ')

    def reset_experiment_list(self):
        self.exp_list = self.PFMF_client.get_experiment_names()
        names = [s[1] for s in self.exp_list]

        self.wExpCombo.setlist(names)

    def get_form_data(self):
        name = self.wExpCombo.component('entryfield').getvalue()
        name = name.rstrip('\n ')
        desc = self.wDescText.getvalue().rstrip('\n ')

        qseq = str(self.wSeqText.getvalue())
        qseq = qseq.replace(' ', '').replace('\n', '').strip().upper()

        qdesc = self.wDtxt.getvalue().rstrip('\n ')

        if _FIXED_RANGE:
            min_len = 6
            max_len = 20
            return (name, desc, qseq, qdesc, min_len, max_len)

        # ********* NOT USED FOR NOW *********************************
        try:
            min_len = int(self.L0.getvalue())
        except ValueError:
            min_len = None
        try:
            max_len = int(self.L1.getvalue())
        except ValueError:
            max_len = None
        return (name, desc, qseq, qdesc, min_len, max_len)
        # ************************************************************

    def reset(self, state=None):
        self.reset_experiment_list()
