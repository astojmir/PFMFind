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
import os.path
from tkMessageBox import askquestion
from tkMessageBox import showerror, showinfo, Message
from tkFileDialog import askopenfilename
from pfmfind.GUI.view import View
from pfmfind.GUI.progress_dialog import ProgressDialog
from pfmfind.search.db import db

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

        self.wBatch = Tkinter.Button(self, text='Upload From File',
                                     command=self._load_batch)
        self.wBatch.pack(side='left', anchor='nw', padx=5, pady=5)

        self.wUpdate = Tkinter.Button(self, text='Update',
                                      command=self._update_database)
        self.wUpdate.pack(side='right', anchor='ne', padx=5, pady=5)

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

    def _validate_input(self, exp_data, new_exp=False, show_errors=True):

        names = set(s[1] for s in self.exp_list)
        a = 'ACDEFGHIKLMNPQRSTVWY'

        name, desc, qseq, qdesc, min_len, max_len = exp_data
        errmsg = None

        # Non-empty experiment description
        if len(desc) == 0:
            errmsg = 'Must have experiment description.'

        if new_exp and errmsg is None:
            # Name must be non-empty
            if len(name) == 0:
                errmsg = 'No experiment name given.'

            # Name may have no more than 40 chars.
            elif len(name) > 40:
                errmsg = 'Experiment name too long.',

            # Name must be unique
            elif name in names:
                errmsg = 'Experiment name must be unique.'

            # Sequence non-empty and valid
            elif len(qseq) == 0:
                errmsg = 'No sequence input.'

            elif len(qseq.translate('#'*256, a)) > 0:
                errmsg = 'Sequence contains invalid letters.'

            # Sequence description non-empty
            elif len(qdesc) == 0:
                errmsg = 'Must have sequence description.'

        if errmsg is not None:
            if show_errors:
                showerror('Input Error', errmsg, parent=self.parent)
            return False
        return True

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

    def _load_batch(self):

        path = askopenfilename(defaultextension='.fas',
                               filetypes=[('FASTA File','.fas'),
                                          ('FASTA File','.fasta')],
                               parent = self.parent,
                               title = 'Choose FASTA-formatted File',
                               initialdir=os.getcwd())
        if not len(path): return
        root_name = os.path.basename(path)

        try:
            batchdb = db(path)
        except:
            showerror('Batch Loading Error',
                      'Could not process selected FASTA file.',
                      parent=self.parent)
            return

        # pb = ProgressDialog(self.parent, text="Uploading query sequences",
        #                     max=batchdb.no_seq+1)
        # pb.activate(geometry = 'centerscreenalways')
        names = set(s[1] for s in self.exp_list)
        num_invalid = 0
        for i in xrange(batchdb.no_seq):
            name = root_name[:33] + '-%6.6d' % i
            if name in names:
                num_invalid += 1
                continue
            qseq = batchdb.get_seq(i)
            qdesc = batchdb.get_def(i)
            desc = name + ' : ' + qdesc
            exp_data = (name, desc, qseq, qdesc, 6, 20)
            if self._validate_input(exp_data, new_exp=True, show_errors=False):
                eid = self.PFMF_client.create_experiment(*exp_data)
                self.reset_experiment_list(set_gui=False)
            else:
                num_invalid += 1
        # pb.deactivate()
        if num_invalid:
            showerror('Batch Loading Error',
                      '%d sequences from the selected'
                      ' FASTA file could not be loaded.' % num_invalid,
                      parent=self.parent)
        self.reset_experiment_list(set_gui=True)

    def reset_experiment_list(self, set_gui=True):
        self.exp_list = self.PFMF_client.get_experiment_names()
        names = [s[1] for s in self.exp_list]

        if set_gui:
            self.wExpCombo.setlist(names)

    def get_form_data(self):

        name = str(self.wExpCombo.component('entryfield').getvalue())
        name = name.rstrip('\n ')
        desc = str(self.wDescText.getvalue().rstrip('\n '))
        qseq = str(self.wSeqText.getvalue())
        qseq = qseq.replace(' ', '').replace('\n', '').strip().upper()
        qdesc = str(self.wDtxt.getvalue().rstrip('\n '))
        return (name, desc, qseq, qdesc, 6, 20)

    def reset(self, state=None):
        self.reset_experiment_list()
