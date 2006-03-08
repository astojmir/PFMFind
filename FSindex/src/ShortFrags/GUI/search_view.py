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


import Tkinter, Pmw, string, threading
from tkMessageBox import askquestion, showwarning
from ShortFrags.GUI.view import View
from ShortFrags.GUI.ScrolledSeq import ScrolledSeq, col_f
from ShortFrags.Expt.SearchServer import RNG_SRCH, KNN_SRCH, REL_SRCH

_ffont = Pmw.logicalfont('Fixed')

def _empty_func():
    pass

class SearchView(Tkinter.Frame, View):
    
    _default_plugins = ['default_matrix', 'default_profile']
    _search_types = [[REL_SRCH, 'E-value', 1.0, float,
                      {'validator':'real', 'min':0.0, 'max':1e05}],
                     [RNG_SRCH, 'Range', 25, int,
                      {'validator':'numeric', 'min':0}],
                     [KNN_SRCH, 'Nearest\nNeigbours', 100, int,
                      {'validator':'numeric', 'min':0}],
                     ]

    def __init__(self, parent, PFMF_client,
                 set_func = lambda l,f,i : 0,
                 enable_func=_empty_func,
                 disable_func=_empty_func):
        Tkinter.Frame.__init__(self,parent)
        self.pack_propagate(0)

        self.params_width = 290


        self.PFMF_client = PFMF_client
        self.set_func = set_func
        self.enable_func = enable_func
        self.disable_func = disable_func
        
        self.state = (None, None, None)
        self.jobs = {}
        self.search_finished = False
        self.qseq = None
        self.thread_lock = threading.Lock()
        self.current_thread = 0
        self.matrix_text = None
        self.search_type_widgets = []
        self.first_search = None
        self.plugin_widgets = []
        self.plugins = []
        self.default_plugin = None
        self.current_plugin = None

        # ************************************************************
        # ******* Search Options on the left *************************
        # ************************************************************

        self.wOptions = Tkinter.Frame(self, width=self.params_width)
        self.wOptions.pack(side='left', anchor='nw', fill='y')
        self.wOptions.pack_propagate(0)

        # ***** Search Arguments *************************************
        
        self.wCutoffParams = Pmw.Group(self.wOptions,
                                     tag_text='Cutoff Options')
        self.wCutoffParams.pack(anchor='nw', fill='x', padx=5, pady=5)

        self.vCutoffParams = Tkinter.IntVar()
        parent = self.wCutoffParams.interior()
        for i, st in enumerate(self._search_types):
            rb = Tkinter.Radiobutton(parent,
                                     text=st[1],
                                     variable=self.vCutoffParams,
                                     value=i)
            rb.grid(row=i+1,column=0, sticky='w')
            ef = Pmw.EntryField(parent,
                                entry_width=8,
                                validate=st[4],
                                value=st[2])
            ef.grid(row=i+1, column=1, sticky='w')

            self.search_type_widgets.append((rb,ef))
        self.vCutoffParams.set(0)

        # ***** Plugin input *****************************************
        
        self.wPluginInput = Pmw.Group(self.wOptions,
                                      tag_text='Matrix Options')
        self.wPluginInput.pack(fill='x', padx=5, pady=5)
        self.wPluginLabel = Tkinter.Label(self.wPluginInput.interior(),
                                          text='Plugin')
        self.wPluginLabel.grid(row=0, column=0, sticky='w',
                               padx=5, pady=5) 
        self.wPluginChooser = \
            Pmw.OptionMenu(self.wPluginInput.interior(),
                           items = [],
                           menubutton_width = 20,
                           command = self.set_plugin_widgets)
        self.wPluginChooser.grid(row=0, column=1, sticky='w')

        # ***** Jobs *************************************************

        self.wJobs = Pmw.Group(self.wOptions, tag_text='Jobs')
        self.wJobs.pack(fill='x', padx=5, pady=5)

        self.wJobs1 = Pmw.RadioSelect(self.wJobs.interior(),
                                      buttontype = 'radiobutton',
                                      labelpos = 'w',
                                      command = self._set_l,
                                      orient='vertical',
                                      pady=0)
        for text in ['Current', 'All', 'Remaining']:
            self.wJobs1.add(text)
        self.wJobs1.grid(row=1,column=0)

        self.wJobs2 = Pmw.RadioSelect(self.wJobs.interior(),
                                      buttontype = 'radiobutton',
                                      labelpos = 'w',
                                      orient='vertical',
                                      pady=0)
        for text in ['Fragments', 'Lengths']:
            self.wJobs2.add(text)
        self.wJobs2.grid(row=1, column=1, sticky='sw')

        self.wJobs1.invoke('Current')
        self.wJobs2.invoke('Lengths')

        self.wJobsButtons = Pmw.ButtonBox(self.wJobs.interior())
        self.wJobsButtons.grid(row=2, column=0, columnspan=2,
                               pady=2)  

        self.wJobsButtons.add('Set', command = self._set_jobs)
        self.wJobsButtons.add('Unset', command = self._unset_jobs)
        self.wJobsButtons.add('Run', command = self._run_searches)
        self.wJobsButtons.alignbuttons()


        # ************************************************************
        # ******* Sequence and Matrix on the right *******************
        # ************************************************************

        self.wSeqMat =  Tkinter.Frame(self)
        self.wSeqMat.pack(side='left', anchor='nw',
                          fill='both', expand = 1) 

        # ***** Full sequence viewer *********************************

        self.wSeqText = ScrolledSeq(self.wSeqMat,
                                    self.PFMF_client.query_sequence,
                                    self.set_func)
        self.wSeqText.pack(anchor='nw', fill='x')
  
        self.wSeqText.tag_config('green yellow', background = 'green yellow')
        self.wSeqText.tag_config('green', background = 'green')
        self.wSeqText.tag_config('red', background = 'red')
        self.wSeqText.tag_config('orange', background = 'orange')
        self.wSeqText.tag_config('royal blue',
                                 background = 'royal blue')

        # ***** Matrix viewer ****************************************

        self.wMatrixText = Pmw.ScrolledText(self.wSeqMat,
                                            usehullsize = 1,
                                            text_wrap='none',
                                            text_font = _ffont,
                                            text_padx = 4,
                                            text_pady = 4,
                                            text_state = 'disabled',
                                            text_cursor = 'left_ptr',
                                            )
        self.wMatrixText.pack(anchor='nw', fill='both', expand = 1)

    def _set_l(self, tag):
        nb = self.wJobs2.index(Pmw.END)
        if tag == 'Current':
            for i in range(nb+1):
                b = self.wJobs2.button(i)
                b.config(state='disabled')
                b.update()
        else:
            for i in range(nb+1):
                b = self.wJobs2.button(i)
                b.config(state='normal')
                b.update()
        
    def _highlight_jobs(self):
        self.wSeqText.tag_remove('green yellow', '1.0', '3.0')
        self.wSeqText.tag_remove('green', '1.0', '3.0')
        self.wSeqText.tag_remove('royal blue', '1.0', '3.0')
        self.wSeqText.tag_remove('orange', '1.0', '3.0')
        this_length = [(l,f) for l,f in self.jobs.iterkeys()
                       if l == self.state[0]]
        other_lengths = [(l,f) for l,f in self.jobs.iterkeys()
                         if l != self.state[0]]

        for l,f in this_length:
            pos = "2.%d" % col_f(f)
            if f in self.PFMF_client.max_iters:
                max_iters = self.PFMF_client.max_iters[f].get(l, -1)
            else:
                max_iters=-1
            if self.jobs[(l,f)][0] == max_iters + 1:
                self.wSeqText.tag_add('green yellow', pos)
            elif self.jobs[(l,f)][0] > max_iters + 1:
                self.wSeqText.tag_add('royal blue', pos)
            else:
                self.wSeqText.tag_add('orange', pos)

        fdict = dict()
        for l,f in other_lengths:
            if f in self.PFMF_client.max_iters:
                max_iters = self.PFMF_client.max_iters[f].get(l, -1)
            else:
                max_iters=-1
            fdict[f] = 0
            idiff = self.jobs[(l,f)][0] - (max_iters + 1)
            if idiff != 0 and fdict[f] != 0:
                fdict[f] = min(fdict[f], idiff)
            else:
                fdict[f] = idiff
                
        for f in fdict.iterkeys():
            pos = "1.%d" % col_f(f)
            if fdict[f] == 0:
                self.wSeqText.tag_add('green yellow', pos)
            elif fdict[f] > 0:
                self.wSeqText.tag_add('royal blue', pos)
            else:
                self.wSeqText.tag_add('orange', pos)

        # ALSO YELLOW IF NO HITS FOR NEXT ITERATIONS (THROUGH MAX_ITER)
        # RED IF INVALID QUERY SEQUENCE
        # THINK OF COMBINATIONS

    def _set_jobs(self):

        # Get search_type and cutoff
        i = self.vCutoffParams.get()
        search_type = self._search_types[i][0]
        v = self.search_type_widgets[i][1].getvalue()
        cutoff = self._search_types[i][3](v)

        # Get plugin options
        plugin, plugin_args = self._get_plugin_options()

        # Get job type
        job_type1 = self.wJobs1.getvalue()
        job_type2 = self.wJobs2.getvalue()

        (length, fragment, iteration) = self.state

        # Now assign jobs
        conflicts = []
        remaining = []
        if job_type1 == 'Current':
            # Only one coordnate
            if (length, fragment) in self.jobs:
                conflicts = [(length, fragment)]
            else:
                remaining = [(length, fragment)]
        else:
            if job_type2 == 'Fragments':
                lengths = [length]
            else:
                lengths = range(self.PFMF_client.min_len,
                                self.PFMF_client.max_len)

            fl = len(self.qseq)+ 1

            remaining = [(l,f) for l in lengths
                         for f in range(len(self.qseq) + 1 - l)
                         if (l,f) not in self.jobs]
            
            if job_type1 == 'All':
                conflicts = [(l,f) for l in lengths
                             for f in range(len(self.qseq) + 1 - l)
                             if (l,f) in self.jobs]
                if len(conflicts) and \
                    askquestion('Conflicts',
                                'Change previously scheduled searches?',
                                parent=self) == 'no':
                    conflicts = []

        for key in conflicts + remaining:
            self.jobs[key] = (iteration,
                              search_type,
                              cutoff,
                              plugin,
                              plugin_args)

        self._highlight_jobs()

    def _unset_jobs(self):

        # Get job type
        job_type1 = self.wJobs1.getvalue()
        job_type2 = self.wJobs2.getvalue()

        (length, fragment, iteration) = self.state

        if job_type1 == 'Remaining':
            return
        elif job_type1 == 'Current':
            self.jobs.pop((length, fragment), None)
        else:
            if job_type2 == 'Fragments':
                lengths = [length]
            else:
                lengths = range(self.PFMF_client.min_len,
                                self.PFMF_client.max_len)
            for l,i in self.jobs.keys():
                if l in lengths:
                    del(self.jobs[(l,i)])

        self._highlight_jobs()

    def _show_matrix(self, text=None):
        plugin, plugin_args = self._get_plugin_options()

        self.wMatrixText.setvalue('Computing matrix...')

        self.current_thread += 1
        thread_id = self.current_thread
        self.matrix_text = None
        thr = threading.Thread(target=self._show_matrix_thread,
                               args=(thread_id,
                                     self.state,
                                     plugin, plugin_args))  
        thr.start()
        self._process_matrix(thr, thread_id)
        
    def _process_matrix(self, thread, thread_id):
        if thread_id != self.current_thread:
            return
        elif thread.isAlive():
            self.after(100, self._process_matrix,
                       thread, thread_id)
        elif self.matrix_text and self.matrix_text[0] == thread_id: 
            self.wMatrixText.setvalue(self.matrix_text[1])

    def _show_matrix_thread(self, thread_id, state,
                            plugin, plugin_args):

        (length, fragment, iteration) = state
        plugin_func = self.plugins[plugin][0]
        HL = self.PFMF_client.select_lif_search(length,
                                                iteration-1,
                                                fragment)
        
        # Check if futher computation is really needed
        if thread_id != self.current_thread:
            return

        PM, matrix_type, ctype = plugin_func(HL, *plugin_args)
        matrix_text = \
            self.PFMF_client.plugin_matrix_str(PM, matrix_type, ctype)
        print_info = self.plugins[plugin][1]
        if print_info:
            info_text = '\n\n' + print_info(HL, *plugin_args)
        else:
            info_text = ""

        if thread_id == self.current_thread:
            self.matrix_text = (thread_id, matrix_text + info_text) 


    def _run_searches(self):
        self.search_finished = False
        self.disable_func()
        for i in xrange(self.wJobsButtons.numbuttons()):
            self.wJobsButtons.button(i).configure(state='disabled')
        thr = threading.Thread(target=self._search_thread)
        thr.start()
        self._update_searches()

    def _search_thread(self):
        self.PFMF_client.run_search_jobs(self.jobs)
        self.search_finished = True

    def _update_searches(self):
        if not self.search_finished:
            self.after(100, self._update_searches)
            return

        self.enable_func()
        for i in xrange(self.wJobsButtons.numbuttons()):
            self.wJobsButtons.button(i).configure(state='normal')
        self._highlight_jobs()

    def set_size(self, height=400, width=300):

        self.configure(height=height, width=width)
        self.update_idletasks()

    def reset(self, state=(None, None, None)):

        if state == self.state: return
        
        self.state = state
        (length, fragment, iteration) = self.state

        if self.qseq != self.PFMF_client.query_sequence:
            self.qseq = self.PFMF_client.query_sequence
            self.wSeqText.reset(state, self.qseq)
        else:
            self.wSeqText.reset(state)
        
        self._highlight_jobs()
        first_search = (iteration == 0)
        if first_search != self.first_search:
            self.first_search = first_search
            if self.first_search:
                self.plugins = self.PFMF_client.start_plugins
                self.default_plugin = self._default_plugins[0]
            else:
                self.plugins = self.PFMF_client.iteration_plugins  
                self.default_plugin = self._default_plugins[1]
            self.wPluginChooser.setitems(self.plugins.keys(),
                                         self.default_plugin)
            self.wPluginChooser.invoke()
        else:
            self._show_matrix()
           
    def set_plugin_widgets(self, plugin):
        self.current_plugin = plugin
        
        # Destroy old widgets
        for w, l, default in self.plugin_widgets:
            w.destroy()
            l.destroy()
        self.plugin_widgets = []
        
        # Create new widgets
        arg_list = self.plugins[plugin][2]

        for i, (label, choice, default_choice) in enumerate(arg_list):

            l = Tkinter.Label(self.wPluginInput.interior(),
                              text=label)
            l.grid(row=i+1, column=0, sticky='w', padx=5, pady=5)            
            
            if isinstance(choice, list):
                width = min(21, max(map(len, choice))) 
                w = Pmw.OptionMenu(self.wPluginInput.interior(),
                                   items = choice,
                                   initialitem = default_choice,
                                   menubutton_width = width,
                                   command = self._show_matrix)
                
            elif isinstance(choice, str):
                w = Pmw.EntryField(self.wPluginInput.interior(),
                                   entry_width = 8,
                                   validate = choice,
                                   value = default_choice,
                                   modifiedcommand = self._show_matrix) 
            else:
                w = Tkinter.Label(self.wPluginInput.interior(),
                                  text='INVALID OPTION')

                self.plugin_widgets.append((None, default_choice)) 

            w.grid(row=i+1, column=1, sticky='w')
            self.plugin_widgets.append((w, l, default_choice)) 
        
        self._show_matrix()

    def _get_plugin_options(self):
        plugin = self.wPluginChooser.getvalue()
        arg_list = self.plugins[plugin][2]
        plugin_args = []
        for i, (label, choice, default_choice) in enumerate(arg_list): 
            w = self.plugin_widgets[i][0]                
            if isinstance(choice, list):
                plugin_args.append(w.getvalue())
            elif isinstance(choice, str):
                if not w.valid():
                    plugin_args.append(default_choice)
                else:
                    val = w.getvalue()
                    if choice == 'real':
                        val = float(val)
                    elif choice == 'numeric' or choice == 'integer':
                        val = int(val)
                    plugin_args.append(val)
            else:
                plugin_args.append(default_choice)
        return plugin, plugin_args 

