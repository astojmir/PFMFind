import Tkinter
from tkMessageBox import askquestion
import tkMessageBox
import Pmw
import os
import sys
import cPickle
import gzip
from ShortFrags.Expt import fragexpt
from ShortFrags.Expt import index
from ShortFrags.Expt import DirichletInfo
from matrix_opts import *
from query_input import *
from search_action import *
from search_params import *
from status import *
from matrix_view import *
from index_view import *
from hits_view import *
from full_view import *
from index_creation import CreateIndexDialog
from Console import Console
from progress_bar import ProgressBar
import threading
from Bio.SubsMat import MatrixInfo

THREADS = 5

class MainGui(Tkinter.Frame):
    def __init__(self, parent=None, dict={}):
        Tkinter.Frame.__init__(self,parent)
        self.pack(fill='both', expand = 1)
        self.parent=parent

        # Get screen resolution and set (fixed) heights and widths
        # Others are variable
        self.width = self.winfo_screenwidth()
        self.height = self.winfo_screenheight()
        self.params_width = 290
        self.msg_bar_height = 20
        
        # State variables 
        self.saved = True
        self.state = {}
        self.jobs = {}

        # ***** Basic Layout *****
        # Menu on the top
        self.balloon = Pmw.Balloon(self)
        self.menu_bar = Pmw.MenuBar(self, balloon = self.balloon, hotkeys=True)
        self.menu_bar.pack(side='top', anchor='nw', fill='x', expand=1)

        # Message Bar and progress bar on the bottom
        self.wInfoPart = Tkinter.Frame(self)
        self.wInfoPart.pack(side='bottom', anchor='s', fill='x', expand = 1)

        self.messageBar = Pmw.MessageBar(self.wInfoPart,
                                         entry_relief='groove',
                                         entry_font = Pmw.logicalfont('Fixed'))
        self.messageBar.pack(side='left', anchor='w', fill='x', expand = 1)
        self.messageBar.message('state', '')
        self.wProgressBar = ProgressBar(self.wInfoPart, width=int(self.width*0.15),
                                        height=self.msg_bar_height)
        self.wProgressBar.pack(side='right', anchor='e')
        self.wProgressBar.pack_forget()

        # Toolbox on the left
        self.workPart = Tkinter.Frame(self)
        self.workPart.pack(fill='both', expand=1)

        self.paramsPart = Tkinter.Frame(self.workPart, width=self.params_width)
        self.paramsPart.pack(side='left', anchor='nw', fill='y')
        self.paramsPart.pack_propagate(0)

        self.toolbox_filler = Tkinter.Frame(self.paramsPart, bd=3,
                                            relief='raised')
        self.toolbox_filler.pack(anchor='n',fill='both', expand=1)

        # Views on the right
        self.mainPart = Tkinter.Frame(self.workPart, bd=3, relief='raised')
        self.mainPart.pack(side='left', anchor='nw', fill='both', expand = 1)


        # ***** End of Basic Layout *****


        # **** Viewing Bar ****
        self.view_bar = Pmw.RadioSelect(self.mainPart,
                                        labelpos = None,
                                        command = self.view_button,
                                        )
        self.view_bar.pack(side='bottom', fill='x', expand = 1)            

        # **** End of Viewing Bar ****


        # **** Toolbox - create but don't pack (hidden) ****
        # Status
        self.statusPart = Tkinter.Frame(self.paramsPart, bd=3, relief='raised')
        self.IndxStatus = StatusShow(self.statusPart,
                                     ufunc = self.update_counters)
        self.IndxStatus.pack()

        # Search Arguments
        self.SrchArgs = Tkinter.Frame(self.paramsPart, bd=3, relief='raised')

        # Search Parameters
        self.MatParams = MatrixOpts(self.SrchArgs)
        self.sm = self.MatParams.add_case(MatrixInfo.available_matrices,
                                          ['None', 'Quasi', 'Avg', 'Max'],
                                          'blosum62',
                                          'None',
                                          lambda tag: self.update_matrix())
        self.pm = self.MatParams.add_case(['PSSM', 'subst mat'],
                                          update_func = lambda tag: self.update_matrix(),
                                          scale_default=1.0,
                                          weights = ['None', 'Henikoff'],
                                          regularisers = DirichletInfo.NAMES,
                                          reg_default = 'recode3.20comp')
        self.MatParams.pack()

        self.CutoffParams = SrchCutoff(self.SrchArgs)
        self.CutoffParams.pack(anchor='w')

        # Actions
        self.no_srch = Tkinter.IntVar()
        self.NoSearch = Tkinter.Checkbutton(self.SrchArgs, text="No Search",
                                            variable=self.no_srch)
        self.NoSearch.pack(anchor='w',pady=5)
        self.no_srch.set(0)
        
        self.Actions = ActionChooser(self.SrchArgs,
                                     set=self.update_jobs,
                                     run=self.run_jobs,
                                     reset=self.run_jobs)
        self.Actions.pack()

        # Filters
        self.filterArgs = Tkinter.Frame(self.paramsPart, bd=3,
                                        relief='raised')
        # **** End of Toolbox ****


        # **** Main displays ****
        self.console_view = Console(self.mainPart, dict=dict)
        self.index_view = IndexView(self.mainPart)
        self.matrix_view = MatrixView(self.mainPart)
        self.wHitsView = HitsView(self.mainPart)
        self.wFullView = FullView(self.mainPart,
                                  msg_func = lambda text:
                                  self.messageBar.message('state',text),
                                  set_func = self.IndxStatus.set_values)

        # **** End of Main displays ****


        # **** Top Menu ****
        self.menu_bar.addmenu('File', 'Close this window or exit')
        self.menu_bar.addmenuitem('File', 'command', 'Create new project',
                            command = self.new, label = 'New')
        self.menu_bar.addmenuitem('File', 'command', 'Open saved project',
                            command = self.open, label = 'Open...')
        self.menu_bar.addmenuitem('File', 'separator')
        self.menu_bar.addmenuitem('File', 'command', 'Save project',
                            command = self.save, label = 'Save')
        self.menu_bar.addmenuitem('File', 'command', 'Save project',
                            command = self.save_as, label = 'Save As...')
        self.menu_bar.addmenuitem('File', 'separator')
        self.menu_bar.addmenuitem('File', 'command', 'Close this window',
                            command = self.close, label = 'Close')
        self.menu_bar.addmenuitem('File', 'separator')
        self.menu_bar.addmenuitem('File', 'command', 'Create Index',
                            command = self.create_index, label = 'Create Index')
        self.menu_bar.addmenuitem('File', 'separator')
        self.menu_bar.addmenuitem('File', 'command', 'Exit the application',
                            command = self.quit, label = 'Quit')

##         self.menu_bar.addmenu('Edit', 'Undo/Redo Searches')
##         self.menu_bar.addmenuitem('Edit', 'command', 'Undo last search',
##                             command = None, label = 'Undo')
##         self.menu_bar.addmenuitem('Edit', 'command', 'Redo previous search',
##                             command = None, label = 'Redo')

        self.menu_bar.addmenu('View', 'View Results')

        self.menu_bar.addmenu('Options', 'Set search and filter options')
        self.menu_bar.addmenuitem('Options', 'command', 'Settings',
                                  command = None, label = 'Settings...')

##         self.menu_bar.addmenu('Tools', 'Search and filter')
        
        self.menu_bar.addmenu('Help', 'User manuals', name = 'help')
        self.menu_bar.addmenuitem('Help', 'command', 'About this application',
                                  command = None, label = 'About...')

        # **** End of Top Menu ****

        # Configure the balloon to displays its status messages in the
        # message bar.
        self.balloon.configure(statuscommand = self.messageBar.helpmessage)

        # Add Views to both top menu and view bar
        # Need to have sorted list of names as well
        self.view_names = ['Console',
                           'Index',
                           'Matrix',
                           'Hits',
                           'Full',
                           ]
        self.views = {'Console': self.console_view,
                      'Index': self.index_view,
                      'Matrix': self.matrix_view,
                      'Hits': self.wHitsView,
                      'Full': self.wFullView,
                      }
        self.viewVar = Tkinter.StringVar()
        for text in self.view_names: 
            self.menu_bar.addmenuitem('View', 'radiobutton', 'View ' + text,
                                command = lambda : self.view_button(None),
                                label = text, variable=self.viewVar,
                                value=text)
        for text in self.view_names:
            self.view_bar.add(text)

        # Height and width of the work part
        # - expands and shrinks according to the screen size
        self.work_width = self.width - self.params_width
        self.menu_bar.update()
        self.wInfoPart.update()
        self.work_height = self.height -\
                           self.msg_bar_height -\
                           self.menu_bar.winfo_height() - 100

        for v in self.views.values():
            v.set_size(self.work_height, self.work_width)
        

        # Disable all
        self.current_view = self.index_view
        self.disable_search()
        self.disable_views()
        
        # Redirect messages to console
        sys.stdout = self.console_view.stdout
        sys.stderr = self.console_view.stderr


    def new(self):
        seqin = QuerySeqInput(self)
        res = seqin.activate(geometry = 'centerscreenalways')
        if res == {}:
            return

        # Initialise query sequence
        # res['clusters_file'],res['index']  
        self.FE = fragexpt.FullExpt()
        self.FE.set_params(res['sequence'],
                           res['description'],
                           res['range'],
                           res['lengths'])

        self.save_fullname = res['save_name']
        self.saved = False

        self.enable_views()

        try:
            self.messageBar.message('state','Loading Index...')
            self.messageBar.update_idletasks()
            self.FE.load_index(res['index'])
            self.enable_search()
        except Exception, inst: 
            showerror('Load Error',  inst.__str__(), parent=self)

        try:
            self.messageBar.message('state','Loading Descriptions...')
            self.messageBar.update_idletasks()
            self.FE.load_descriptions(res['clusters_file'])
        except Exception, inst: 
            showerror('Load Error',  inst.__str__(), parent=self)
        
        self.messageBar.message('state','')

    def open(self):
        self.save_fullname = \
        askopenfilename(defaultextension='.ftb',
                        filetypes=[('Fragment Experiment','.ftb')],
                        parent=self,
                        title='Open Experiment',
                        initialdir=os.getcwd()
                        )
        if not len(self.save_fullname):
            return

        try:
            self.messageBar.message('state','Loading Experiment...') 
            self.messageBar.update_idletasks()
            self.FE = fragexpt.load_FE(self.save_fullname)
            self.enable_views()
            self.saved = True
        except Exception, inst:
            showerror('Open Error',  inst.__str__(), parent=self)
            return

        try:
            if self.FE.index_filename != None:
                self.messageBar.message('state','Loading Index...')
                self.messageBar.update_idletasks()
                self.FE.load_index(self.FE.index_filename)
                self.enable_search()
            elif self.FE.fasta_filename != None:
                self.messageBar.message('state','Loading FASTA filename...')
                self.messageBar.update_idletasks()
                self.FE.load_fastadb(self.FE.fasta_filename)
        except Exception, inst: 
            showerror('Open Error',  inst.__str__(), parent=self)


        try:
            if self.FE.desc_filename != None:
                self.messageBar.message('state','Loading Descriptions...')
                self.messageBar.update_idletasks()
                self.FE.load_descriptions(self.FE.desc_filename)
        except Exception, inst: 
            showerror('Open Error',  inst.__str__(), parent=self)
            
        self.messageBar.message('state','')
        
    def save(self):
        self.messageBar.message('state', 'Saving Experiment...')
        try:
            self.FE.save(self.save_fullname)
        except (RuntimeError, IOError), inst:
            showerror('Save Error', inst.__str__(), parent=self)
        self.messageBar.message('state','')
        self.saved = True

    def save_as(self):
        head, tail = os.path.split(self.save_fullname)
        path = asksaveasfilename(defaultextension='.ftb',
                                 filetypes=[('Fragment Experiment','.ftb')],
                                 parent = self,
                                 title = 'Choose Save Filename',
                                 initialdir = head,
                                 )
        if not len(path):
            return
        self.save_fullname = path
        self.save()

    def close(self):
        msg = 'Do you want to save changes?'
        # Have we saved?
        if not self.saved:
            res= tkMessageBox.Message(parent=self,
                                      icon= "warning",
                                      type= "yesnocancel",
                                      title= 'Fragment Toolbox',
                                      message= msg,
                                      ).show()
            if res == 'cancel':
                return
            elif res == 'yes':
                self.save()                
        # Cleanup all members and hide toolbox
        self.disable_search()
        self.disable_views()

    def create_index(self):
        dg = CreateIndexDialog(self)

        while 1:
            res = dg.activate(geometry = 'centerscreenalways')
            if res == {}:
                break
            try:
                self.messageBar.message('state', 'Creating Index...')
                self.messageBar.update_idletasks()
                self.update()
                I = index.FSIndex(res['fasta_name'], res['pttn'], res['use_sa'])
                self.messageBar.message('state', 'Saving Index...')
                self.messageBar.update_idletasks()
                I.save(res['index_name'])
                self.messageBar.message('state', '')
                break
            except Exception, inst:
                showerror('Index Creation Error', inst.__str__(), parent=self)


    def quit(self):
        self.close()
        self.parent.destroy()

    def view_button(self, button):
        if button == None:
            button = self.viewVar.get()
            self.view_bar.setvalue(button)
        else:
            self.viewVar.set(button)

        new_view = self.views[button]
        if new_view != self.current_view:
            self.current_view.pack_forget()
            new_view.pack(anchor='nw', fill='both', expand=1)
            new_view.pack_propagate(0)
            self.current_view = new_view
        self.work_width = self.width - self.params_width
        self.menu_bar.update_idletasks()
        self.wInfoPart.update_idletasks()
        self.work_height = self.height -\
                           self.msg_bar_height -\
                           self.menu_bar.winfo_height() - 100
        new_view.set_size(self.work_height, self.work_width)
        new_view.reset(self.state)

    def update_counters(self):
        v = self.IndxStatus.get_values()
        
        self.state['iter'] = v[2]
        self.state['coords'] = v[0:2]
        c = v[0:2]

        case = 1
        if v[2] == 0: case = 0

        # Update all
        if c in self.jobs and self.jobs[c]['iter'] == v[2]:
            self.MatParams.set_values(case, self.jobs[c]['matrix'],
                                      self.jobs[c]['conv_type'])
            self.CutoffParams.set_value(self.jobs[c]['cutoff'])
            # TO DO: update views
        else:
            self.MatParams.set_values(case)
            self.CutoffParams.restore_defaults()

        cur_view = self.viewVar.get()
        if cur_view == 'Hits' or cur_view == 'Full':
            self.current_view.reset(self.state)
            
    def update_matrix(self):
        res = self.MatParams.get_data()
        self.state.update(res)
        if res['matrix'] == 'subst mat':
            self.MatParams.wReg.component('menubutton').configure(state='disabled')
        else:
            self.MatParams.wReg.component('menubutton').configure(state='normal')
            
        if self.viewVar.get() == 'Matrix':
            self.current_view.reset(self.state)

    def update_jobs(self):
        act = self.Actions.get_values()
        coords = self.state['coords']
        no_srch = self.no_srch.get()

        conf = []
        remn = []
        if act[0] == 'Current':
            # Only one coordnate
            if coords in self.jobs:
                conf = [coords]
            else:
                remn = [coords]
        else:
            # Get the list of lengths
            if act[1] == 'Fragments':
                lengths = [coords[0]]
            else:
                lengths = range(*self.FE.length_range)

            fl = len(self.FE.rqseq)+ 1

            # Generate the list of conflicts and OK coords
            if not no_srch:
                remn = [(l,f) for l in lengths for f in range(fl-l)
                        if (l,f) not in self.jobs]
            if act[0] == 'All' or no_srch:
                conf = [(l,f) for l in lengths for f in range(fl-l)
                        if (l,f) in self.jobs]
                if not no_srch and len(conf) and \
                askquestion('Conflicts', 'Change previously scheduled searches?',
                         parent=self) == 'no':
                    conf = []
        
        # Now either set or delete jobs
        if no_srch:
            for key in conf:
                del(self.jobs[key])
            self.update_counters()
        else:
            val = {}
            val['cutoff'] = self.CutoffParams.get_value()
            val['matrix'] = self.state['matrix']
            val['conv_type'] = self.state['conv_type']
            val['iter'] = self.state['iter']
            val['scale'] = self.state['scale']
            val['weights'] = self.state['weight']
            val['reg'] = self.state['reg']
            for key in conf + remn:
                self.jobs[key] = val

    def run_jobs(self):
        # Check for conflicts
        for k,v in self.jobs.items():
            iters = self.FE.get_iters(*k)
            if len(iters) > v['iter']:
                if askquestion('Conflicts', 'Overwrite previous searches?',
                               parent=self) == 'no':
                    return
                else:
                    break

        self.wProgressBar.updateProgress(0, len(self.jobs))
        self.wProgressBar.pack()
        self.messageBar.message('state', 'Searching...')
        self.FE.serial_search(self.jobs, self.wProgressBar.updateProgress)
        self.messageBar.message('state', '')
        self.wProgressBar.pack_forget()
        self.IndxStatus.set_values(*self.state['coords'])
        self.saved = False
                                
    def enable_views(self):
        self.toolbox_filler.pack_forget()
        self.statusPart.pack(side='top', fill='both')
        self.filterArgs.pack(anchor='n',fill='both', expand=1)

        # Enable menus
        # File menu
        mb = self.menu_bar.component('File' + '-menu')
        for text in  ['Save', 'Save As...', 'Close']: 
            i = mb.index(text)
            mb.entryconfigure(i, state='normal')
        for text in  ['New', 'Open...']: 
            i = mb.index(text)
            mb.entryconfigure(i, state='disabled')

        # View menu
        mb = self.menu_bar.component('View' + '-menu')
        for text in ['Hits', 'Full']: 
            i = self.view_bar.index(text)
            self.view_bar.button(i).config(state='normal')
            i = mb.index(text)
            mb.entryconfigure(i, state='normal')

        self.state['FE'] = self.FE

        # Update counters
        self.IndxStatus.reset_params(self.FE.length_range,
                                     len(self.FE.rqseq),
                                     lambda l,f: len(self.FE.get_iters(l,f)),
                                     self.update_counters,
                                     self.FE.rqseq)                
        self.update_counters()
        
    def enable_search(self):
        if self.FE.index_filename == None:
            return

        # Show toolbox
        self.filterArgs.pack_forget()
        self.SrchArgs.pack(anchor='n', fill='x')
        self.filterArgs.pack(anchor='n',fill='both', expand=1)
        
        # View menu
        mb = self.menu_bar.component('View' + '-menu')
        for text in ['Index', 'Matrix']: 
            i = self.view_bar.index(text)
            self.view_bar.button(i).config(state='normal')
            i = mb.index(text)
            mb.entryconfigure(i, state='normal')
            
    def disable_search(self):
        self.view_bar.invoke('Console')
        self.SrchArgs.pack_forget()
        # View menu
        mb = self.menu_bar.component('View' + '-menu')
        for text in ['Index', 'Matrix']: 
            i = self.view_bar.index(text)
            self.view_bar.button(i).config(state='disabled')
            i = mb.index(text)
            mb.entryconfigure(i, state='disabled')

    def disable_views(self):
        self.view_bar.invoke('Console')

        # Delete variables
        self.state = {}
        self.jobs = {}
        self.saved = True
        vars = ['FE',
                'save_fullname',
                ]
        for v in vars:
            if v in self.__dict__:
                del(self.__dict__[v])

        # Hide toolbox
        self.statusPart.pack_forget()
        self.filterArgs.pack_forget()
        self.toolbox_filler.pack(anchor='n',fill='both', expand=1)

        # Disable menus
        # File menu
        mb = self.menu_bar.component('File' + '-menu')
        for text in  ['Save', 'Save As...', 'Close']: 
            i = mb.index(text)
            mb.entryconfigure(i, state='disabled')
        for text in  ['New', 'Open...']: 
            i = mb.index(text)
            mb.entryconfigure(i, state='normal')
        # View menu
        mb = self.menu_bar.component('View' + '-menu')
        for text in ['Hits', 'Full']: 
            i = self.view_bar.index(text)
            self.view_bar.button(i).config(state='disabled')
            i = mb.index(text)
            mb.entryconfigure(i, state='disabled')
        
        # Clear views
        self.wFullView.clear()
