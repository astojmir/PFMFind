import Tkinter, tkMessageBox, Pmw
from tkMessageBox import askquestion
import os, sys

from ShortFrags.Expt import fragexpt

from status import StatusShow
from index_view import IndexView
from search_view import SearchView
from hits_view import HitsView
from keyword_view import KeywordView
from experiment_view import ExperimentView
from settings_view import SettingsView
from Console import Console

_ffont = Pmw.logicalfont('Fixed')

_SETTINGS = 0
_EXPERIMENT = 1
_SEARCH = 2
_RESULTS = 3
_INDEX = 4

from view import View
class DummyView(Tkinter.Frame, View):
    def __init__(self, parent, PFMF_client):
        Tkinter.Frame.__init__(self, parent)

class PFMFindGUI(Tkinter.Frame):
    def __init__(self, parent=None, config_file=None, globals_dict=None):
        Tkinter.Frame.__init__(self, parent)
        self.pack(fill='both', expand = 1)
        self.parent=parent


        # State variables 
        self.PFMF_client = fragexpt.PFMFindClient()
        self.state = (None, None, None)

        # Try to get the client running with a given config
        if config_file:
            fp = file(config_file, 'r')
            self.PFMF_client.read_config(fp)
            fp.close()
        
 
        # ************************************************************
        # ******* Menu on the top ************************************
        # ************************************************************

        self.wBalloon = Pmw.Balloon(self)
        self.wMenuBar = Pmw.MenuBar(self, balloon=self.wBalloon,
                                    hotkeys=True) 
        self.wMenuBar.pack(side='top', anchor='nw', fill='x',
                           expand=1) 

        self.wMenuBar.addmenu('File', 'Exit')
        self.wMenuBar.addmenuitem('File', 'command',
                                  'Exit the application',
                                  command = self.quit, label = 'Quit')
        self.wMenuBar.addmenu('View', 'Views')

        # ************************************************************
        # ******* Status bar *****************************************
        # ************************************************************

        self.wStatus = StatusShow(self, self.PFMF_client,
                                  ufunc=self._update_state,
                                  set_func=self._set_experiment,
                                  unset_func=self._unset_experiment)
        self.wStatus.pack(side='top', anchor='nw', fill='x', expand=1,
                          padx=10, pady=2)


        # ************************************************************
        # ******* Notebook with views ********************************
        # ************************************************************

        self.wNote = Pmw.NoteBook(self,
                                  createcommand=self._create_view,
                                  raisecommand=self._raise_view)
        self.wNote.pack(fill = 'both', expand=1, padx=10, pady=5) 
        self.vView = Tkinter.StringVar()

        # ************************************************************
        # ******* Message bar on the bottom **************************
        # ************************************************************

        self.wMessageBar = Pmw.MessageBar(self,
                                          entry_relief='groove',
                                          entry_font = _ffont)
        self.wMessageBar.pack(side='top', anchor='w',
                              fill='x', expand=1)  
        self.wMessageBar.message('state', '')
        self.wBalloon.configure(statuscommand =
                                self.wMessageBar.helpmessage) 
        self.msg_bar_height = 20

        # ************************************************************
        # ******* Set widget sizes ***********************************
        # ************************************************************

        self.work_width = self.winfo_screenwidth()
        self.work_height = self.winfo_screenheight() -\
                           self.msg_bar_height -\
                           self.wMenuBar.winfo_height() -\
                           self.wStatus.winfo_height() - 160

        self.wNote.component('hull').configure(height=self.work_height,
                                               width=self.work_width)

        # ************************************************************
        # ******* Create views ***************************************
        # ************************************************************

        # views: name, type, class, additional args as dict
        # page_PFMF_client=self.PFMF_client is passed automatically 
        
        self.views = [('Settings', _SETTINGS, SettingsView,
                       {'page_update_func': self.update_tab_state}),
                      ('Experiment', _EXPERIMENT, ExperimentView,
                       {'page_set_exp_func': self.wStatus.reset}),
                      ('Index', _INDEX, IndexView, {}),
                      ('Search', _SEARCH, SearchView,
                       {'page_set_func': self.wStatus.set_values}),
                      ('Hits', _RESULTS, HitsView,
                       {'page_set_func': self.wStatus.set_values}), 
                      ('Keywords', _RESULTS, KeywordView,
                       {'page_set_func': self.wStatus.set_values,
                        'page_msg_func': self._msg_func}),
                      ]
        if globals_dict:
            self.views.append(('Console', _SETTINGS, Console,
                               {'page_dict': globals_dict}))

        self.wViews = {}
        dummy_page = self.wNote.add('Dummy')
        
        for vname, vtype, vclass, vkwargs in self.views:
            self.wViews[vname] = \
                (self.wNote.add(vname, page_pyclass=vclass, 
                                page_PFMF_client=self.PFMF_client,
                                tab_state='disabled',
                                **vkwargs), vtype)
            self.wMenuBar.addmenuitem('View', 'radiobutton', vname,
                                      command=self._activate_view,
                                      label=vname,
                                      variable=self.vView, 
                                      value=vname, state='disabled')
        self.wNote.delete(0)


        # ************************************************************
        # ******* Enable some views **********************************
        # ************************************************************
 
        self.wStatus.reset()
        self._enable_tabs(_SETTINGS)
        self.update_tab_state()

    def _create_view(self, name):
        if name not in self.wViews: return
        w = self.wViews[name][0]
        w.set_size(self.work_height, self.work_width)

    def _raise_view(self, name):
        if name not in self.wViews: return
        w = self.wViews[name][0]
        w.set_size(self.work_height, self.work_width)
        w.reset(self.state)        
        self.vView.set(name)
        
    def _msg_func(self, text):
        self.wMessageBar.message('state', text)

    def _update_state(self, length, fragment, iteration):
        self.state = (length, fragment, iteration)
        name = self.wNote.getcurselection()
        self.wViews[name][0].reset(self.state)

    def _enable_tabs(self, vtype):
        mb = self.wMenuBar.component('View' + '-menu')
        for name, (w, wvtype) in self.wViews.items():
            if wvtype == vtype:
                self.wNote.tab(name).configure(state='normal')
                i = mb.index(name)
                mb.entryconfigure(i, state='normal')

    def _disable_tabs(self, vtype):
        mb = self.wMenuBar.component('View' + '-menu')
        for name, (w, wvtype) in self.wViews.items():
            if wvtype == vtype:
                self.wNote.tab(name).configure(state='disabled')
                i = mb.index(name)
                mb.entryconfigure(i, state='disabled')
                
    def _activate_view(self):
        name = self.vView.get()
        self.wNote.selectpage(name)


    def _set_experiment(self):
        if self.PFMF_client.attached:
            self._enable_tabs(_SEARCH)
        else:
            self._disable_tabs(_SEARCH)
            
        self._enable_tabs(_RESULTS)

    def _unset_experiment(self):
        self._disable_tabs(_SEARCH)
        self._disable_tabs(_RESULTS)

    def update_tab_state(self):
        if self.PFMF_client.conn:
            self._enable_tabs(_EXPERIMENT)
            if self.PFMF_client.attached:
                self._enable_tabs(_INDEX)
                if self.PFMF_client.cur_expt:
                    self._enable_tabs(_SEARCH)
                else:
                    self._disable_tabs(_SEARCH)
            else:
                self._disable_tabs(_INDEX)
                self._disable_tabs(_SEARCH)
        else:
            self._disable_tabs(_EXPERIMENT)
            self._disable_tabs(_INDEX)
            self._disable_tabs(_SEARCH)
            self._disable_tabs(_RESULTS)

    def quit(self):
        self.parent.destroy()

            
