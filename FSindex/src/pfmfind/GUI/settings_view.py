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


import Tkinter, tkFont, Pmw, os, os.path, sys
from tkFileDialog import askopenfilename, askdirectory,\
     asksaveasfilename
from tkMessageBox import showerror, showinfo, Message,\
     askquestion
from pfmfind.GUI.view import View
from pfmfind.GUI.optionmenu import OptionMenu

_FIXED_RANGE = True # Allow only lengths from 6 to 20


class SettingsView(Tkinter.Frame, View):
    def __init__(self, parent, PFMF_client, update_func):
        self.parent = parent
        Tkinter.Frame.__init__(self, parent)
        self.pack_propagate(0)

        self.PFMF_client = PFMF_client
        self.update_func = update_func

        # ************************************************************
        # ******* Database Settings **********************************
        # ************************************************************

        self.wDbGrp = Pmw.Group(self, tag_text='Database Settings')
        self.wDbGrp.pack(anchor='nw', fill='x', padx=5, pady=7)

        w = Tkinter.Frame(self.wDbGrp.interior())
        w.pack(anchor='w')

        Tkinter.Label(w, text="PostgreSQL Driver:").grid(row=0,\
            column=0, padx=5, pady=5, sticky='w')
        self.wDriverMenu = OptionMenu(w, items=['psycopg2', 'pgdb'],
                                      menubutton_width=20,
                                      menu_font=self.ffont,
                                      menubutton_font=self.ffont)
        self.wDriverMenu.grid(row=0, column=1, padx=5, pady=5, sticky='w')

        Tkinter.Label(w, text="Database:").grid(row=0, column=2, padx=5,
                                              pady=5, sticky='w')
        self.wDbEntry = Pmw.EntryField(w, entry_width = 25,
                                       entry_font=self.ffont)
        self.wDbEntry.grid(row=0, column=3, padx=5, pady=5,
                           sticky='w')

        Tkinter.Label(w, text="Host:").grid(row=2, column=0, padx=5,
                                            pady=5, sticky='w')
        self.wHostEntry = Pmw.EntryField(w, entry_width = 30,
                                            entry_font=self.ffont)
        self.wHostEntry.grid(row=2, column=1, padx=5, pady=5,
                             sticky='w')

        Tkinter.Label(w, text="Port:").grid(row=2, column=2, padx=5,
                                            pady=5, sticky='w')
        self.wPortEntry = Pmw.EntryField(w, entry_width =10,
                                            entry_font=self.ffont)
        self.wPortEntry.grid(row=2, column=3, padx=5, pady=5,
                             sticky='w')

        Tkinter.Label(w, text="User:").grid(row=4, column=0, padx=5,
                                            pady=5, sticky='w')
        self.wUserEntry = Pmw.EntryField(w, entry_width =25,
                                            entry_font=self.ffont)
        self.wUserEntry.grid(row=4, column=1, padx=5, pady=5,
                             sticky='w')

        Tkinter.Label(w, text="Password:").grid(row=4, column=2,
                                                padx=5, pady=5,
                                                sticky='w')
        self.wPswdEntry = Pmw.EntryField(w, entry_width =25,
                                         entry_show="*",
                                         entry_font=self.ffont)
        self.wPswdEntry.grid(row=4, column=3, padx=5, pady=5,
                             sticky='w')

        Tkinter.Label(w, text="Dataset schema:").grid(row=6, column=0,
                                                      padx=5, pady=5,
                                                      sticky='w')

        self.wDSEntry = Pmw.EntryField(w, entry_width =25,
                                       entry_font=self.ffont)
        self.wDSEntry.grid(row=6, column=1, padx=5, pady=5,
                           sticky='w')
        Tkinter.Label(w, text="PFMFind schema:").grid(row=6, column=2,
                                                      padx=5, pady=5,
                                                      sticky='w')
        self.wPSEntry = Pmw.EntryField(w, entry_width =25,
                                       entry_font=self.ffont)
        self.wPSEntry.grid(row=6, column=3, padx=5, pady=5,
                           sticky='w')

        self.wDbButFrm = Tkinter.Frame(w)
        self.wDbButFrm.grid(row=8, column=0, columnspan=2,
                            padx=5, pady=10, sticky='we')

        self.wDbClrBut = Tkinter.Button(self.wDbButFrm,
                                        text='Reset',
                                        command=self._reset_db)
        self.wDbClrBut.pack(side='left', padx=10)
        self.wDbConBut = Tkinter.Button(self.wDbButFrm,
                                        text='Connect',
                                        command=self._connect_db)
        self.wDbConBut.pack(side='left', padx=10)
        self.wDbDisBut = Tkinter.Button(self.wDbButFrm,
                                        text='Disconnect',
                                        command=self._disconnect_db)
        self.wDbDisBut.pack(side='left', padx=10)


        # ************************************************************
        # ******* Index Settings *************************************
        # ************************************************************

        self.wIndexGrp = Pmw.Group(self, tag_text='Index Settings')
        self.wIndexGrp.pack(anchor='nw', padx=5, pady=7, fill='x')

        w = Tkinter.Frame(self.wIndexGrp.interior())
        w.pack(anchor='w')

        Tkinter.Label(w, text="Host:").grid(row=0,
                                            column=0, padx=5,
                                            pady=5, sticky='w')
        self.wIxHostEntry = Pmw.EntryField(w, entry_width=30,
                                           entry_font=self.ffont)
        self.wIxHostEntry.grid(row=0, column=1, padx=5, pady=5,
                               sticky='w')

        Tkinter.Label(w, text="Port:").grid(row=0, column=2, padx=5,
                                            pady=5, sticky='w')
        self.wIxPortEntry = Pmw.EntryField(w, entry_width=10,
                                           entry_font=self.ffont)
        self.wIxPortEntry.grid(row=0, column=3, padx=5, pady=5,
                               sticky='w')


        self.wIxButFrm = Tkinter.Frame(w)
        self.wIxButFrm.grid(row=2, column=0, columnspan=2,
                            padx=5, pady=10, sticky='we')

        self.wIxClrBut = Tkinter.Button(self.wIxButFrm,
                                        text='Reset',
                                        command=self._reset_ix)
        self.wIxClrBut.pack(side='left', padx=10)
        self.wIxConBut = Tkinter.Button(self.wIxButFrm,
                                        text='Connect',
                                        command=self._connect_ix)
        self.wIxConBut.pack(side='left', padx=10)
        self.wIxDisBut = Tkinter.Button(self.wIxButFrm,
                                        text='Disconnect',
                                        command=self._disconnect_ix)
        self.wIxDisBut.pack(side='left', padx=10)



        # ************************************************************
        # ******* Plugin Settings ************************************
        # ************************************************************

        self.wPluginGrp = Pmw.Group(self, tag_text='Plugin Settings')
        self.wPluginGrp.pack(anchor='nw', padx=5, pady=7, fill='x',
                             expand=1)

        w = Tkinter.Frame(self.wPluginGrp.interior())
        w.pack(anchor='w')

        Tkinter.Label(w, text='Custom Plugin Path:').grid(row=0,\
            column=0, padx=5, pady=5, sticky='w')

        self.wPlgPathEntry = Pmw.EntryField(w, entry_width=60,
                                            entry_font=self.ffont)
        self.wPlgPathEntry.grid(row=0, column=1, padx=5, pady=5,
                                sticky='w')

        self.wPlgPathBut = Tkinter.Button(w, text='Choose...',\
            width=5, command=self._set_plugin_path)
        self.wPlgPathBut.grid(row=0, column=3, padx=5, pady=5,
                              sticky='w')

        self.wPlButFrm = Tkinter.Frame(w)
        self.wPlButFrm.grid(row=2, column=0, columnspan=2,
                            padx=5, pady=10, sticky='we')

        self.wPlClrBut = Tkinter.Button(self.wPlButFrm,
                                        text='Reset',
                                        command=self._reset_pl)
        self.wPlClrBut.pack(side='left', padx=10)
        self.wPlSetBut = Tkinter.Button(self.wPlButFrm,
                                        text='Set',
                                        command=self._set_pl)
        self.wPlSetBut.pack(side='left', padx=5)

        # ************************************************************
        # ******* Load/Save buttons **********************************
        # ************************************************************

        self.wLoad = Tkinter.Button(self, text='Load Settings...',
                                    command=self._load)
        self.wLoad.pack(side='left', anchor='nw', padx=5, pady=5)

        self.wSave = Tkinter.Button(self, text='Save Current Settings...',
                                    command=self._save)
        self.wSave.pack(side='right', anchor='ne', padx=5, pady=5)


        # Update the form
        self._db_form = {'db': self.wDbEntry, 'host': self.wHostEntry,
                         'port': self.wPortEntry, 'user': self.wUserEntry,
                         'password': self.wPswdEntry}


        self.update()

    def _reset_db(self):
        if self.PFMF_client.driver:
            self.wDriverMenu.setvalue(self.PFMF_client.driver)
        else:
            self.wDriverMenu.invoke(0)

        for k, v in self._db_form.iteritems():
            if k in self.PFMF_client.dbargs:
                self._db_form[k].setvalue(self.PFMF_client.dbargs[k])
            else:
                self._db_form[k].setvalue('')

        if self.PFMF_client.db_schema:
            self.wDSEntry.setvalue(self.PFMF_client.db_schema)
        else:
            self.wDSEntry.setvalue('')

        if self.PFMF_client.PFMF_schema:
            self.wPSEntry.setvalue(self.PFMF_client.PFMF_schema)
        else:
            self.wPSEntry.setvalue('')

    def _connect_db(self):
        dbargs = dict()

        s = self.wPortEntry.getvalue().strip()
        if len(s):
            try:
                port = int(s)
            except ValueError:
                showerror('Input Error', 'Port must be a number.',
                          parent=self.parent)
                return
            dbargs['port'] = s

        dbargs['driver'] = self.wDriverMenu.getvalue()

        s = self.wDbEntry.getvalue().strip()
        if len(s):
            dbargs['db'] = s

        s = self.wHostEntry.getvalue().strip()
        if len(s):
            dbargs['host'] = s

        s = self.wUserEntry.getvalue().strip()
        if len(s):
            dbargs['user'] = s
        s = self.wPswdEntry.getvalue().strip()
        if len(s):
            dbargs['password'] = s

        try:
            self.PFMF_client.open(**dbargs)
        except:
            showerror('Connection Error', 'Could not connect to' \
                      ' PostgreSQL database.',
                      parent=self.parent)
            return

        schemata = dict()
        s = self.wDSEntry.getvalue().strip()
        if len(s):
            schemata['db_schema'] = s
        s = self.wPSEntry.getvalue().strip()
        if len(s):
            schemata['PFMF_schema'] = s
        try:
            self.PFMF_client.set_schema(**schemata)
        except:
            showerror('Connection Error', 'Could not initialise the' \
                      ' given\n or default database schemata.',
                      parent=self.parent)
            self.PFMF_client.close()
            return

        self.update()

    def _disconnect_db(self):
        self.PFMF_client.close()
        self.update()


    def _reset_ix(self):
        if self.PFMF_client.host:
            self.wIxHostEntry.setvalue(self.PFMF_client.host)
        else:
            self.wIxHostEntry.setvalue('')
        if self.PFMF_client.port:
            self.wIxPortEntry.setvalue(str(self.PFMF_client.port))
        else:
            self.wIxPortEntry.setvalue('')

    def _connect_ix(self):
        host = self.wIxHostEntry.getvalue().strip()
        port = self.wIxPortEntry.getvalue().strip()

        try:
            port = int(port)
        except ValueError:
            showerror('Input Error', 'Port must be a number.',
                      parent=self.parent)
            return

        if not self.PFMF_client.attach(host, port):
            showerror('Connection Error', 'Could not connect to' \
                      ' FSIndex server\n at %s:%s.' % (host, port),
                      parent=self.parent)
            return
        self.update()

    def _disconnect_ix(self):
        self.PFMF_client.detach()
        self.update()

    def _set_plugin_path(self):
        path = askdirectory(mustexist=1,
                            parent = self.parent,
                            title = 'Choose Plugin Directory',
                            initialdir=os.getcwd(),
                            )
        if path == ():
            return
        self.wPlgPathEntry.setvalue(path)

    def _reset_pl(self):
        if self.PFMF_client.plugin_dir:
            self.wPlgPathEntry.setvalue(self.PFMF_client.plugin_dir)
        else:
            self.wPlgPathEntry.setvalue('')

    def _set_pl(self):
        path = self.wPlgPathEntry.getvalue()
        if not os.path.isdir(path):
            showerror('Input Error', 'Invalid plugin path.',
                      parent=self.parent)
            return
        self.PFMF_client.init_plugins(path)

    def _load(self):
        path = askopenfilename(defaultextension='.xml',
                               filetypes=[('XML File','.xml')],
                               parent = self.parent,
                               title = 'Choose Settings File',
                               initialdir=os.getcwd())
        if not len(path): return
        fp = file(path, 'r')
        self.PFMF_client.read_config(fp)
        fp.close()
        self.update()

    def _save(self):
        path = asksaveasfilename(defaultextension='.xml',
                               filetypes=[('XML File','.xml')],
                               parent = self.parent,
                               title = 'Choose Settings File',
                               initialdir=os.getcwd())
        if path == (): return
        fp = file(path, 'w')
        self.PFMF_client.write_config(fp)
        fp.close()


    def update(self):
        self._reset_db()
        self._reset_ix()
        self._reset_pl()

        if self.PFMF_client.conn:
            self.wDbConBut.configure(state='disabled')
            self.wDbDisBut.configure(state='normal')
            self.wLoad.configure(state='disabled')

            if self.PFMF_client.host and self.PFMF_client.port:
                self.wDbDisBut.configure(state='disabled')
                self.wIxConBut.configure(state='disabled')
                self.wIxDisBut.configure(state='normal')
            else:
                self.wIxConBut.configure(state='normal')
                self.wIxDisBut.configure(state='disabled')

        else:
            self.wDbConBut.configure(state='normal')
            self.wDbDisBut.configure(state='disabled')
            self.wLoad.configure(state='normal')
            self.wIxConBut.configure(state='disabled')
            self.wIxDisBut.configure(state='disabled')

        self.update_func()
