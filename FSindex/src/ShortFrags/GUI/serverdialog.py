import Tkinter
import Pmw
from string import upper, translate, strip, replace
from tkMessageBox import showerror, showinfo, Message
from string import upper, translate, strip, replace
from tkFileDialog import *
import os

class ServerDialog(Pmw.Dialog):
    def __init__(self, parent=None):

        Pmw.Dialog.__init__(self, parent, 
                            buttons = ('Accept', 'Cancel', 'Reset'),
                            defaultbutton = 'Accept',
                            title = 'Enter Server Location',
                            command = self.command)
        w = self.interior()
        self.parent = parent
       
        self.V0 = {'validator':'integer', 'min':0, 'max':65535}

        # Host
        host = Tkinter.Frame(w)
        self.host = Pmw.EntryField(w, labelpos = 'w', label_text = 'Hostname:',
                                   entry_width = 25, value='depot')
        self.host.pack(side='left', pady=5)
        self.port = Pmw.EntryField(w, labelpos = 'w', label_text = 'Port:',
                                   entry_width = 6, validate = self.V0, value=50007)
        self.port.pack(side='left', pady=5)


    def command(self, result):
        if result == 'Reset':
            self.reset()
            return
        elif result == 'Cancel':
            self.deactivate({})
            return        
        # result == 'Accept'
        res = {}

        res['host'] = self.host.getvalue()
        res['port'] = int(self.port.getvalue())
        self.deactivate(res)

    def reset(self):
        self.host.setvalue('depot')
        self.port.setvalue(50007)

if __name__ == '__main__':
    root = Pmw.initialise()
    root.title("Short Fragment Toolbox")
    ServerDialog()
    root.mainloop()

