import Tkinter
import Pmw

class ActionChooser(Tkinter.Frame):
    def __init__(self, parent=None, run=None, set=None, reset=None):
        Tkinter.Frame.__init__(self,parent)

        w = Pmw.Group(self, tag_text='Apply to')
        w.grid(row=1, column=0, columnspan=2, sticky='nwse')
        cw = w.interior()
        
        self.lradio = Pmw.RadioSelect(cw,
                                      buttontype = 'radiobutton',
                                      labelpos = 'w',
                                      command = self.set_l,
                                      orient='vertical',
                                      pady=0,
                                      )
        for text in ['Current', 'All', 'Remaining']:
            self.lradio.add(text)
        self.lradio.grid(row=1,column=0)

        self.rradio = Pmw.RadioSelect(cw,
                                      buttontype = 'radiobutton',
                                      labelpos = 'w',
                                      command = self.set_r,
                                      orient='vertical',
                                      pady=0)
        for text in ['Fragments', 'Lengths']:
            self.rradio.add(text)
        self.rradio.grid(row=1, column=1, sticky='sw')
        self.lradio.invoke('Current')
        self.rradio.invoke('Lengths')

        self.buttons = Pmw.ButtonBox(self)
        self.buttons.grid(row=2,column=0,columnspan=2,pady = 2)

        self.buttons.add('Set', command = set)
        self.buttons.add('Search', command = run)
        self.buttons.add('Reset', command = reset)


    def set_l(self, tag):
        self.which_var = tag
        if tag == 'Current':
            nb = self.rradio.index(Pmw.END)
            for i in range(nb+1):
                b = self.rradio.button(i)
                b.config(state='disabled')
                b.update()
        else:
            nb = self.rradio.index(Pmw.END)
            for i in range(nb+1):
                b = self.rradio.button(i)
                b.config(state='normal')
                b.update()

    def set_r(self, tag):
        self.what_var = tag

    def get_values(self):
        return (self.which_var, self.what_var)

    def disable(self):
        nb = self.lradio.index(Pmw.END)
        for i in range(nb+1):
            b = self.lradio.button(i)
            b.config(state='disabled')
            b.update()
        nb = self.rradio.index(Pmw.END)
        for i in range(nb+1):
            b = self.rradio.button(i)
            b.config(state='disabled')
            b.update()
        nb = self.buttons.index(Pmw.END)
        for i in range(nb+1):
            b = self.buttons.button(i)
            b.config(state='disabled')
            b.update()

    def enable(self):
        nb = self.lradio.index(Pmw.END)
        for i in range(nb+1):
            b = self.lradio.button(i)
            b.config(state='normal')
            b.update()
        nb = self.rradio.index(Pmw.END)
        for i in range(nb+1):
            b = self.rradio.button(i)
            b.config(state='normal')
            b.update()
        nb = self.action.index(Pmw.END)
        for i in range(nb+1):
            b = self.action.button(i)
            b.config(state='normal')
            b.update()
        nb = self.buttons.index(Pmw.END)
        for i in range(nb+1):
            b = self.buttons.button(i)
            b.config(state='normal')
            b.update()
        self.lradio.invoke('Current')


if __name__ == '__main__':
    root = Pmw.initialise()
    root.title("Short Fragment Toolbox")
    ActionChooser().pack()
    root.mainloop()
