import Tkinter
import Pmw

class TextDisplay(Pmw.ScrolledText):
    def __init__(self, parent=None, label=None, text=None):
        fixedFont = Pmw.logicalfont('Fixed')
        Pmw.ScrolledText.__init__(self, parent,
                # borderframe = 1,
                labelpos = 'n',
                label_text=label,
                usehullsize = 1,
                text_wrap='none',
                text_font = fixedFont,
                text_padx = 4,
                text_pady = 4)

        if text:
            self.insert('end', text)
        self.configure(text_state = 'disabled')

    def set_size(self, height=400, width=300):
        self.configure(hull_height=height, hull_width=width)
        
    def set_text(self, label=None, text=None):
        if label:
            self.configure(label_text=label)
        if text:
            self.configure(text_state = 'normal')
            self.setvalue(text)
            self.configure(text_state = 'disabled')
        
    def insert_text(self, text):
        self.configure(text_state = 'normal')
        self.insert('end', text)
        self.configure(text_state = 'disabled')
