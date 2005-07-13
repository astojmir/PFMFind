class View:
    """
    Meant to be used as mixin class with something
    derived from Tkinter.Frame (or set_size needs
    to be overwritten).
    """
    def __init__(self):
        pass
    
    def set_size(self, height=400, width=300):
        self.configure(height=height, width=width)
        
    def reset(self, state=(None, None, None)):
        pass
