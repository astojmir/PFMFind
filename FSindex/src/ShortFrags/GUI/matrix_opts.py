import Tkinter
import Pmw

class MatrixOpts(Tkinter.Frame):
    def __init__(self, parent=None):
        Tkinter.Frame.__init__(self,parent)
        self.pack()
        self.cases = []
        self.current = -1

        # Setup variables
        self.mat_var = Tkinter.StringVar()
        # Matrix Menu
        self.mat_menu = Pmw.OptionMenu(self,
                                       labelpos = 'w',
                                       label_text = 'Score\nMatrix:',
                                       menubutton_textvariable = self.mat_var,
                                       items = [],
                                       menubutton_width = 10,
                                       command = None)
        self.mat_menu.grid(row=0, column=0, padx = 3, pady = 5, sticky='w')

        # Conversion selection
        self.conv_select = Pmw.RadioSelect(self,
                                           buttontype = 'radiobutton',
                                           labelpos = 'n',
                                           label_text = 'Conversion:',
                                           hull_borderwidth = 2,
                                           hull_relief = 'ridge')

        # Scale
        self.wScale = Pmw.EntryField(self,
                                     labelpos = 'w',
                                     entry_width = 8,
                                     label_text = 'Scale:',
                                     validate = {'validator' : 'real',
                                                 'min' : 0.0},
                                     modifiedcommand = None)

        # Weights and regularisers
        self.vWeights = Tkinter.StringVar()
        self.vReg = Tkinter.StringVar()
        
        self.wWeights = Pmw.OptionMenu(self,
                                       labelpos = 'w',
                                       label_text = 'Weighting:',
                                       menubutton_textvariable = self.vWeights,
                                       items = [],
                                       menubutton_width = 10,
                                       command = None)
        self.wReg = Pmw.OptionMenu(self,
                                   labelpos = 'w',
                                   label_text = 'Regulariser:',
                                   menubutton_textvariable = self.vReg,
                                   items = [],
                                   menubutton_width = 20,
                                   command = None)
                
    def add_case(self, matrices, convs=[], mat_default=None,
                 conv_default=None, update_func=None,
                 scale_default=None, weights = [], weight_default=None,
                 regularisers = [], reg_default=None):

        if mat_default == None:
            mat_default = matrices[0]
            
        if conv_default == None and len(convs):
            conv_default = convs[0]

        if weight_default == None and len(weights):
            weight_default = weights[0]

        if reg_default == None and len(regularisers):
            reg_default = regularisers[0]

        C = {'matrix_list': matrices,
             'conv_list': convs,
             'matrix_default': mat_default,
             'conv_default': conv_default,
             'update_func': update_func,
             'scale_default': scale_default,
             'weights': weights,
             'weight_default': weight_default,
             'regularisers': regularisers,
             'reg_default': reg_default,
             }
        self.cases.append(C)
        return len(self.cases)-1

    def _change_cases(self):
        C = self.cases[self.current]

        self.mat_menu.setitems(C['matrix_list'])
        self.mat_menu.configure(command = C['update_func'])

        if C['scale_default'] != None:
            self.wScale.configure(modifiedcommand =
                                  lambda : C['update_func'](None))
            self.wScale.setvalue(str(C['scale_default']))
            self.wScale.grid(row=0, column=1, pady = 5, sticky='w')
        else:
            self.wScale.grid_forget()
            

        self.conv_select.deleteall()
        if len(C['conv_list']) > 0:
            self.conv_select.configure(command = C['update_func'])
            self.conv_select.grid(row=2, column=0, columnspan=2, padx = 3, pady = 1)
            for text in C['conv_list']:
                self.conv_select.add(text)
        else: 
            self.conv_select.grid_forget()

        if len(C['weights']) > 0:
            self.wWeights.setitems(C['weights'])
            self.wWeights.configure(command = C['update_func'])
            self.wWeights.grid(row=3, column=0, columnspan=2, padx = 3,
                               pady = 5, sticky='w')
        else: 
            self.wWeights.grid_forget()

        if len(C['regularisers']) > 0:
            self.wReg.setitems(C['regularisers'])
            self.wReg.configure(command = C['update_func'])
            self.wReg.grid(row=4, column=0, columnspan=2, padx = 3,
                           pady = 5, sticky='w')
        else: 
            self.wReg.grid_forget()

    def set_values(self, case=None, matrix=None, conv=None, scale=None,
                   weight=None, reg=None):
        if case != self.current and case != None:
            self.current = case
            self._change_cases()

        C = self.cases[self.current]

        if matrix == None:
            matrix = C['matrix_default']
        self.mat_var.set(matrix)

        if len(C['conv_list']):
            if conv == None:
                conv = C['conv_default']
            self.conv_select.setvalue(conv)

        if C['scale_default'] != None:
            if scale == None:
                scale = C['scale_default']
            self.wScale.setvalue(str(scale))

        if len(C['weights']):
            if weight == None:
                weight = C['weight_default']
            self.vWeights.set(weight)

        if len(C['regularisers']):
            if reg == None:
                reg = C['reg_default']
            self.vReg.set(reg)

        C['update_func'](None)

        
    def get_data(self):
        C = self.cases[self.current]
        res = {}
        res['matrix'] = self.mat_var.get()
        if len(C['conv_list']):
            res['conv_type'] = self.conv_select.getvalue()
        else:
            res['conv_type'] = -1
            
        if self.cases[self.current]['scale_default'] != None:
            val_str = self.wScale.getvalue()
            if val_str == '':
                res['scale'] = self.cases[self.current]['scale_default']
            else:
                res['scale'] = float(val_str)
        else:
            res['scale'] = None

        if len(C['weights']):
            res['weight'] = self.vWeights.get()
        else:
            res['weight'] = None

        if len(C['regularisers']):
            res['reg'] = self.vReg.get()
        else:
            res['reg'] = None
            
        return res
