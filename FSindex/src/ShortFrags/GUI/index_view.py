from ShortFrags.GUI.view import View
from ShortFrags.GUI.text_display import TextDisplay

class IndexView(TextDisplay, View):
    def __init__(self, parent=None):
        TextDisplay.__init__(self, parent, label='Index')
        View.__init__(self)

    def reset(self, state):
        text = state['FE'].search_client.index_data_str()
        self.set_text('Index', text)
