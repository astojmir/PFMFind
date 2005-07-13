from ShortFrags.GUI.view import View
from ShortFrags.GUI.text_display import TextDisplay

class IndexView(TextDisplay, View):
    def __init__(self, parent, PFMF_client):
        TextDisplay.__init__(self, parent, label='Index')
        View.__init__(self)
        self.PFMF_client = PFMF_client
        
    def reset(self, state):
        text = self.PFMF_client.index_data_str()
        self.set_text('Index', text)
