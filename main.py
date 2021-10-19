from pmFrame import PMFrame
import wx

class PeptideMatcher(wx.App):

    def OnInit(self):
        self.frame = PMFrame(None, wx.ID_ANY, '')
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True


if __name__ == '__main__':
    peptide_matcher = PeptideMatcher(0)
    peptide_matcher.MainLoop()
