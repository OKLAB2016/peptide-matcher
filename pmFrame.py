from gui import MyFrame
from peptideMatcher import PeptideMatcher
import wx, xlsxwriter

class PMFrame(MyFrame):

    def on_button_fasta(self, event):
        dlg = wx.FileDialog(self, message='Choose database fasta file', defaultFile='', style=(wx.FD_OPEN | wx.FD_CHANGE_DIR))
        if dlg.ShowModal() == wx.ID_OK:
            fastaPath = dlg.GetPath()
            self.text_fasta.SetValue(fastaPath)
            peptidesPath = self.text_peptides.GetValue()
            self.button_run.Enable(peptidesPath != '')
        dlg.Destroy()

    def on_button_peptides(self, event):
        dlg = wx.FileDialog(self, message='Choose file with peptides', defaultFile='', style=(wx.FD_OPEN | wx.FD_CHANGE_DIR))
        if dlg.ShowModal() == wx.ID_OK:
            peptidesPath = dlg.GetPath()
            self.text_peptides.SetValue(peptidesPath)
            fastaPath = self.text_fasta.GetValue()
            self.button_run.Enable(fastaPath != '')
        dlg.Destroy()

    def on_button_run(self, event):
        self.grid_matches.ClearGrid()
        num_rows = self.grid_matches.GetNumberRows()
        if num_rows > 0:
            self.grid_matches.DeleteRows(numRows=num_rows)
        fasta = self.text_fasta.GetValue()
        secstruct_included = self.radio_box_secstruct.GetSelection()
        peptides = self.text_peptides.GetValue()
        flanks = int(self.spin_flanks.GetValue())
        progress_dialog = wx.ProgressDialog('Matching...', 'Matching the peptides...', parent=self, style=(wx.PD_APP_MODAL | wx.PD_AUTO_HIDE))
        progress_dialog.Pulse()
        peptide_matcher = PeptideMatcher(fasta, secstruct_included, peptides, flanks, self.grid_matches, progress_dialog)

        try:
            peptide_matcher.run()
            self.button_save.Enable(True)
        except Exception as e:
            try:
                msg_dialog = wx.MessageDialog(self, str(e), 'Error', wx.OK | wx.ICON_ERROR)
                msg_dialog.ShowModal()
                msg_dialog.Destroy()
            finally:
                e = None
                del e

        progress_dialog.Destroy()

    def on_button_save(self, event):
        dlg = wx.FileDialog(self, message='Choose file to save the output', wildcard='MS Excel 2007 Spreadsheets (*.xlsx)|*.xlsx', style=(wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT))
        if dlg.ShowModal() == wx.ID_OK:
            outputPath = dlg.GetPath()
            table = self.grid_matches.GetTable()
            row_num = table.GetRowsCount()
            col_num = table.GetColsCount()
            workbook = xlsxwriter.Workbook(outputPath)
            worksheet = workbook.add_worksheet()
            cell_format = workbook.add_format({'bold': True})
            worksheet.set_row(0, None, cell_format)
            for col in range(col_num):
                val = self.grid_matches.GetColLabelValue(col)
                worksheet.write(0, col, val)

            for row in range(row_num):
                for col in range(col_num):
                    val = table.GetValue(row, col)
                    worksheet.write(row + 1, col, val)

            workbook.close()
        dlg.Destroy()
