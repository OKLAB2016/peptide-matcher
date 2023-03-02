from peptide_matcher.gui import BasicFrame
from peptide_matcher import PeptideMatcher, wrap_logos, wrap_scores
import wx
import xlsxwriter
import traceback
from time import time

class PMFrame(BasicFrame):

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

        time_start = time()
        with open(peptides) as peptides_fp:
            peptide_matcher = PeptideMatcher(peptides_fp, fasta, secstruct_included, flanks)
            try:
                row = 0
                for output in peptide_matcher.run():
                    peptide = output['peptide']
                    peplen = len(peptide)
                    if output['matches']:
                        n_logo = wrap_logos(output['n_logos'])
                        c_logo = wrap_logos(output['c_logos'])
                        for match in output['matches']:
                            self.grid_matches.AppendRows()
                            self.grid_matches.SetCellValue(row, 0, peptide)
                            self.grid_matches.SetCellValue(row, 1, str(peplen))
                            self.grid_matches.SetCellValue(row, 2, match['record_id'])
                            self.grid_matches.SetCellValue(row, 3, str(match['start']))
                            self.grid_matches.SetCellValue(row, 4, str(match['end']))
                            self.grid_matches.SetCellValue(row, 5, str(match['c_term']))
                            self.grid_matches.SetCellValue(row, 6, match['n_flank'])
                            self.grid_matches.SetCellValue(row, 7, match['c_flank'])
                            self.grid_matches.SetCellValue(row, 8, n_logo)
                            self.grid_matches.SetCellValue(row, 9, c_logo)
                            if secstruct_included:
                                self.grid_matches.SetCellValue(row, 10, match['sst_n_term'])
                                self.grid_matches.SetCellValue(row, 11, match['sst_pept'])
                                self.grid_matches.SetCellValue(row, 12, match['sst_c_term'])
                                self.grid_matches.SetCellValue(row, 13, match['tm_n_term'])
                                self.grid_matches.SetCellValue(row, 14, match['tm_pept'])
                                self.grid_matches.SetCellValue(row, 15, match['tm_c_term'])
                                self.grid_matches.SetCellValue(row, 16, wrap_scores(match['conf_n_term']))
                                self.grid_matches.SetCellValue(row, 17, wrap_scores(match['conf_pept']))
                                self.grid_matches.SetCellValue(row, 18, wrap_scores(match['conf_c_term']))
                                self.grid_matches.SetCellValue(row, 19, wrap_scores(match['acc_n_term']))
                                self.grid_matches.SetCellValue(row, 20, wrap_scores(match['acc_pept']))
                                self.grid_matches.SetCellValue(row, 21, wrap_scores(match['acc_c_term']))
                        row += 1
                    else:
                        self.grid_matches.AppendRows()
                        self.grid_matches.SetCellValue(row, 0, peptide)
                        self.grid_matches.SetCellValue(row, 1, str(peplen))
                        self.grid_matches.SetCellValue(row, 2, 'No match')
                        row += 1
                    progress_dialog.Pulse()

                self.button_save.Enable(True)
            except Exception as e:
                try:
                    msg_dialog = wx.MessageDialog(self, '%s:\n%s' % (str(e), traceback.format_exc()), 'Error', wx.OK | wx.ICON_ERROR)
                    msg_dialog.ShowModal()
                    msg_dialog.Destroy()
                finally:
                    e = None
                    del e

        print("Run time: " + str(time() - time_start))
        progress_dialog.Destroy()

    def on_button_save(self, event):
        dlg = wx.FileDialog(self, message='Choose file to save the output', wildcard='MS Excel 2007 Spreadsheets (*.xlsx)|*.xlsx', style=(wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT))
        if dlg.ShowModal() == wx.ID_OK:
            outputPath = dlg.GetPath()
            table = self.grid_matches.GetTable()
            row_num = table.GetNumberRows()
            col_num = table.GetNumberCols()
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
