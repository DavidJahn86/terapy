# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:35:05 2015

@author: jahndav
"""
from PyQt4 import QtGui
import sys
from os import path
from ui_formatdialog import Ui_Dialog

class FormatDialog(QtGui.QDialog):
    def __init__(self):
        super(FormatDialog, self).__init__()        

        self.ui=Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.cb_dataFormat.currentIndexChanged.connect(self.changeFormat)
   
    def setFilenames(self,fn):
        self.filenames=fn
        for fn in self.filenames:
            self.ui.lw_inputfiles.addItem(QtGui.QListWidgetItem(path.split(fn)[1]))
        
#        self.autoDetectFileFormat()
    
    def getDataFormat(self):
        dataformat={'time_factor':float(self.ui.le_timescaling.text()),
                'time_col':self.ui.sb_timecol.value(),
                'X_col':self.ui.sb_dataxcol.value(),
                'Y_col':self.ui.sb_dataycol.value(),
                'dec_sep':self.ui.cb_separator.currentText(),
                'skiprows':self.ui.sb_skiprows.value()}    

        return dataformat
    
    def doAveraging(self):
        return self.ui.cb_averaging.isChecked()

    def autoDetectFileFormat(self):
        pass
        
        
    def changeFormat(self):
        "0: Marburg"
        "1: INRIM"
        "2: custom"
        
        if self.ui.cb_dataFormat.currentIndex()==0:
            self.ui.cb_separator.setCurrentIndex(0)
            self.ui.le_timescaling.setText("1")
            self.ui.sb_dataxcol.setValue(1)
            self.ui.sb_dataycol.setValue(2)
            self.ui.sb_skiprows.setValue(0)
            self.ui.sb_timecol.setValue(0)
            
        if self.ui.cb_dataFormat.currentIndex()==1:
            self.ui.cb_separator.setCurrentIndex(1)
            self.ui.le_timescaling.setText("1")
            self.ui.sb_dataxcol.setValue(3)
            self.ui.sb_dataycol.setValue(5)
            self.ui.sb_skiprows.setValue(0)
            self.ui.sb_timecol.setValue(2)
            
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = FormatDialog()
    window.show()
    sys.exit(app.exec_())
