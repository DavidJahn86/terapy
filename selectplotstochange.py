# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 08:29:45 2015

@author: david
"""
from PyQt4 import QtGui
import sys

from ui_selectplotstochange import Ui_Dialog

class PlotsToChange(QtGui.QDialog):
    def __init__(self,availablePlots,label_description):
        super(PlotsToChange, self).__init__()        
        
        self.ui=Ui_Dialog()
        self.ui.setupUi(self)
        for tp in availablePlots:
            self.ui.cb_whichtochange.addItem(tp)
        
        self.ui.lb_description.setText(label_description)
        
    def getPlotsToChange(self):
        return self.ui.cb_whichtochange.currentIndex()
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = PlotsToChange(['Plot 1'],'Select curve for which the Dynamic range will be plotted')
    window.show()
    sys.exit(app.exec_())
