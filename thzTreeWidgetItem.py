# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 13:01:05 2015

@author: jahndav
"""

from PyQt4 import QtGui

class THzTreeWidgetItem(QtGui.QTreeWidgetItem):
    def __init__(self):
        super(THzTreeWidgetItem,self).__init__()
        self.tdData=None
        self.fdData=None
        self.tdline=[]
        self.fdline=[]
        self.color=None
        
    def removeLines(self):
        for line in self.tdline:
            line.remove()

        for line in self.fdline:
            line.remove()
        
        print self.fdline
        
    def __del__(self):
        print "item destructed"
        
        