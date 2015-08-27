# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'selectplotstochange.ui'
#
# Created: Sat Aug 15 08:46:34 2015
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(400, 140)
        self.verticalLayout = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.lb_description = QtGui.QLabel(Dialog)
        self.lb_description.setObjectName(_fromUtf8("lb_description"))
        self.verticalLayout.addWidget(self.lb_description)
        self.cb_whichtochange = QtGui.QComboBox(Dialog)
        self.cb_whichtochange.setObjectName(_fromUtf8("cb_whichtochange"))
        self.cb_whichtochange.addItem(_fromUtf8(""))
        self.cb_whichtochange.addItem(_fromUtf8(""))
        self.verticalLayout.addWidget(self.cb_whichtochange)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.lb_description.setText(_translate("Dialog", "spacer", None))
        self.cb_whichtochange.setItemText(0, _translate("Dialog", "All Selected", None))
        self.cb_whichtochange.setItemText(1, _translate("Dialog", "All", None))

