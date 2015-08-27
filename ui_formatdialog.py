# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'formatdialog.ui'
#
# Created: Thu Aug 13 12:53:18 2015
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
        Dialog.resize(481, 311)
        self.verticalLayout_2 = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.splitter = QtGui.QSplitter(Dialog)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.groupBox_2 = QtGui.QGroupBox(self.splitter)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.lw_inputfiles = QtGui.QListWidget(self.groupBox_2)
        self.lw_inputfiles.setObjectName(_fromUtf8("lw_inputfiles"))
        self.verticalLayout.addWidget(self.lw_inputfiles)
        self.label_7 = QtGui.QLabel(self.groupBox_2)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.verticalLayout.addWidget(self.label_7)
        self.cb_dataFormat = QtGui.QComboBox(self.groupBox_2)
        self.cb_dataFormat.setObjectName(_fromUtf8("cb_dataFormat"))
        self.cb_dataFormat.addItem(_fromUtf8(""))
        self.cb_dataFormat.addItem(_fromUtf8(""))
        self.cb_dataFormat.addItem(_fromUtf8(""))
        self.verticalLayout.addWidget(self.cb_dataFormat)
        self.cb_averaging = QtGui.QCheckBox(self.groupBox_2)
        self.cb_averaging.setChecked(True)
        self.cb_averaging.setObjectName(_fromUtf8("cb_averaging"))
        self.verticalLayout.addWidget(self.cb_averaging)
        self.groupBox = QtGui.QGroupBox(self.splitter)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.formLayout = QtGui.QFormLayout(self.groupBox)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.FieldsStayAtSizeHint)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.cb_separator = QtGui.QComboBox(self.groupBox)
        self.cb_separator.setMaximumSize(QtCore.QSize(45, 16777215))
        self.cb_separator.setObjectName(_fromUtf8("cb_separator"))
        self.cb_separator.addItem(_fromUtf8(""))
        self.cb_separator.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.cb_separator)
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_6)
        self.sb_timecol = QtGui.QSpinBox(self.groupBox)
        self.sb_timecol.setObjectName(_fromUtf8("sb_timecol"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.sb_timecol)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_3)
        self.sb_dataxcol = QtGui.QSpinBox(self.groupBox)
        self.sb_dataxcol.setProperty("value", 1)
        self.sb_dataxcol.setObjectName(_fromUtf8("sb_dataxcol"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.sb_dataxcol)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.label_5)
        self.sb_dataycol = QtGui.QSpinBox(self.groupBox)
        self.sb_dataycol.setProperty("value", 2)
        self.sb_dataycol.setObjectName(_fromUtf8("sb_dataycol"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.sb_dataycol)
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_4)
        self.sb_skiprows = QtGui.QSpinBox(self.groupBox)
        self.sb_skiprows.setObjectName(_fromUtf8("sb_skiprows"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.sb_skiprows)
        self.le_timescaling = QtGui.QLineEdit(self.groupBox)
        self.le_timescaling.setMaximumSize(QtCore.QSize(45, 16777215))
        self.le_timescaling.setObjectName(_fromUtf8("le_timescaling"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.le_timescaling)
        self.verticalLayout_2.addWidget(self.splitter)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        self.cb_separator.setCurrentIndex(0)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.groupBox_2.setTitle(_translate("Dialog", "Files to be imported:", None))
        self.label_7.setText(_translate("Dialog", "Use the following format:", None))
        self.cb_dataFormat.setItemText(0, _translate("Dialog", "Marburg data format", None))
        self.cb_dataFormat.setItemText(1, _translate("Dialog", "Inrim data format", None))
        self.cb_dataFormat.setItemText(2, _translate("Dialog", "custom data format", None))
        self.cb_averaging.setText(_translate("Dialog", "Average Data", None))
        self.groupBox.setTitle(_translate("Dialog", "Preferences:", None))
        self.label.setText(_translate("Dialog", "Rescale Time Axis by:", None))
        self.label_2.setText(_translate("Dialog", "Decimal Separator:", None))
        self.cb_separator.setItemText(0, _translate("Dialog", ",", None))
        self.cb_separator.setItemText(1, _translate("Dialog", ".", None))
        self.label_6.setText(_translate("Dialog", "Time Column:", None))
        self.label_3.setText(_translate("Dialog", "Data Column (X):", None))
        self.label_5.setText(_translate("Dialog", "Data Column (Y):", None))
        self.label_4.setText(_translate("Dialog", "Skip rows: ", None))
        self.le_timescaling.setText(_translate("Dialog", "1", None))

