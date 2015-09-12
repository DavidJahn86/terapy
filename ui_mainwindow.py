# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindow.ui'
#
# Created: Wed Sep  9 15:29:27 2015
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

class Ui_TeraView(object):
    def setupUi(self, TeraView):
        TeraView.setObjectName(_fromUtf8("TeraView"))
        TeraView.resize(1247, 475)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(TeraView.sizePolicy().hasHeightForWidth())
        TeraView.setSizePolicy(sizePolicy)
        self.centralWidget = QtGui.QWidget(TeraView)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralWidget.sizePolicy().hasHeightForWidth())
        self.centralWidget.setSizePolicy(sizePolicy)
        self.centralWidget.setObjectName(_fromUtf8("centralWidget"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.centralWidget)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralWidget)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.horizontalLayout_5 = QtGui.QHBoxLayout(self.tab)
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.splitter_2 = QtGui.QSplitter(self.tab)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName(_fromUtf8("splitter_2"))
        self.layoutWidget = QtGui.QWidget(self.splitter_2)
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.spectrumCanvas = MatplotlibWidget(self.layoutWidget)
        self.spectrumCanvas.setObjectName(_fromUtf8("spectrumCanvas"))
        self.verticalLayout.addWidget(self.spectrumCanvas)
        self.splitter = QtGui.QSplitter(self.splitter_2)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.fileTree = QtGui.QTreeWidget(self.splitter)
        self.fileTree.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.fileTree.setObjectName(_fromUtf8("fileTree"))
        self.horizontalLayout_5.addWidget(self.splitter_2)
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.tabWidget.addTab(self.tab_3, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.tabWidget)
        self.horizontalLayout_3.addLayout(self.horizontalLayout)
        TeraView.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(TeraView)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1247, 25))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuPlots = QtGui.QMenu(self.menuBar)
        self.menuPlots.setObjectName(_fromUtf8("menuPlots"))
        self.menuSpectrum_Analysis = QtGui.QMenu(self.menuPlots)
        self.menuSpectrum_Analysis.setObjectName(_fromUtf8("menuSpectrum_Analysis"))
        self.menuData_Operations = QtGui.QMenu(self.menuBar)
        self.menuData_Operations.setObjectName(_fromUtf8("menuData_Operations"))
        self.menuTime_Domain = QtGui.QMenu(self.menuData_Operations)
        self.menuTime_Domain.setObjectName(_fromUtf8("menuTime_Domain"))
        TeraView.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(TeraView)
        self.mainToolBar.setObjectName(_fromUtf8("mainToolBar"))
        TeraView.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.mainStatus = QtGui.QStatusBar(TeraView)
        self.mainStatus.setSizeGripEnabled(True)
        self.mainStatus.setObjectName(_fromUtf8("mainStatus"))
        TeraView.setStatusBar(self.mainStatus)
        self.actionLoad = QtGui.QAction(TeraView)
        self.actionLoad.setObjectName(_fromUtf8("actionLoad"))
        self.actionClose = QtGui.QAction(TeraView)
        self.actionClose.setObjectName(_fromUtf8("actionClose"))
        self.actionPlot_uncertainty_intervals = QtGui.QAction(TeraView)
        self.actionPlot_uncertainty_intervals.setObjectName(_fromUtf8("actionPlot_uncertainty_intervals"))
        self.actionPlot_SNR = QtGui.QAction(TeraView)
        self.actionPlot_SNR.setObjectName(_fromUtf8("actionPlot_SNR"))
        self.actionPlot_Dynamic_Range = QtGui.QAction(TeraView)
        self.actionPlot_Dynamic_Range.setObjectName(_fromUtf8("actionPlot_Dynamic_Range"))
        self.actionNoise_Floor_calculated = QtGui.QAction(TeraView)
        self.actionNoise_Floor_calculated.setObjectName(_fromUtf8("actionNoise_Floor_calculated"))
        self.actionNoise_Floor_from_fft = QtGui.QAction(TeraView)
        self.actionNoise_Floor_from_fft.setObjectName(_fromUtf8("actionNoise_Floor_from_fft"))
        self.actionFrequency_Domain = QtGui.QAction(TeraView)
        self.actionFrequency_Domain.setObjectName(_fromUtf8("actionFrequency_Domain"))
        self.actionZero_Padding = QtGui.QAction(TeraView)
        self.actionZero_Padding.setObjectName(_fromUtf8("actionZero_Padding"))
        self.actionWindowing = QtGui.QAction(TeraView)
        self.actionWindowing.setObjectName(_fromUtf8("actionWindowing"))
        self.actionScale_Data = QtGui.QAction(TeraView)
        self.actionScale_Data.setObjectName(_fromUtf8("actionScale_Data"))
        self.menuFile.addAction(self.actionLoad)
        self.menuFile.addAction(self.actionClose)
        self.menuSpectrum_Analysis.addAction(self.actionPlot_uncertainty_intervals)
        self.menuSpectrum_Analysis.addAction(self.actionPlot_SNR)
        self.menuSpectrum_Analysis.addAction(self.actionPlot_Dynamic_Range)
        self.menuSpectrum_Analysis.addAction(self.actionNoise_Floor_calculated)
        self.menuSpectrum_Analysis.addAction(self.actionNoise_Floor_from_fft)
        self.menuPlots.addAction(self.menuSpectrum_Analysis.menuAction())
        self.menuTime_Domain.addAction(self.actionZero_Padding)
        self.menuTime_Domain.addAction(self.actionWindowing)
        self.menuTime_Domain.addAction(self.actionScale_Data)
        self.menuData_Operations.addAction(self.menuTime_Domain.menuAction())
        self.menuData_Operations.addAction(self.actionFrequency_Domain)
        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuPlots.menuAction())
        self.menuBar.addAction(self.menuData_Operations.menuAction())

        self.retranslateUi(TeraView)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(TeraView)

    def retranslateUi(self, TeraView):
        TeraView.setWindowTitle(_translate("TeraView", "MainWindow", None))
        self.fileTree.headerItem().setText(0, _translate("TeraView", "Plot", None))
        self.fileTree.headerItem().setText(1, _translate("TeraView", "Name", None))
        self.fileTree.headerItem().setText(2, _translate("TeraView", "Time Domain Info", None))
        self.fileTree.headerItem().setText(3, _translate("TeraView", "Frequency Domain Info", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("TeraView", "Spectrum Analysis", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("TeraView", "Image", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("TeraView", "Optical Constants", None))
        self.menuFile.setTitle(_translate("TeraView", "File", None))
        self.menuPlots.setTitle(_translate("TeraView", "Plots", None))
        self.menuSpectrum_Analysis.setTitle(_translate("TeraView", "Spectrum Analysis", None))
        self.menuData_Operations.setTitle(_translate("TeraView", "Data Operations", None))
        self.menuTime_Domain.setTitle(_translate("TeraView", "Time Domain", None))
        self.actionLoad.setText(_translate("TeraView", "Import", None))
        self.actionClose.setText(_translate("TeraView", "Close", None))
        self.actionPlot_uncertainty_intervals.setText(_translate("TeraView", "plot uncertainty intervals", None))
        self.actionPlot_SNR.setText(_translate("TeraView", "plot SNR", None))
        self.actionPlot_Dynamic_Range.setText(_translate("TeraView", "plot Dynamic Range", None))
        self.actionNoise_Floor_calculated.setText(_translate("TeraView", "Noise Floor calculated", None))
        self.actionNoise_Floor_from_fft.setText(_translate("TeraView", "Noise Floor from fft", None))
        self.actionFrequency_Domain.setText(_translate("TeraView", "Frequency Domain", None))
        self.actionZero_Padding.setText(_translate("TeraView", "Zero Padding", None))
        self.actionWindowing.setText(_translate("TeraView", "Windowing", None))
        self.actionScale_Data.setText(_translate("TeraView", "Modify Data", None))

from matplotlibwidget import MatplotlibWidget
