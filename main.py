import sys
#for dynamic loading of the GUI

from PyQt4 import QtGui, uic, QtCore
import numpy as np
import matplotlib.pyplot as plt
import TeraData
from uncertainties import unumpy
from os import path
import matplotlib.gridspec as gridspec
import glob

from thzTreeWidgetItem import THzTreeWidgetItem
from formatdialog import FormatDialog
from selectplotstochange import PlotsToChange
from ui_mainwindow import Ui_TeraView

#for matplotlib widget


class MyWindow(QtGui.QMainWindow):
    def __init__(self):
        super(MyWindow, self).__init__()        
        
        self.ui=Ui_TeraView()
        self.ui.setupUi(self)
        self.initializeSpectrumCanvas()
            
        self.ui.actionClose.triggered.connect(QtGui.qApp.quit)
        self.ui.actionLoad.triggered.connect(self.loadFile)

        self.ui.actionPlot_Dynamic_Range.triggered.connect(self.plotDR)
        self.ui.actionPlot_SNR.triggered.connect(self.plotSNR)
        self.ui.actionPlot_uncertainty_intervals.triggered.connect(self.plotuncertainty)

        #Actions in Data Operations Menu
        self.ui.actionScale_Data.triggered.connect(self.showDataManipulation)
        self.ui.actionWindowing.triggered.connect(self.showWindowing) 
        self.ui.actionZero_Padding.triggered.connect(self.showZeroPadding)
        self.ui.actionInterpolate_Data.triggered.connect(self.showTDInterpolation)
        self.ui.actionAverage_Datasets.triggered.connect(self.showAveraging)
        
        ##connect all the buttons in preferences menus
        self.ui.preferences.hide()

        #this is the first combobox of this type, but i will use it at several places again
        graphsModel=self.ui.cbwhichGraphs.model()
        self.ui.cbwhichGraphs_windowing.setModel(graphsModel)
        self.ui.cbwhichGraphs_zeropadding.setModel(graphsModel)
        self.ui.cbwhichGraphs_interpolation.setModel(graphsModel)
        self.ui.cbwhichGraphs_5.setModel(graphsModel)
        self.ui.cbwhichGraphs_6.setModel(graphsModel)
        self.ui.cbwhichGraphs_7.setModel(graphsModel)
        
        #Buttons  in Spectrum analysis
        self.ui.pbApplyTDManipulation.clicked.connect(self.applyTDManipulation)        
        self.ui.pbPreviewTDManipulation.clicked.connect(self.dataManipulationTemporarilyUpdatePlot) 
        self.ui.pbCancelTDManipulation.clicked.connect(self.cancelPreferences)
        
        self.ui.pbApplyWindowing.clicked.connect(self.applyWindowing)
        self.ui.pbPreviewWindowing.clicked.connect(self.previewWindowing)
        self.ui.pbCancelWindowing.clicked.connect(self.cancelPreferences)
        
        self.ui.pbApplyPadding.clicked.connect(self.applyZeroPadding)
        self.ui.pbPreviewPadding.clicked.connect(self.previewZeroPadding)
        self.ui.pbCancelPadding.clicked.connect(self.cancelPreferences)        
        
        self.ui.pbApplyInterpolation.clicked.connect(self.applyTDInterpolation)
        self.ui.pbPreviewInterpolation.clicked.connect(self.previewTDInterpolation)
        self.ui.pbCancelInterpolation.clicked.connect(self.cancelPreferences)

        self.ui.pbApplyAveraging.clicked.connect(self.applyAveraging)
        self.ui.pbPreviewAveraging.clicked.connect(self.previewAveraging)
        self.ui.pbCancelAveraging.clicked.connect(self.cancelPreferences)
        
        #TreeView Behaviour        
        self.ui.fileTree.itemChanged.connect(self.updateSpectrumAnalysisPlot)
        self.ui.fileTree.itemDoubleClicked.connect(self.onTreeWidgetItemDoubleClicked)
        self.ui.fileTree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.fileTree.customContextMenuRequested.connect(self.onCustomContextMenu)
                
        self.removeItemAction=QtGui.QAction(self)
        self.removeItemAction.triggered.connect(self.removeCurrentlySelectedItem)
        self.removeItemAction.setText('Remove Selected Entry')
        
        self.TopLevelItemMenu=QtGui.QMenu(self)
        self.ChildMenu=QtGui.QMenu(self)
                
        self.TopLevelItemMenu.addAction(self.ui.menuSpectrum_Analysis.menuAction())
        self.TopLevelItemMenu.addSeparator()
        self.TopLevelItemMenu.addAction(self.removeItemAction)        
        self.ChildMenu.addAction(self.removeItemAction)
         
        #for testing:
        params={'time_factor':1,
            'time_col':0,
            'X_col':1,
            'Y_col':2,
            'dec_sep':'.',
            'skiprows':0}
        files=glob.glob('*Reference*.txt')
        tdData=TeraData.TimeDomainData.importMultipleFiles(files,params)
        tdData.setDataSetName('test')
        x=self.fillTree('test',tdData)          
         
        self.show()
        
    def showAveraging(self):
        self.ui.preferences.setCurrentIndex(4)
        self.ui.lwAverage.clear()
        
        for ix in range(self.ui.fileTree.topLevelItemCount()):
            it=self.ui.fileTree.topLevelItem(ix)
            x=QtGui.QListWidgetItem()
            x.setCheckState(QtCore.Qt.Checked)
            x.setText(it.text(1))
            col=QtGui.QColor(int(it.color[0]*255),int(it.color[1]*255),int(it.color[2]*255),int(it.color[3]*255))
            x.setBackgroundColor(col)
            self.ui.lwAverage.addItem(x)

        self.ui.leNewName_averaging.setText('Average')
        self.ui.preferences.show()

    def previewAveraging(self):
        
        av_where=self.ui.cbAverageWhere.currentIndex() #0 for time domain 1 for fd averaging
        
        tdata=[]        
        fdata=[]
        for ix in range(self.ui.lwAverage.count()):
            if self.ui.lwAverage.item(ix).checkState():
                tdata.append(self.ui.fileTree.topLevelItem(ix).tdData)
                fdata.append(self.ui.fileTree.topLevelItem(ix).fdData)
                  
        if av_where==0 and len(tdata)>1:
            av_tdata=TeraData.TimeDomainData.averageTimeDomainDatas(tdata)
            lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)
       
            if not lastitem.tempitem:
                tempItem=self.fillTree('__tempItem__',av_tdata)
                tempItem.tempitem=True
                tempItem.color=plt.cm.colors.colorConverter.to_rgba('0.75')
                col=QtGui.QColor(int(tempItem.color[0]*255),int(tempItem.color[1]*255),int(tempItem.color[2]*255),int(tempItem.color[3]*255))
                self.refreshCanvas()
            else:
                self.temporaryUpdate(lastitem,av_tdata)            
        
        #averaging in frequency domain not yet implemented here
        self.refreshCanvas()
   
    def applyAveraging(self):
        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.previewAveraging()
        
        self.insertTemporaryCopy()
        self.ui.preferences.hide()
        
    def applyTDInterpolation(self):
        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.previewTDInterpolation()
        
        where=self.ui.cbwhichGraphs_interpolation.currentIndex()
        
        if self.ui.rbManipulateCopy_interpolation.isChecked():
            #should be the temporary data 
            self.insertTemporaryCopy()
            self.ui.preferences.hide()
        else:
            self.insertTemporaryCopy(self.ui.fileTree.topLevelItem(where))
            self.ui.preferences.hide()                    
     

    
    def showTDInterpolation(self):
        self.ui.preferences.setCurrentIndex(3)
        self.ui.le_newName_interpolation.setText("Copy of " + self.ui.cbwhichGraphs_interpolation.currentText())
        self.ui.preferences.show()

    def previewTDInterpolation(self):
        
        step=self.ui.dsbinterpolationStep.value()*1e-15        
        
        where=self.ui.cbwhichGraphs_windowing.currentIndex()        
        tempItem=self.createTemporaryCopy(self.ui.fileTree.topLevelItem(where))
        
        #manipulate the temporary data
        tdata=tempItem.tdData
        mint=min(tdata.getTimeAxisRef())
        maxt=max(tdata.getTimeAxisRef())
        desiredLength=int((maxt-mint)/step)
        tdata=tdata.getInterpolatedTimeDomainData(desiredLength,mint,maxt)
        self.temporaryUpdate(tempItem,tdata)
        self.refreshCanvas()
    
    def applyTDInterpolation(self):
        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.previewTDInterpolation()
        
        where=self.ui.cbwhichGraphs_interpolation.currentIndex()
        
        if self.ui.rbManipulateCopy_interpolation.isChecked():
            #should be the temporary data 
            self.insertTemporaryCopy()
            self.ui.preferences.hide()
        else:
            self.insertTemporaryCopy(self.ui.fileTree.topLevelItem(where))
            self.ui.preferences.hide()            
        
        
    def showZeroPadding(self):
        self.ui.preferences.setCurrentIndex(2)
        self.ui.le_newName_padding.setText("Copy of " + self.ui.cbwhichGraphs_zeropadding.currentText())
        self.ui.preferences.show()

    def previewZeroPadding(self):
        
        mode=self.ui.cbZeroPadding.currentIndex()
        
        where=self.ui.cbwhichGraphs_windowing.currentIndex()        
        no_zeros=self.ui.sbNumberZeros.value()
        
        tempItem=self.createTemporaryCopy(self.ui.fileTree.topLevelItem(where))
        
        #manipulate the temporary data
        tdata=tempItem.tdData
        if mode==0:
            tdata=tdata.zeroPaddData(no_zeros,paddmode='zero',where='end')   
        elif mode==1:
            tdata=tdata.zeroPaddData(no_zeros,paddmode='zero',where='start')   
        elif mode==2:
            tdata=tdata.zeroPaddData(no_zeros,paddmode='gaussian',where='end')    
        elif mode ==3:
            tdata=tdata.zeroPaddData(no_zeros,paddmode='gaussian',where='start')   

        self.temporaryUpdate(tempItem,tdata)
        self.refreshCanvas()
    
    def applyZeroPadding(self):
        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.previewZeroPadding()
        
        where=self.ui.cbwhichGraphs_zeropadding.currentIndex()
        
        if self.ui.rbManipulateCopy_Padding.isChecked():
            #should be the temporary data 
            self.insertTemporaryCopy()
            self.ui.preferences.hide()
        else:
            self.insertTemporaryCopy(self.ui.fileTree.topLevelItem(where))
            self.ui.preferences.hide()    
    
    def showWindowing(self):
        self.ui.preferences.setCurrentIndex(1)
        self.ui.le_newName_windowing.setText("Copy of " + self.ui.cbwhichGraphs_windowing.currentText())
        self.ui.preferences.show()

    def previewWindowing(self):
        windowfunction=self.ui.cbwindowfunction.currentText()
        risingEdgelen=self.ui.dsb_risingEdge.value()*1e-12
        where=self.ui.cbwhichGraphs_windowing.currentIndex()        
        
        tempItem=self.createTemporaryCopy(self.ui.fileTree.topLevelItem(where))
        
        #manipulate the temporary data
        tdata=tempItem.tdData
        tdata=tdata.getWindowedData(risingEdgelen,windowfunction)
        self.temporaryUpdate(tempItem,tdata)
        self.refreshCanvas()
    
    def applyWindowing(self):
        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.previewWindowing()
        
        where=self.ui.cbwhichGraphs_windowing.currentIndex()
        
        if self.ui.rbManipulateCopy_windowing.isChecked():
            #should be the temporary data 
            self.insertTemporaryCopy()
            self.ui.preferences.hide()
        else:
            self.insertTemporaryCopy(self.ui.fileTree.topLevelItem(where))
            self.ui.preferences.hide()
            
    def showDataManipulation(self):
        #check if this is the first manipulation of the original data or if a tempdata entry is already existing
        self.ui.preferences.setCurrentIndex(0)
        self.ui.leNewName.setText("Copy of " + self.ui.cbwhichGraphs.currentText())
        self.ui.preferences.show()
        
    def dataManipulationTemporarilyUpdatePlot(self):    
        timeshift=self.ui.dsbTimeShift.value()*1e-12
        factor=self.ui.dsbAmplitudeFactor.value()
        where=self.ui.cbwhichGraphs.currentIndex()
        
        tempItem=self.createTemporaryCopy(self.ui.fileTree.topLevelItem(where))
        
        #manipulate the temporary data
        tdata=tempItem.tdData
        uefield=tdata.getUEfield()
        uefield*=factor
        tdata=TeraData.TimeDomainData(tdata.getTimeAxisRef()+timeshift,unumpy.nominal_values(uefield),unumpy.std_devs(uefield),tdata.getDataSetName())
        self.temporaryUpdate(tempItem,tdata)
        self.refreshCanvas()
    
    def applyTDManipulation(self):

        #make sure that the data of the last item is correct   
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)    
        if not lastitem.tempitem:
            self.dataManipulationTemporarilyUpdatePlot()
        
        where=self.ui.cbwhichGraphs.currentIndex()
        
        if self.ui.rbManipulateCopy.isChecked():
            #should be the temporary data 
            self.insertTemporaryCopy()
            self.ui.preferences.hide()
        else:
            self.insertTemporaryCopy(self.ui.fileTree.topLevelItem(where))
            self.ui.preferences.hide()
  
    def cancelPreferences(self):
        #also take back the preview in plots
        root=self.ui.fileTree.invisibleRootItem()
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)
        if lastitem.tempitem:        
            root.removeChild(lastitem)
        self.ui.preferences.hide()       
    
    def createTemporaryCopy(self,originalitem):
        #check if last item is already a temporary item
        #temporary items are not added to combo box
        #temporary items are by default grey colored
        if not self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1).tempitem:
            x=THzTreeWidgetItem()
            
            x.refreshCanvas=self.refreshCanvas
            x.setFlags(x.flags() | QtCore.Qt.ItemIsEditable)
            x.tdData=originalitem.tdData
            x.fdData=TeraData.FrequencyDomainData.fromTimeDomainData(originalitem.tdData)
            x.setCheckState(0,QtCore.Qt.Checked)
            x.setText(1,originalitem.text(1))
            x.color= plt.cm.colors.colorConverter.to_rgba('0.75')
            col=QtGui.QColor(int(x.color[0]*255),int(x.color[1]*255),int(x.color[2]*255),int(x.color[3]*255))
            x.setBackgroundColor(0,col)        
            x.tempitem=True            
            self.doTdFdPlot(x)
            self.ui.fileTree.addTopLevelItem(x)
            tempItem=x
            tempItem.setHidden(True)
        else:
            tempItem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)
        return tempItem
    
    def temporaryUpdate(self,tempItem,newtdData,newfdData=None):
        tempItem.tdData=newtdData
        if newfdData is None:        
            tempItem.fdData=TeraData.FrequencyDomainData.fromTimeDomainData(newtdData)
        self.doTdFdPlot(tempItem)

    def insertTemporaryCopy(self,originalItem=None):
        '''overwrites originalItem if given, else inserts a copy'''
        lastitem=self.ui.fileTree.topLevelItem(self.ui.fileTree.topLevelItemCount()-1)
        if originalItem is None:
            lastitem.color=plt.cm.brg(np.random.rand(1)[0])
            col=QtGui.QColor(int(lastitem.color[0]*255),int(lastitem.color[1]*255),int(lastitem.color[2]*255),int(lastitem.color[3]*255))
            lastitem.setBackgroundColor(0,col)
            lastitem.tempitem=False
            lastitem.setHidden(False)
            lastitem.setText(1,"Copy of " + lastitem.text(1))
            lastitem.tdData.setDataSetName("Copy of " + lastitem.text(1))
            self.updateDetails(lastitem)
            self.doTdFdPlot(lastitem)
        else:
            originalItem.tdData=lastitem.tdData
            originalItem.fdData=lastitem.fdData
            self.updateDetails(originalItem)
            self.ui.fileTree.invisibleRootItem().removeChild(lastitem)
            self.doTdFdPlot(originalItem)
    
    def refreshCanvas(self):
        #no legends so far
#        for ax in self.ui.spectrumCanvas.axes:
#            ax.relim()
#            ax.autoscale_view()
        self.ui.spectrumCanvas.draw()
    
    def removeCurrentlySelectedItem(self):
        msg=QtGui.QMessageBox()
        msg.setText("Remove Item permanently?")
        msg.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
        msg.setDefaultButton(QtGui.QMessageBox.Ok)
        
        if msg.exec_()==QtGui.QMessageBox.Ok:
            index=self.ui.fileTree.selectedIndexes()[0].row()
            self.ui.cbwhichGraphs.removeItem(index)
            root=self.ui.fileTree.invisibleRootItem()
            for item in self.ui.fileTree.selectedItems():
                (item.parent() or root).removeChild(item)

    def onCustomContextMenu(self,point):
        item=self.ui.fileTree.itemAt(point)
        if item != None and item.isSelected():
            if item.parent()==None:
                self.TopLevelItemMenu.exec_(self.ui.fileTree.mapToGlobal(point))
            else:
                self.ChildMenu.exec_(self.ui.fileTree.mapToGlobal(point))
        
    def initializeSpectrumCanvas(self):
        
        gs = gridspec.GridSpec(2, 2)
        
        #initialize the three main plots
        self.ui.spectrumCanvas.axes.append(self.ui.spectrumCanvas.figure.add_subplot(gs[0,:]))
        self.ui.spectrumCanvas.axes.append(self.ui.spectrumCanvas.figure.add_subplot(gs[1,0]))
        self.ui.spectrumCanvas.axes.append(self.ui.spectrumCanvas.figure.add_subplot(gs[1,1]))
      
        ax=self.ui.spectrumCanvas.axes
        
        #add twinned axes for subplots
        ax_SNRTD=ax[0].twinx()
        ax_SNRTD.axes.get_yaxis().set_visible(False)       
        self.ui.spectrumCanvas.axes.append(ax_SNRTD)
        
        ax_SNRFD=ax[1].twinx()
        ax_SNRFD.axes.get_yaxis().set_visible(False)       
        self.ui.spectrumCanvas.axes.append(ax_SNRFD)
                
        
        ax[0].set_xlabel('time in ps')
        ax[0].set_ylabel('Amplitude')        
        
        ax[1].set_xlabel('frequency in THz')
        ax[1].set_ylabel('Amplitude, dB Scale')        
        ax[1].set_xlim([0,10])
        ax[1].set_ylim([-90,0])

        ax[2].set_xlabel('frequency in THz')
        ax[2].set_ylabel('Phase')

        
    def plotSNR(self):
        self.plotspecial('SNR')        

    def plotDR(self):
        self.plotspecial('Dynamic Range')
    
    def plotuncertainty(self):
        self.plotspecial('uncertainty')

    def plotspecial(self,what):
        curvelist=[]
        for i in range(self.ui.fileTree.topLevelItemCount()):
            curvelist.append(self.ui.fileTree.topLevelItem(i).text(1))

        myplotdialog=PlotsToChange(curvelist,"Select Data for which you want the " + what + " plotted: ")
        if myplotdialog.exec_()==QtGui.QDialog.Rejected:
            return

        where=myplotdialog.getPlotsToChange()
        
        #first try to add just to the first entry a new child        
        if where>1:
            for i in range(self.ui.fileTree.topLevelItem(where-2).childCount()):
                if self.ui.fileTree.topLevelItem(where-2).child(i).text(1).split(" ")[-1]==what:
                    return 0
            self.addSubplot(self.ui.fileTree.topLevelItem(where-2),what)
        else:
            for row in range(self.ui.fileTree.topLevelItemCount()):
                if where==1 or self.ui.fileTree.topLevelItem(row).checkState(0):
                    found=False
                    for i in range(self.ui.fileTree.topLevelItem(row).childCount()):
                        if self.ui.fileTree.topLevelItem(row).child(i).text(1).split(" ")[-1]==what:
                            found=True
                            break
                    if not found:
                        self.addSubplot(self.ui.fileTree.topLevelItem(row),what)

    def addSubplot(self,tlw,what):
        x=THzTreeWidgetItem()
        x.setCheckState(0,QtCore.Qt.Checked)
        leg_label=tlw.text(1)+" " + what
        x.setText(1,leg_label)
        x.refreshCanvas=self.refreshCanvas
      
        if what=='SNR':
            ax_SNRTD=self.ui.spectrumCanvas.figure.axes[3]
            ax_SNRTD.axes.get_yaxis().set_visible(True)
            ax_SNRFD=self.ui.spectrumCanvas.figure.axes[4]
            ax_SNRFD.axes.get_yaxis().set_visible(True)
                        
            x.tdline=ax_SNRTD.plot(tlw.tdData.getTimeAxis()*1e12,tlw.tdData.getSNR(),'--',label=leg_label,color=tlw.color)
            x.fdlineabs=ax_SNRFD.plot(tlw.fdData.getFrequenciesRef()/1e12,20*np.log10(tlw.fdData.getSNR()),'--',label=leg_label,color=tlw.color)
#            x.fdlineabs=self.ui.spectrumCanvas.figure.axes[1].plot(tlw.fdData.getfreqsGHz()/1e3,20*np.log10(tlw.fdData.getSNR()),label=leg_label)
            x.setCheckState(0,QtCore.Qt.Checked)
            x.color=tlw.color
            x.setText(2,'test')
            x.setText(3,'test')
            pattern=QtCore.Qt.Dense3Pattern
        
        if what=='uncertainty':
            
            #plot absolute and phase along with uncertainties
            f=1
            #uabs=unumpy.uarray(tlw.fdData.getFAbs(),tlw.fdData.getFAbsUnc())
            #scale the uncertainty for the 20*log10 plot
            #uabs=20*unumpy.log10(uabs)
            #uabs-=np.amax(uabs)
            #u_a=unumpy.nominal_values(uabs)
            #u_s=unumpy.std_devs(uabs)
            x.tdline.append(self.ui.spectrumCanvas.figure.axes[0].plot(tlw.tdData.getTimeAxisRef()*1e12,tlw.tdData.getEfield()+tlw.tdData.getUncertainty(),linestyle='--',color=tlw.color)[0])
            x.tdline.append(self.ui.spectrumCanvas.figure.axes[0].plot(tlw.tdData.getTimeAxisRef()*1e12,tlw.tdData.getEfield()-tlw.tdData.getUncertainty(),linestyle='--',color=tlw.color,label=leg_label)[0])
     
            #x.fdlineabs.append(self.ui.spectrumCanvas.figure.axes[1].plot(tlw.fdData.getfreqsGHz()[1:]/1e3,u_a[1:]+u_s[1:],linestyle='--',color=tlw.color)[0])
            #.fdlineabs.append(self.ui.spectrumCanvas.figure.axes[1].plot(tlw.fdData.getfreqsGHz()[1:]/1e3,u_a[1:]-u_s[1:],linestyle='--',color=tlw.color,label=leg_label)[0])

            #x.fdlinephase.append(self.ui.spectrumCanvas.figure.axes[2].plot(tlw.fdData.getfreqsGHz()[1:]/1e3,tlw.fdData.getFPh()[1:]+tlw.fdData.getFPhUnc()[1:],linestyle='--',color=tlw.color)[0])
            #x.fdlinephase.append(self.ui.spectrumCanvas.figure.axes[2].plot(tlw.fdData.getfreqsGHz()[1:]/1e3,tlw.fdData.getFPh()[1:]-tlw.fdData.getFPhUnc()[1:],linestyle='--',color=tlw.color)[0])

            x.color=tlw.color
            x.setText(2,'test')
            x.setText(3,'test')
            pattern=QtCore.Qt.Dense3Pattern
        
        if what=='Dynamic Range':             
            ax_DRTD=self.ui.spectrumCanvas.figure.axes[3]
            ax_DRTD.axes.get_yaxis().set_visible(True)
            ax_DRFD=self.ui.spectrumCanvas.figure.axes[4]
            ax_DRFD.axes.get_yaxis().set_visible(True)
            
            x.tdline=ax_DRTD.plot(tlw.tdData.getTimeAxisRef()*1e12,tlw.tdData.getDynamicRange(),'-.',label=leg_label,color=tlw.color)
            x.fdlineabs=ax_DRFD.plot(tlw.fdData.getFrequenciesRef()/1e12,tlw.fdData.getDynamicRange(),'-.',label=leg_label,color=tlw.color)
            x.color=tlw.color            
            x.setText(2,'test')
            x.setText(3,'test')
            pattern=QtCore.Qt.Dense6Pattern
            
        col=QtGui.QColor(int(x.color[0]*255),int(x.color[1]*255),int(x.color[2]*255),int(x.color[3]*255))
        brush=QtGui.QBrush(col,pattern)        
        x.setBackground(0,brush)
        tlw.addChild(x)
        self.refreshCanvas()
        
    def windowing(self,data):
        if self.ui.cb_windowing.currentIndex()>0:
            start_window=self.ui.dsb_zeropaddfrom.value()
            end_window=self.ui.dsb_zeropaddfrom.value()
            
            
    def zeroPadding(self,dataset): 
              
        if self.ui.cb_zeropadding.isChecked():
            no_zeros=self.ui.sb_nozeros.value()
            dataset.tdData.zeroPaddData(no_zeros)
#            dataset.tdline[0].remove()
            dataset.fdData=TeraData.FdData(dataset.tdData)
            self.updateDetails(dataset)
            self.doTdFdPlot(dataset)
            print("zeropadding applied")
            
    def loadFile(self):
        filenames=QtGui.QFileDialog.getOpenFileNames()
        if len(filenames)==0:
            return 
            
        myformatdialog=FormatDialog()
        
        myformatdialog.setFilenames(filenames)
        self.ui.mainStatus.showMessage("Load Files")
        if myformatdialog.exec_()==QtGui.QDialog.Rejected:
            return
        fileformat=myformatdialog.getDataFormat()
        if myformatdialog.doAveraging():
            display_filename=path.split(str(filenames[0]))[1]
            for i in range(1,len(filenames)):
                display_filename+=" \n" + path.split(str(filenames[i]))[1]
            tdData=TeraData.TimeDomainData.importMultipleFiles(filenames,fileformat)
            tdData.setDataSetName(display_filename)
            x=self.fillTree(display_filename,tdData) 
        else:          
            for fn in filenames:
                display_filename=path.split(str(fn))[1]
                tdData=TeraData.TimeDomainData.fromFile(fn,fileformat)
                tdData=tdData._removeLinearDrift()
                tdData.setDataSetName(display_filename)
                x=self.fillTree(display_filename,tdData)
        self.ui.mainStatus.clearMessage()
        return filenames
        
    def fillTree(self,display_filename,tdData):
        x=THzTreeWidgetItem()
        x.refreshCanvas=self.refreshCanvas
        x.setFlags(x.flags() | QtCore.Qt.ItemIsEditable)
        x.tdData=tdData
        x.fdData=TeraData.FrequencyDomainData.fromTimeDomainData(tdData)
        x.setCheckState(0,QtCore.Qt.Checked)
        x.setText(1,display_filename)
        
        self.updateDetails(x)
        
        self.doTdFdPlot(x)
        self.ui.fileTree.addTopLevelItem(x)        
        col=QtGui.QColor(int(x.color[0]*255),int(x.color[1]*255),int(x.color[2]*255),int(x.color[3]*255))
        x.setBackgroundColor(0,col)
        self.ui.cbwhichGraphs.addItem(x.text(1))
        self.ui.cbwhichGraphs.setItemData(self.ui.cbwhichGraphs.count()-1,col,QtCore.Qt.BackgroundRole)
        return x        
    
    def updateSpectrumAnalysisPlot(self,item,column):
        #something else happened so return
        if column>1:
            return
        
        #the name of the plot changed, so change legend entry!
        #plot name not needed right now
        if column==1:
            
            #self.ui.cbwhichGraphs.setItemText(self.ui.fileTree.selectedIndexes()[0].row()+2,item.text(1))
            item.tdline[0].set_label(item.text(1))
            item.fdlineabs[0].set_label(item.text(1))
            item.fdlinephase[0].set_label(item.text(1))
            
            for i in range(item.childCount()):
                oldlabel=item.child(i).text(1)
                item.child(i).setText(1,item.text(1) +" " +oldlabel.split(" ")[-1])                
                item.child(i).tdline[0].set_label(item.text(1) +" " +oldlabel.split(" ")[-1])
                item.child(i).fdlineabs[0].set_label(item.text(1) + " " + oldlabel.split(" ")[-1])
                item.child(i).fdlinephase[0].set_label(item.text(1) + " " + oldlabel.split(" ")[-1])
                
      #the plots status changed
        if column ==0:
            if item.checkState(0):
                for line in item.tdline+item.fdlineabs+item.fdlinephase:
                    line.set_label(item.text(1))
                    line.set_visible(True)
                    line.axes.get_yaxis().set_visible(True)
            else:
                for line in item.tdline+item.fdlineabs+item.fdlinephase:
                    line.set_label(None)
                    line.set_visible(False)
                    found=False
                    #switch axes labeling off, if this was the last visible line                    
                    for l in line.axes.lines:
                        if l.get_visible()==True:
                            found=True
                            break
                    if not found:
                        line.axes.get_yaxis().set_visible(False)
                    
        self.refreshCanvas()
        
    def updateDetails(self,thztreeitem):
        d=thztreeitem.fdData.getBandwidth()
        thztreeitem.setText(2,'dt='+'{:3.2f}'.format(thztreeitem.tdData.getTimeStep()*1e15) + 'fs' +'\n'
                                'Pulsewidth=' + '{:3.2f}'.format(thztreeitem.tdData.getPeakWidth()*1e12) + 'ps')
        thztreeitem.setToolTip(2,'Pulse Position: ' + '{:3.3f}'.format(thztreeitem.tdData.getPeakPosition()*1e12) + ' ps\n'+ 
                                'DataLength=' + '{:d}'.format(thztreeitem.tdData.getSamplingPoints()) + 'Points\n' +
                                'TimeWindowLength=' + '{:3.2f}'.format((thztreeitem.tdData.getTimeAxisRef()[-1]-thztreeitem.tdData.getTimeAxisRef()[0])*1e12)+'Ps')
        
        thztreeitem.setText(3,'df=' +'{:3.2f}'.format(thztreeitem.fdData.getfbins()/1e9) + 'GHz' +'\n'
                                'Bandwidth='+'{:3.1f}'.format((d[1]-d[0])/1e12) + 'THz')
    
    def doTdFdPlot(self,thztreeitem):
        #in case no color yet defined
        if thztreeitem.color==None:
            thztreeitem.color=plt.cm.brg(np.random.rand(1)[0])
        
        #in case of update, remove existing lines from plot, 
        for line in thztreeitem.tdline:
            line.set_xdata(thztreeitem.tdData.getTimeAxisRef()*1e12)
            line.set_ydata(thztreeitem.tdData.getEfield())
            line.set_color(thztreeitem.color)
            line.set_label(thztreeitem.text(1))
            
        for line in thztreeitem.fdlineabs:
            
            absdata=20*np.log10(abs(thztreeitem.fdData.getSpectrumRef()))
            line.set_xdata(thztreeitem.fdData.getFrequenciesRef()*1e-12)
            line.set_ydata(absdata-np.amax(absdata))
            line.set_color(thztreeitem.color)
            line.set_label(thztreeitem.text(1))
            
        for line in thztreeitem.fdlinephase:
            
            line.set_xdata(thztreeitem.fdData.getFrequenciesRef()*1e-12)
            line.set_ydata(thztreeitem.fdData.getPhasesRef())
            line.set_color(thztreeitem.color)            
            line.set_label(thztreeitem.text(1))

        if len(thztreeitem.tdline)==0:
            
            thztreeitem.tdline=self.ui.spectrumCanvas.figure.axes[0].plot(thztreeitem.tdData.getTimeAxisRef()*1e12,thztreeitem.tdData.getEfield(),color=thztreeitem.color,label=thztreeitem.text(1))
            self.ui.spectrumCanvas.figure.axes[0].set_xlim([thztreeitem.tdData.getTimeAxisRef()[0]*1e12,thztreeitem.tdData.getTimeAxisRef()[-1]*1e12])
            

            absdata=20*np.log10(abs(thztreeitem.fdData.getSpectrumRef()))
            absdata-=np.amax(absdata)
            thztreeitem.fdlineabs=self.ui.spectrumCanvas.figure.axes[1].plot(thztreeitem.fdData.getFrequenciesRef()*1e-12,absdata,color=thztreeitem.color,label=thztreeitem.text(1))
            thztreeitem.fdlinephase=self.ui.spectrumCanvas.figure.axes[2].plot(thztreeitem.fdData.getFrequenciesRef()*1e-12,thztreeitem.fdData.getPhasesRef(),color=thztreeitem.color,label=thztreeitem.text(1))
        
        self.refreshCanvas()

    def onTreeWidgetItemDoubleClicked(self,item,column): 
        if column==1:
            self.ui.fileTree.editItem(item,column)
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = MyWindow()
    sys.exit(app.exec_())