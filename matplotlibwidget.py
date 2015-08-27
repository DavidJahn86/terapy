__version__="1.0.0"

from PyQt4.QtGui import QSizePolicy
from PyQt4.QtCore import QSize

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as Canvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure

from matplotlib import rcParams
rcParams['font.size']=9

class MatplotlibWidget(Canvas):
    #some help information here

    def __init__(self,parent=None,title='',xlabel='',ylabel='',xlim=None,
                 ylim=None,xscale='linear',yscale='linear',width=4,height=3,dpi=100,hold=False):
        self.figure=Figure(figsize=(width,height),dpi=dpi)
        self.axes=[]
#        
#        self.axes.set_title(title)
#        self.axes.set_xlabel(xlabel)
#        self.axes.set_ylabel(ylabel)
#        
#        if xscale is not None:
#            self.axes.set_xscale(xscale)
#        if yscale is not None:
#            self.axes.set_yscale(yscale)
#        if xlim is not None:
#            self.axes.set_xlim(*xlim)
#        if ylim is not None:
#            self.axes.set_ylim(*ylim)
##        self.axes.hold(hold)
        
        Canvas.__init__(self,self.figure)
        self.setParent(parent)
          
#        Canvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        Canvas.updateGeometry(self)
        
        self.navigationToolbar=NavigationToolbar(self,self)
        
    def sizeHint(self):
        w,h=self.get_width_height()
        return QSize(w,h)
        
    def minimumSizeHint(self):
        return QSize(10,10)
        
if __name__=='__main__':
    import sys
    from PyQt4.QtGui import QMainWindow, QApplication
    from numpy import linspace
    
    class ApplicationWindow(QMainWindow):
        def __init__(self):
            QMainWindow.__init__(self)
            self.mplwidget=MatplotlibWidget(self,title='Example',xlabel='Linear scale',
                                            ylabel='Log scale',hold=True,yscale='log')
            
            self.mplwidget.setFocus()
            self.setCentralWidget(self.mplwidget)
            self.mplwidget.axes.append(self.mplwidget.figure.add_subplot(111))
            self.plot(self.mplwidget.axes)
            x=linspace(-10,10)            
            newline=self.mplwidget.axes[0].plot(x,x**4)
            oldline=self.mplwidget.axes[0].lines[0]
            self.updateplot(self.mplwidget.axes[0],oldline,newline)
            
        def plot(self,axes):
            
            x=linspace(-10,10)
            axes[0].plot(x,x**2)
            axes[0].plot(x,x**3)
    
        def updateplot(self,axes,oldline,newline):
            oldline.remove()
            axes.add_line(newline)
        
    app=QApplication(sys.argv)
    win=ApplicationWindow()
    win.show()
    sys.exit(app.exec_())
    
    