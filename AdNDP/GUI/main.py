from PyQt5.Qt import *
from resource.ui.Ui_mainWindow import Ui_MainWindow
from AdNDP import AdNDP

class Window(QMainWindow,Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

    def readInitFile(self):
        path = QFileDialog.getOpenFileName(self,"请选择初始化文件","./resource/profile","Init Files (*.ini);;All Files (*)")[0]
        #path = 'C:/Source/Graduation Design/AdNDP/AdNDP/GUI/resource/profile/AdNDP.ini'
        self.AdNDP = AdNDP(path,self.output_td)

    def analysis(self):
        self.AdNDP.partition()
        self.AdNDP.output()


    def exportOutputFile(self):
        pass




if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())