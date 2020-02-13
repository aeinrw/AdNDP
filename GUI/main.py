from PyQt5.Qt import QWidget,QFileDialog,QApplication,QProgressBar,QMessageBox
from resource.ui.Ui_mainWindow import Ui_Form
from AdNDP import AdNDP

class Window(QWidget,Ui_Form):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

    def readInitFile(self):
        path = QFileDialog.getOpenFileName(self,"请选择初始化文件","../example","Init Files (*.ini);;All Files (*)")[0]
        #path = 'C:/Source/Graduation Design/AdNDP/AdNDP/GUI/resource/profile/AdNDP.ini'
        if path != '':
            self.AdNDP = AdNDP(path,self.output_td)
            self.readInitFile_btn.setEnabled(False)
            self.analysis_btn.setEnabled(True)
        else:
            QMessageBox.critical(self,"Error!","Please pick a file!")

    def adndpAnalysis(self):
        self.analysis_btn.setEnabled(False)
        self.AdNDP.partition()
        self.AdNDP.output()
        self.clear_btn.setEnabled(True)

    def exportOutputFile(self):
        # self.exportFile_btn.setEnabled(False)
        # self.readInitFile_btn.setEnabled(True)
        pass

    def clear(self):
        del self.AdNDP
        self.output_td.clear()
        self.readInitFile_btn.setEnabled(True)
        self.clear_btn.setEnabled(False)





if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())