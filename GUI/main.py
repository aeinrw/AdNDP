from PyQt5.Qt import QWidget, QFileDialog, QApplication, QProgressBar, QMessageBox, QMainWindow
from resource.ui.Ui_mainWindowUi import Ui_MainWindow
from AdNDP import AdNDP


class Window(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(QMainWindow, self).__init__()
        self.setupUi(self)

    def openFile(self):
        path = QFileDialog.getOpenFileName(
            self, "请选择初始化文件", "../example", "Init Files (*.ini);;All Files (*)")[0]
        if path != '':
            self.AdNDP = AdNDP(path, self.output_td)
        else:
            QMessageBox.critical(self, "Error!", "Please pick a file!")

    def startAnalysis(self):
        self.AdNDP.partition()
        self.AdNDP.output()

    def exportFile(self):
        pass

    def clearScreen(self):
        del self.AdNDP
        self.output_td.clear()


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())
