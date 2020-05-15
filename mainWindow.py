from PyQt5.QtWidgets import QFileDialog, QMessageBox, QMainWindow
from PyQt5.QtCore import QThread, pyqtSlot, pyqtSignal

from Ui_mainWindow import Ui_MainWindow
from setting import SettingWindow
from about import AboutWindow

from AdNDP import AdNDP


class WorkThread(QThread):

    finishSignal = pyqtSignal()

    def __init__(self, AdNDP):
        super(WorkThread, self).__init__()
        self.AdNDP = AdNDP

    def run(self):
        self.AdNDP.partition()
        self.AdNDP.output()
        self.finishSignal.emit()


class Window(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(QMainWindow, self).__init__()
        self.setupUi(self)

        self.statusBar.showMessage("请选择文件")

        self.AdNDP = AdNDP()
        self.AdNDP.resultSignal.connect(self.resultTextEdit.append)
        self.AdNDP.logSignal.connect(self.logTextEdit.append)
        self.AdNDP.informationSignal.connect(self.informationTextEdit.append)
        self.AdNDP.setMaximumSignal.connect(self.progressBar.setMaximum)
        self.AdNDP.setValueSignal.connect(self.progressBar.setValue)

    @pyqtSlot()
    def on_openFileAct_triggered(self):
        path = QFileDialog.getOpenFileName(
            self, "请选择初始化文件", "./example", "NBO Files (*.log);;All Files (*)")[0]
        # path = "./example/Li5+.log"
        if path != '':
            self.AdNDP.readFile(path)
            self.startAnalysisAct.setEnabled(True)
            self.settingAct.setEnabled(True)
        else:
            QMessageBox.warning(self, "错误", "请选择文件")

    @pyqtSlot()
    def on_startAnalysisAct_triggered(self):
        self.statusBar.showMessage("分析中......")
        self.WorkThread = WorkThread(self.AdNDP)
        self.WorkThread.start()
        self.WorkThread.finishSignal.connect(self.on_workThread_finishSignal)

    def on_workThread_finishSignal(self):
        QMessageBox.about(self, "提示", "分析完成")
        self.exportFileAct.setEnabled(True)
        self.statusBar.showMessage(
            "分析结束，共找到" + str(self.AdNDP.nboAmnt) + "个轨道")

    @pyqtSlot()
    def on_aboutAct_triggered(self):
        self.aboutWindow = AboutWindow()
        self.aboutWindow.show()

    @pyqtSlot()
    def on_exportFileAct_triggered(self):
        savePath = QFileDialog.getSaveFileName(
            self, "请选择保存位置", "./example", "log格式 (*.log)")[0]
        if savePath != '':
            self.AdNDP.nboPlotMolden(savePath)
            QMessageBox.about(self, "提示", "保存成功")
        else:
            QMessageBox.warning(self, "提示", "请重新选择保存位置")

    @pyqtSlot()
    def on_settingAct_triggered(self):
        self.settingWindow = SettingWindow(self.AdNDP.NAt)
        self.settingWindow.resultSignal.connect(self.AdNDP.setThreshold)
        self.settingWindow.show()

    @pyqtSlot()
    def on_clearScreenAct_triggered(self):
        self.resultTextEdit.clear()
        self.informationTextEdit.clear()
        self.logTextEdit.clear()

        self.startAnalysisAct.setEnabled(False)
        self.exportFileAct.setEnabled(False)
        self.settingAct.setEnabled(False)
        self.statusBar.showMessage("清除完成，请重新选择文件......")


if __name__ == '__main__':
    import sys
    from PyQt5.QtWidgets import QApplication
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())
