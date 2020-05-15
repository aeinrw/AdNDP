from Ui_setting import Ui_Form
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QCheckBox
from PyQt5.QtCore import pyqtSlot, pyqtSignal, Qt


class SettingWindow(QWidget, Ui_Form):

    resultSignal = pyqtSignal(list)

    def __init__(self, n):
        super(SettingWindow, self).__init__()
        self.setupUi(self)

        self.n = n
        self.verticalLayout.setAlignment(Qt.AlignCenter)
        self.checkBoxList = []
        for i in range(self.n):
            self.checkBoxList.append(
                QCheckBox('{:d}c-2e'.format(i+1), self.groupBox))
            self.checkBoxList[i].setChecked(True)
            self.verticalLayout.addWidget(self.checkBoxList[i])

    @pyqtSlot()
    def on_pushButton_clicked(self):
        result = []
        for i in range(self.n):
            result.append(not self.checkBoxList[i].isChecked())
        print(result)
        self.resultSignal.emit(result)


if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    app = QApplication(sys.argv)
    window = SettingWindow(12)
    window.show()
    sys.exit(app.exec())
