from Ui_setting import Ui_Form
from PyQt5.QtWidgets import QWidget, QCheckBox, QHBoxLayout, QDoubleSpinBox, QSizePolicy
from PyQt5.QtCore import pyqtSlot, pyqtSignal, Qt, QObject


class MyWidget(QHBoxLayout):
    def __init__(self, i):
        super(MyWidget, self).__init__()

        self.checkBox = QCheckBox('{:2d}c-2e'.format(i + 1))
        self.doubleSpinBox = QDoubleSpinBox()
        self.addWidget(self.checkBox)
        self.addWidget(self.doubleSpinBox)
        self.set()

    def set(self):

        self.checkBox.setChecked(True)
        self.doubleSpinBox.setMaximum(2.0)
        self.doubleSpinBox.setMinimum(0.0)
        self.doubleSpinBox.setSingleStep(0.1)
        self.doubleSpinBox.setValue(0.1)
        self.checkBox.clicked['bool'].connect(self.doubleSpinBox.setEnabled)


class SettingWindow(QWidget, Ui_Form):

    resultSignal = pyqtSignal(list)

    def __init__(self, n):
        super(SettingWindow, self).__init__()
        self.setupUi(self)

        self.n = n
        self.verticalLayout.setAlignment(Qt.AlignCenter)
        self.widgetList = []
        for i in range(self.n):
            self.widgetList.append(MyWidget(i))
            self.verticalLayout.addLayout(self.widgetList[i])

    @pyqtSlot()
    def on_pushButton_clicked(self):
        result = []
        for i in range(self.n):
            flag = self.widgetList[i].checkBox.isChecked()
            if flag == True:
                result.append(self.widgetList[i].doubleSpinBox.value())
            else:
                result.append(0)
        self.resultSignal.emit(result)


if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    app = QApplication(sys.argv)
    window = SettingWindow(12)
    window.show()
    sys.exit(app.exec())
