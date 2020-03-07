# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'c:\Source\AdNDP\GUI\resource\ui\mainWindowUi.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.readInitFile_btn = QtWidgets.QPushButton(self.centralwidget)
        self.readInitFile_btn.setObjectName("readInitFile_btn")
        self.verticalLayout.addWidget(self.readInitFile_btn)
        self.analysis_btn = QtWidgets.QPushButton(self.centralwidget)
        self.analysis_btn.setObjectName("analysis_btn")
        self.verticalLayout.addWidget(self.analysis_btn)
        self.exportFile_btn = QtWidgets.QPushButton(self.centralwidget)
        self.exportFile_btn.setObjectName("exportFile_btn")
        self.verticalLayout.addWidget(self.exportFile_btn)
        self.clear_btn = QtWidgets.QPushButton(self.centralwidget)
        self.clear_btn.setObjectName("clear_btn")
        self.verticalLayout.addWidget(self.clear_btn)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.output_td = QtWidgets.QTextEdit(self.centralwidget)
        self.output_td.setObjectName("output_td")
        self.horizontalLayout.addWidget(self.output_td)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.readInitFile_btn.clicked.connect(MainWindow.readInitFile)
        self.analysis_btn.clicked.connect(MainWindow.adndpAnalysis)
        self.exportFile_btn.clicked.connect(MainWindow.exportOutputFile)
        self.clear_btn.clicked.connect(MainWindow.clear)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.readInitFile_btn.setText(_translate("MainWindow", "选择NBO文件"))
        self.analysis_btn.setText(_translate("MainWindow", "进行AdNDP分析"))
        self.exportFile_btn.setText(_translate("MainWindow", "导出文件"))
        self.clear_btn.setText(_translate("MainWindow", "清除"))
