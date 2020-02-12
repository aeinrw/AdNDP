# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(546, 327)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.btn_selectNboFile = QtWidgets.QPushButton(self.centralwidget)
        self.btn_selectNboFile.setGeometry(QtCore.QRect(60, 80, 56, 17))
        self.btn_selectNboFile.setObjectName("btn_selectNboFile")
        self.textEdit_output = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_output.setGeometry(QtCore.QRect(180, 50, 201, 171))
        self.textEdit_output.setObjectName("textEdit_output")
        self.btn_AdNDP = QtWidgets.QPushButton(self.centralwidget)
        self.btn_AdNDP.setGeometry(QtCore.QRect(50, 120, 81, 17))
        self.btn_AdNDP.setObjectName("btn_AdNDP")
        self.btn_exportFile = QtWidgets.QPushButton(self.centralwidget)
        self.btn_exportFile.setGeometry(QtCore.QRect(60, 160, 56, 17))
        self.btn_exportFile.setObjectName("btn_exportFile")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 546, 18))
        self.menubar.setObjectName("menubar")
        self.menu = QtWidgets.QMenu(self.menubar)
        self.menu.setObjectName("menu")
        self.menu_2 = QtWidgets.QMenu(self.menubar)
        self.menu_2.setObjectName("menu_2")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menu_2.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "AdNDP"))
        self.btn_selectNboFile.setText(_translate("MainWindow", "选择NBO文件"))
        self.btn_AdNDP.setText(_translate("MainWindow", "进行AdNDP分析"))
        self.btn_exportFile.setText(_translate("MainWindow", "导出文件"))
        self.menu.setTitle(_translate("MainWindow", "文件"))
        self.menu_2.setTitle(_translate("MainWindow", "保存"))
