# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'c:\Source\Graduation Design\AdNDP\AdNDP\GUI\resource\ui\mainWindow.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(566, 357)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setMinimumSize(QtCore.QSize(546, 291))
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.readInitFile_btn = QtWidgets.QPushButton(self.centralwidget)
        self.readInitFile_btn.setObjectName("readInitFile_btn")
        self.verticalLayout.addWidget(self.readInitFile_btn)
        self.analysis_btn = QtWidgets.QPushButton(self.centralwidget)
        self.analysis_btn.setObjectName("analysis_btn")
        self.verticalLayout.addWidget(self.analysis_btn)
        self.exportFile_btn = QtWidgets.QPushButton(self.centralwidget)
        self.exportFile_btn.setObjectName("exportFile_btn")
        self.verticalLayout.addWidget(self.exportFile_btn)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.output_td = QtWidgets.QTextEdit(self.centralwidget)
        self.output_td.setMaximumSize(QtCore.QSize(460, 179))
        self.output_td.setObjectName("output_td")
        self.horizontalLayout.addWidget(self.output_td)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 566, 18))
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
        self.readInitFile_btn.pressed.connect(MainWindow.readInitFile)
        self.analysis_btn.pressed.connect(MainWindow.analysis)
        self.exportFile_btn.pressed.connect(MainWindow.exportOutputFile)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "AdNDP"))
        self.readInitFile_btn.setText(_translate("MainWindow", "选择NBO文件"))
        self.analysis_btn.setText(_translate("MainWindow", "进行AdNDP分析"))
        self.exportFile_btn.setText(_translate("MainWindow", "导出文件"))
        self.menu.setTitle(_translate("MainWindow", "文件"))
        self.menu_2.setTitle(_translate("MainWindow", "保存"))
