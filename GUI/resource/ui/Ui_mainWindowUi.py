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
        MainWindow.resize(740, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.output_td = QtWidgets.QTextEdit(self.centralwidget)
        self.output_td.setObjectName("output_td")
        self.horizontalLayout.addWidget(self.output_td)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 740, 18))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuSetting = QtWidgets.QMenu(self.menubar)
        self.menuSetting.setObjectName("menuSetting")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.menu = QtWidgets.QMenu(self.menubar)
        self.menu.setObjectName("menu")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setEnabled(False)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.openFile_action = QtWidgets.QAction(MainWindow)
        self.openFile_action.setObjectName("openFile_action")
        self.exportFile_action = QtWidgets.QAction(MainWindow)
        self.exportFile_action.setObjectName("exportFile_action")
        self.clearScreen_action = QtWidgets.QAction(MainWindow)
        self.clearScreen_action.setObjectName("clearScreen_action")
        self.startAnalysis_action = QtWidgets.QAction(MainWindow)
        self.startAnalysis_action.setObjectName("startAnalysis_action")
        self.stopAnalysis_action = QtWidgets.QAction(MainWindow)
        self.stopAnalysis_action.setObjectName("stopAnalysis_action")
        self.about_action = QtWidgets.QAction(MainWindow)
        self.about_action.setObjectName("about_action")
        self.menuFile.addAction(self.openFile_action)
        self.menuFile.addAction(self.exportFile_action)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.clearScreen_action)
        self.menuHelp.addAction(self.about_action)
        self.menu.addAction(self.startAnalysis_action)
        self.menu.addAction(self.stopAnalysis_action)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menuSetting.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.openFile_action.triggered.connect(MainWindow.openFile)
        self.startAnalysis_action.triggered.connect(MainWindow.startAnalysis)
        self.clearScreen_action.triggered.connect(MainWindow.clearScreen)
        self.exportFile_action.triggered.connect(MainWindow.exportFile)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "AdNDP"))
        self.menuFile.setTitle(_translate("MainWindow", "文件"))
        self.menuSetting.setTitle(_translate("MainWindow", "设置"))
        self.menuHelp.setTitle(_translate("MainWindow", "帮助"))
        self.menu.setTitle(_translate("MainWindow", "分析"))
        self.openFile_action.setText(_translate("MainWindow", "打开"))
        self.exportFile_action.setText(_translate("MainWindow", "导出"))
        self.clearScreen_action.setText(_translate("MainWindow", "清除"))
        self.startAnalysis_action.setText(_translate("MainWindow", "开始分析"))
        self.stopAnalysis_action.setText(_translate("MainWindow", "结束分析"))
        self.about_action.setText(_translate("MainWindow", "关于"))
