# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'c:\Source\AdNDP\GUI\resource\ui\mainWindow.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.setEnabled(True)
        Form.resize(600, 400)
        Form.setMinimumSize(QtCore.QSize(600, 400))
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(Form)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.readInitFile_btn = QtWidgets.QPushButton(Form)
        self.readInitFile_btn.setEnabled(True)
        self.readInitFile_btn.setObjectName("readInitFile_btn")
        self.verticalLayout.addWidget(self.readInitFile_btn)
        self.analysis_btn = QtWidgets.QPushButton(Form)
        self.analysis_btn.setEnabled(False)
        self.analysis_btn.setObjectName("analysis_btn")
        self.verticalLayout.addWidget(self.analysis_btn)
        self.exportFile_btn = QtWidgets.QPushButton(Form)
        self.exportFile_btn.setEnabled(False)
        self.exportFile_btn.setObjectName("exportFile_btn")
        self.verticalLayout.addWidget(self.exportFile_btn)
        self.clear_btn = QtWidgets.QPushButton(Form)
        self.clear_btn.setEnabled(False)
        self.clear_btn.setObjectName("clear_btn")
        self.verticalLayout.addWidget(self.clear_btn)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.output_td = QtWidgets.QTextEdit(Form)
        self.output_td.setEnabled(True)
        self.output_td.setReadOnly(True)
        self.output_td.setObjectName("output_td")
        self.horizontalLayout.addWidget(self.output_td)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(Form)
        self.readInitFile_btn.pressed.connect(Form.readInitFile)
        self.exportFile_btn.clicked.connect(Form.exportOutputFile)
        self.analysis_btn.clicked.connect(Form.adndpAnalysis)
        self.clear_btn.clicked.connect(Form.clear)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "AdNDP分析程序"))
        self.readInitFile_btn.setText(_translate("Form", "选择NBO文件"))
        self.analysis_btn.setText(_translate("Form", "进行AdNDP分析"))
        self.exportFile_btn.setText(_translate("Form", "导出文件"))
        self.clear_btn.setText(_translate("Form", "清除"))
