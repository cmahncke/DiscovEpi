#################################################
# Cedric Mahncke
# cedric.mahncke@leibniz-liv.de
# DiscovEpi
# Version 1.1
# Automatically retrieve protein data,
# predict corresponding epitopes and produce
# an epitope map for whole proteomes.
# 08.07.2024: Product.
#################################################

import os
import sys
import re
from pathlib import Path
import pandas as pd
import webbrowser
from PySide6.QtWidgets import *
from PySide6.QtCore import QObject, QThread, Signal, Slot
from PySide6.QtGui import QIcon
import uniProt_gui as Unp
import NetMHCpan_gui as Nmp
import addtional_guis as add
import heatmap_gui as Hmp

#import faulthandler
#faulthandler.enable()


class MainWindow4(QMainWindow):
    ustart = Signal()
    nstart = Signal()
    hstart = Signal()
    uabort = Signal()
    nabort = Signal()
    habort = Signal()
    uempty = Signal()
    nempty = Signal()
    hempty = Signal()
    unocon = Signal()
    nnocon = Signal()

    def __init__(self):
        super(MainWindow4, self).__init__()
        cur_dir = os.path.dirname(os.path.realpath(__file__))
        self.protein_data = None
        self.prediction_data = None
        self.heatmap_data = None
        self.del_redundant = None
        self.progress_value = None
        self.UNPprogressBar = QProgressBar(self)
        self.NMPprogressBar = QProgressBar(self)
        self.HMPprogressBar = QProgressBar(self)

        self.setWindowTitle("DiscovEpi")

        self.tabs = QTabWidget(self)
        self.tabs.setDocumentMode(True)

        ugrid = QGridLayout(self)
        ugrid.addWidget(QLabel("Organism"), 1, 0)
        ugrid.addWidget(QLabel("Location"), 2, 0)
        ugrid.addWidget(QLabel("Only Reviewed"), 3, 0)
        ugrid.addWidget(QLabel("Progress"), 6, 0)
        ugrid.addWidget(QLabel("Directory"), 8, 0)

        self.org = QComboBox()
        self.org.setEditable(True)
        self.org.setMaximumWidth(550)
        org_lst = pd.read_table(cur_dir+'/data/orgs.txt', header=None)
        self.org.addItems(org_lst[0])
        self.org.setEditText("SARS-CoV-2")
        self.loc = QComboBox()
        self.loc.setEditable(True)
        self.loc.setMaximumWidth(550)
        loc_lst = pd.read_table(cur_dir+'/data/locs.txt', header=None)
        self.loc.addItems(loc_lst[0])
        self.loc.setCurrentText("Membrane")
        self.rev = QCheckBox("")
        self.rev.setChecked(True)
        self.ulog = QLabel("Logs:")
        self.umsg = QLineEdit()
        self.umsg.setReadOnly(True)

        uBtnLayout = QHBoxLayout()
        self.nmpBtn = QPushButton("Continue to NetMHCpan")
        self.nmpBtn.setEnabled(False)
        self.subUP = QPushButton("Submit")
        self.udirBtn = QPushButton("Go to Directory")
        self.udirBtn.setEnabled(False)
        self.uhelpBtn = QPushButton("Help")
        self.uhelpBtn.setIcon(QIcon(cur_dir+'/data/help.png'))
        self.uCnlBtn = QPushButton("Cancel")
        uBtnLayout.addWidget(self.subUP)
        uBtnLayout.addWidget(self.nmpBtn)
        uBtnLayout.addWidget(self.udirBtn)
        uBtnLayout.addWidget(self.uhelpBtn)
        uBtnLayout.addWidget(self.uCnlBtn)

        verticalSpacer = QSpacerItem(300, 20, QSizePolicy.Expanding)
        self.unpdir = QLineEdit(add.get_folder() + "/Epitope_Prediction/", self)
        ugrid.addWidget(self.org, 1, 1)
        ugrid.addWidget(self.loc, 2, 1)
        ugrid.addWidget(self.rev, 3, 1)
        ugrid.addLayout(uBtnLayout, 5, 1)
        ugrid.addWidget(self.UNPprogressBar, 6, 1)
        ugrid.addItem(verticalSpacer, 7, 1)
        ugrid.addWidget(self.unpdir, 8, 1)
        ugrid.addWidget(self.ulog, 9, 0)
        ugrid.addWidget(self.umsg, 9, 1)

        ngrid = QGridLayout(self)
        ngrid.addWidget(QLabel("Allele"), 1, 0)
        ngrid.addWidget(QLabel("Epitope Length [aa]"), 2, 0)
        ngrid.addWidget(QLabel("Threshold"), 3, 0)
        ngrid.addWidget(QLabel("Progress"), 5, 0)
        ngrid.addWidget(QLabel("Directory"), 7, 0)

        self.allele = QComboBox()
        self.allele.setEditable(True)
        mhc_lst = pd.read_table(cur_dir+'/data/mhcs.txt', header=None)
        self.allele.addItems(mhc_lst[0])
        self.allele.setCurrentText("HLA-A*02:01")
        self.len = QComboBox()
        self.len.setEditable(True)
        self.len.addItems(["8", "9", "10", "11", "12", "13", "14"])
        self.len.setCurrentText("9")
        self.thd = QLineEdit("3.0", self)
        self.nlog = QLabel("Logs:")
        self.nmsg = QLineEdit()
        self.nmsg.setReadOnly(True)

        nBtnLayout = QHBoxLayout()
        self.ndirBtn = QPushButton("Go to Directory")
        self.subNMP = QPushButton("Submit")
        self.subNMP.setEnabled(False)
        self.hmpBtn = QPushButton("Continue to epitope map")
        self.hmpBtn.setEnabled(False)
        self.nhelpBtn = QPushButton("Help")
        self.nhelpBtn.setIcon(QIcon(cur_dir+'/data/help.png'))
        self.nCnlBtn = QPushButton("Cancel")
        nBtnLayout.addWidget(self.subNMP)
        nBtnLayout.addWidget(self.hmpBtn)
        nBtnLayout.addWidget(self.ndirBtn)
        nBtnLayout.addWidget(self.nhelpBtn)
        nBtnLayout.addWidget(self.nCnlBtn)

        self.nmpdir = QLineEdit(add.get_folder() + "/Epitope_Prediction/", self)
        ngrid.addWidget(self.allele, 1, 1)
        ngrid.addWidget(self.len, 2, 1)
        ngrid.addWidget(self.thd, 3, 1)
        ngrid.addLayout(nBtnLayout, 4, 1)
        ngrid.addWidget(self.NMPprogressBar, 5, 1)
        ngrid.addItem(verticalSpacer, 6, 1)
        ngrid.addWidget(self.nmpdir, 7, 1)
        ngrid.addWidget(self.nlog, 8, 0)
        ngrid.addWidget(self.nmsg, 8, 1)

        hgrid = QGridLayout(self)
        hgrid.addWidget(QLabel("Number of Proteins: "), 1, 0)
        hgrid.addWidget(QLabel("Protein Cutoff [aa]"), 2, 0)
        hgrid.addWidget(QLabel("Without cutoff the longest protein will determine the range of amino acids to be shown.\n"
                               "To increase resolution on shorter proteins set a length cutoff."), 3, 1, 1, 2)
        hgrid.addWidget(QLabel("Progress"), 5, 0)
        hgrid.addWidget(QLabel("Directory"), 7, 0)
        self.top = QLineEdit("50", self)
        self.hco = QLineEdit("1000", self)
        self.hlog = QLabel("Logs:")
        self.hmsg = QLineEdit()
        self.hmsg.setReadOnly(True)

        hBtnLayout = QHBoxLayout()
        self.hdirBtn = QPushButton("Go to Directory")
        self.hdirBtn.setEnabled(False)
        self.subHMP = QPushButton("Submit")
        self.subHMP.setEnabled(False)
        self.hhelpBtn = QPushButton("Help")
        self.hhelpBtn.setIcon(QIcon(cur_dir+'/data/help.png'))
        self.hCnlBtn = QPushButton("Cancel")
        hBtnLayout.addWidget(self.subHMP)
        hBtnLayout.addWidget(self.hdirBtn)
        hBtnLayout.addWidget(self.hhelpBtn)
        hBtnLayout.addWidget(self.hCnlBtn)

        self.hmpdir = QLineEdit(add.get_folder() + "/Epitope_Prediction/", self)
        hgrid.addWidget(self.top, 1, 1)
        hgrid.addWidget(self.hco, 2, 1)
        hgrid.addLayout(hBtnLayout, 4, 1, 1, 2)
        hgrid.addWidget(self.HMPprogressBar, 5, 1)
        hgrid.addItem(verticalSpacer, 6, 1)
        hgrid.addWidget(self.hmpdir, 7, 1)
        hgrid.addWidget(self.hlog, 8, 0)
        hgrid.addWidget(self.hmsg, 8, 1)

        uwidget = QWidget(self)
        uwidget.setLayout(ugrid)
        nwidget = QWidget(self)
        nwidget.setLayout(ngrid)
        hwidget = QWidget(self)
        hwidget.setLayout(hgrid)
        self.tabs.addTab(uwidget, "UniProt")
        self.tabs.addTab(nwidget, "NetMHCpan")
        self.tabs.addTab(hwidget, "Epitope Map")

        self.setCentralWidget(self.tabs)

        self.nmpBtn.clicked.connect(self.on_nmpBtn_clicked)
        self.hmpBtn.clicked.connect(self.on_hmpBtn_clicked)
        self.subUP.clicked.connect(self.on_submitUP_clicked)
        self.subNMP.clicked.connect(self.on_submitNMP_clicked)
        self.subHMP.clicked.connect(self.on_submitHMP_clicked)
        self.unpdir.editingFinished.connect(self.change_directory)
        self.nmpdir.editingFinished.connect(self.change_directory)
        self.hmpdir.editingFinished.connect(self.change_directory)
        self.udirBtn.clicked.connect(self.on_directoryButton_clicked)
        self.ndirBtn.clicked.connect(self.on_directoryButton_clicked)
        self.hdirBtn.clicked.connect(self.on_directoryButton_clicked)
        self.uhelpBtn.clicked.connect(self.on_helpButton_clicked)
        self.nhelpBtn.clicked.connect(self.on_helpButton_clicked)
        self.hhelpBtn.clicked.connect(self.on_helpButton_clicked)

        self.unp = unpWorker(self)
        self.nmp = nmpWorker(self)
        self.hmp = hmpWorker(self)
        self.ustart.connect(self.unp_started)
        self.nstart.connect(self.nmp_started)
        self.hstart.connect(self.hmp_started)
        self.unocon.connect(self.unp_nocon)
        self.nnocon.connect(self.nmp_nocon)
        self.uempty.connect(self.unp_empty)
        self.nempty.connect(self.nmp_empty)
        self.hempty.connect(self.hmp_empty)
        self.uabort.connect(self.unp_aborted)
        self.nabort.connect(self.nmp_aborted)
        self.habort.connect(self.hmp_aborted)
        self.unp.finished.connect(self.unp_finished)
        self.nmp.finished.connect(self.nmp_finished)
        self.hmp.finished.connect(self.hmp_finished)
        self.uCnlBtn.clicked.connect(self.unp.requestInterruption)
        self.nCnlBtn.clicked.connect(self.nmp.requestInterruption)
        self.hCnlBtn.clicked.connect(self.hmp.requestInterruption)

    def on_nmpBtn_clicked(self):
        self.tabs.setCurrentIndex(1)

    def on_hmpBtn_clicked(self):
        self.tabs.setCurrentIndex(2)

    def on_submitUP_clicked(self):
        del_red_dialog = add.DeleteRedundantDialog(self)
        del_red_dialog.exec()

        if del_red_dialog.answer == 1:
            self.del_redundant = True
        elif del_red_dialog.answer == 0:
            self.del_redundant = False
        elif del_red_dialog.answer == -1:
            return

        Path(self.unpdir.text()).mkdir(parents=True, exist_ok=True)

        self.unp.start()

    def on_submitNMP_clicked(self):

        Path(self.nmpdir.text()).mkdir(parents=True, exist_ok=True)
        self.nmp.start()

    def change_directory(self):
        if self.sender() == self.unpdir:
            self.nmpdir.setText(self.unpdir.text())
            self.hmpdir.setText(self.unpdir.text())
        elif self.sender() == self.nmpdir:
            self.unpdir.setText(self.nmpdir.text())
            self.hmpdir.setText(self.nmpdir.text())
        elif self.sender() == self.hmpdir:
            self.unpdir.setText(self.hmpdir.text())
            self.nmpdir.setText(self.hmpdir.text())

    def on_submitHMP_clicked(self):

        Path(self.hmpdir.text()).mkdir(parents=True, exist_ok=True)

        self.hmp.start()

    def on_directoryButton_clicked(self):
        try:
            webbrowser.open('file:///' + self.unpdir.text())
        except:
            try:
                webbrowser.open('file:///' + self.nmpdir.text())
            except:
                try:
                    webbrowser.open('file:///' + self.hmpdir.text())
                except:
                    print("Error: No directory found.")

    def on_helpButton_clicked(self):
        cur_dir = os.path.dirname(os.path.realpath(__file__))
        webbrowser.open(cur_dir+'/data/README.md')

    def unp_started(self):
        self.umsg.setText("UniProt data retrieval in progress.")
        self.umsg.setStyleSheet("")
        self.udirBtn.setEnabled(False)
        self.nmpBtn.setEnabled(False)
        self.subNMP.setEnabled(False)
        self.update_UNPprogress(0)
        print("UniProt started")

    def nmp_started(self):
        self.nmsg.setText("NetMHCpan prediction in progress.")
        self.nmsg.setStyleSheet("")
        self.hmpBtn.setEnabled(False)
        self.ndirBtn.setEnabled(False)
        self.subHMP.setEnabled(False)
        self.update_NMPprogress(0)
        print("NetMHCpan started")

    def hmp_started(self):
        self.hmsg.setText("Epitope map creation in progress.")
        self.hmsg.setStyleSheet("")
        self.hdirBtn.setEnabled(False)
        self.update_HMPprogress(0)
        print("Epitope map started")

    def unp_finished(self):
        self.update_UNPprogress(100)
        self.umsg.setText("UniProt data successfully retrieved.")
        self.umsg.setStyleSheet("color: green")
        self.udirBtn.setEnabled(True)
        self.nmpBtn.setEnabled(True)
        self.subNMP.setEnabled(True)
        print("UniProt finished")

    def nmp_finished(self):
        self.update_NMPprogress(100)
        self.nmsg.setText("NetMHCpan prediction successfully completed.")
        self.nmsg.setStyleSheet("color: green")
        self.hmpBtn.setEnabled(True)
        self.ndirBtn.setEnabled(True)
        self.subHMP.setEnabled(True)
        print("NetMHCpan finished")

    def hmp_finished(self):
        self.update_HMPprogress(100)
        self.hmsg.setText("Epitope map successfully created.")
        self.hmsg.setStyleSheet("color: green")
        self.hdirBtn.setEnabled(True)
        print("Epitope maps finished")

    def unp_aborted(self):
        self.umsg.setText("UniProt data retrieval aborted.")
        self.umsg.setStyleSheet("color: red")
        self.udirBtn.setEnabled(False)
        self.nmpBtn.setEnabled(False)
        self.subNMP.setEnabled(False)
        self.update_UNPprogress(0)
        print("UniProt aborted")

    def nmp_aborted(self):
        self.nmsg.setText("NetMHCpan prediction aborted.")
        self.nmsg.setStyleSheet("color: red")
        self.hmpBtn.setEnabled(False)
        self.ndirBtn.setEnabled(False)
        self.subHMP.setEnabled(False)
        self.update_NMPprogress(0)
        print("NetMHCpan aborted")

    def hmp_aborted(self):
        self.hmsg.setText("Epitope map creation aborted.")
        self.hmsg.setStyleSheet("color: red")
        self.hdirBtn.setEnabled(False)
        self.update_HMPprogress(0)
        print("Epitope maps aborted")

    def unp_empty(self):
        self.umsg.setText("Dataset is empty. Check your input and try again.")
        self.umsg.setStyleSheet("color: red")
        self.udirBtn.setEnabled(False)
        self.nmpBtn.setEnabled(False)
        self.subNMP.setEnabled(False)
        self.update_UNPprogress(0)
        print("UniProt empty")

    def nmp_empty(self):
        self.nmsg.setText("Dataset is empty. Check your input and try again.")
        self.nmsg.setStyleSheet("color: red")
        self.hmpBtn.setEnabled(False)
        self.ndirBtn.setEnabled(False)
        self.subHMP.setEnabled(False)
        self.update_NMPprogress(0)
        print("NetMHCpan empty")

    def hmp_empty(self):
        self.hmsg.setText("Error in epitope map. Check your input and try again.")
        self.hmsg.setStyleSheet("color: red")
        self.hdirBtn.setEnabled(False)
        self.update_HMPprogress(0)
        print("Epitope map empty")

    def unp_nocon(self):
        self.umsg.setText("Connection error. Check your network and try again.")
        self.umsg.setStyleSheet("color: red")
        self.udirBtn.setEnabled(False)
        self.nmpBtn.setEnabled(False)
        self.subNMP.setEnabled(False)
        self.update_UNPprogress(0)
        print("UniProt connection error")

    def nmp_nocon(self):
        self.nmsg.setText("Connection error. Check your network and try again.")
        self.nmsg.setStyleSheet("color: red")
        self.hmpBtn.setEnabled(False)
        self.ndirBtn.setEnabled(False)
        self.subHMP.setEnabled(False)
        self.update_NMPprogress(0)
        print("NetMHCpan connection error")

    @Slot(int)
    def update_UNPprogress(self, value):
        self.progress_value = value
        self.UNPprogressBar.setValue(value)
        print(value)

    @Slot(int)
    def update_NMPprogress(self, value):
        self.progress_value = value
        self.NMPprogressBar.setValue(value)
        print(value)

    @Slot(int)
    def update_HMPprogress(self, value):
        self.progress_value = value
        self.HMPprogressBar.setValue(value)


class MySignals(QObject):
    signal_int = Signal(int)


class unpWorker(QThread):
    def __init__(self, parent=QThread):
        QThread.__init__(self, parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_UNPprogress)
        self.result = None
        self.org = parent.org
        self.loc = parent.loc
        self.rev = parent.rev
        self.del_redundant = parent.del_redundant
        self.dir = parent.unpdir
        self.parent = parent
        self.starting = parent.ustart
        self.abort = parent.uabort
        self.empty = parent.uempty
        self.nocon = parent.unocon


    def run(self):
        self.starting.emit()

        self.blockSignals(True)

        self.parent.protein_data = Unp.exec_uniprot([re.sub(r"[\t\n]", "", self.org.currentText()), self.loc.currentText(),
                                                     self.rev.isChecked()], self.del_redundant, self.dir.text(),
                                                    self.signals.signal_int, self)
        
        if self.parent.protein_data == -1:
            self.empty.emit()
        elif self.parent.protein_data == -2:
            self.nocon.emit()
        elif self.parent.protein_data == -3:
            self.abort.emit()
        else:
            self.blockSignals(False)


class nmpWorker(QThread):
    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_NMPprogress)
        self.unp_result = parent.protein_data
        self.allele = parent.allele
        self.len = parent.len
        self.thd = parent.thd
        self.dir = parent.nmpdir
        self.parent = parent
        self.starting = parent.nstart
        self.empty = parent.nempty
        self.nocon = parent.nnocon
        self.abort = parent.nabort

    def run(self):
        self.starting.emit()
        self.blockSignals(True)
        self.parent.prediction_data = Nmp.exec_netmhcpan(self.parent.protein_data, [self.allele.currentText(),
                                                                                    self.len.currentText(),
                                                                                    self.thd.text()], self.dir.text(),
                                                         self.signals.signal_int, self)
        if self.parent.prediction_data == -1:
            self.empty.emit()
        elif self.parent.prediction_data == -2:
            self.nocon.emit()
        elif self.parent.prediction_data == -3:
            self.abort.emit()
        else:
            self.blockSignals(False)


class hmpWorker(QThread):
    def __init__(self, parent=None):
        QThread.__init__(self,parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_HMPprogress)
        self.unp_result = parent.protein_data
        self.nmp_result = parent.prediction_data
        self.allele = parent.allele
        self.top = parent.top
        self.hco = parent.hco
        self.dir = parent.hmpdir
        self.parent = parent
        self.starting = parent.hstart
        self.abort = parent.habort
        self.empty = parent.hempty

    def run(self):

        self.starting.emit()
        self.blockSignals(True)

        self.parent.heatmap_data = Hmp.produce_heatmap(self.parent.protein_data, self.parent.prediction_data,
                                                       self.dir.text(), self.parent.top.text(), self.parent.hco.text(),
                                                       self.signals.signal_int, self)

        if self.parent.heatmap_data == -1:
            self.empty.emit()
        elif self.parent.heatmap_data == -2:
            self.abort.emit()
        else:
            self.blockSignals(False)

def open_win():
    app = QApplication()

    window = MainWindow4()
    window.show()

    sys.exit(app.exec())


#open_win()
