#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscoEpi
# A Tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 05.05.2023: Product.
#################################################

import time
import re
from pathlib import Path
import webbrowser
from PySide6.QtWidgets import *
from PySide6.QtCore import QObject, QThread, Signal, Slot
import uniProt_gui as Unp
import NetMHCpan_gui as Nmp
import addtional_guis as add
import heatmap_gui as Hmp


class MainWindow4(QMainWindow):
    def __init__(self):
        super(MainWindow4, self).__init__()
        self.protein_data = None
        self.prediction_data = None
        self.heatmap_data = None
        self.del_redundant = None
        self.progress_value = None
        self.UNPprogressBar = QProgressBar()
        self.NMPprogressBar = QProgressBar()
        self.HMPprogressBar = QProgressBar()

        self.setWindowTitle("Epitope Prediction")

        self.tabs = QTabWidget(self)
        self.tabs.setTabPosition(QTabWidget.North)
        self.tabs.setDocumentMode(True)

        ugrid = QGridLayout(self)
        ugrid.addWidget(QLabel("Organism"), 1, 0)
        ugrid.addWidget(QLabel("Location"), 2, 0)
        ugrid.addWidget(QLabel("Format"), 3, 0)
        ugrid.addWidget(QLabel("Directory"), 9, 0)

        self.org = QLineEdit("SARS-CoV-2")
        self.loc = QLineEdit("Membrane")
        self.fmt = QLineEdit("Tsv")
        self.rev = QCheckBox("Reviewed")
        self.rev.setChecked(True)

        uBtnLayout = QHBoxLayout()
        self.nmpBtn = QPushButton("Continue to NetMHCpan")
        self.nmpBtn.setEnabled(False)
        self.subUP = QPushButton("Submit")
        self.udirBtn = QPushButton("Go to Directory")
        self.udirBtn.setEnabled(False)
        uBtnLayout.addWidget(self.subUP)
        uBtnLayout.addWidget(self.nmpBtn)
        uBtnLayout.addWidget(self.udirBtn)

        verticalSpacer = QSpacerItem(300, 20, QSizePolicy.Expanding)
        self.dir = QLineEdit(add.get_folder()+"\Epitope_Prediction\\")
        ugrid.addWidget(self.org, 1, 1)
        ugrid.addWidget(self.loc, 2, 1)
        ugrid.addWidget(self.fmt, 3, 1)
        ugrid.addWidget(self.rev, 4, 1)
        ugrid.addLayout(uBtnLayout, 6, 1)
        ugrid.addWidget(self.UNPprogressBar, 7, 1)
        ugrid.addItem(verticalSpacer, 8, 1)
        ugrid.addWidget(self.dir, 9, 1)

        ngrid = QGridLayout(self)
        ngrid.addWidget(QLabel("Allele"), 1, 0)
        ngrid.addWidget(QLabel("Epitope"), 2, 0)
        ngrid.addWidget(QLabel("Threshold"), 3, 0)

        self.allele = QLineEdit("HLA-A*02:01")
        self.len = QLineEdit("9")
        self.thd = QLineEdit("3.0")

        nBtnLayout = QHBoxLayout()
        self.ndirBtn = QPushButton("Go to Directory")
        self.ndirBtn.setEnabled(False)
        self.subNMP = QPushButton("Submit")
        self.subNMP.setEnabled(False)
        self.hmpBtn = QPushButton("Continue to Heatmaps")
        self.hmpBtn.setEnabled(False)
        nBtnLayout.addWidget(self.subNMP)
        nBtnLayout.addWidget(self.hmpBtn)
        nBtnLayout.addWidget(self.ndirBtn)

        ngrid.addWidget(self.allele, 1, 1)
        ngrid.addWidget(self.len, 2, 1)
        ngrid.addWidget(self.thd, 3, 1)
        ngrid.addLayout(nBtnLayout, 4, 1)
        ngrid.addWidget(self.NMPprogressBar, 5, 1)

        hgrid = QGridLayout(self)
        hgrid.addWidget(QLabel("Heatmap Length Cutoff"), 1, 0)
        hgrid.addWidget(QLabel("Without cutoff every protein will be displayed in full length."), 2, 0, 1, 2)
        self.hco = QLineEdit("")

        hBtnLayout = QHBoxLayout()
        self.hdirBtn = QPushButton("Go to Directory")
        self.hdirBtn.setEnabled(False)
        self.subHMP = QPushButton("Submit")
        self.subHMP.setEnabled(False)
        hBtnLayout.addWidget(self.subHMP)
        hBtnLayout.addWidget(self.hdirBtn)

        hgrid.addWidget(self.hco, 1, 1)
        hgrid.addLayout(hBtnLayout, 3, 0, 1, 2)
        hgrid.addWidget(self.HMPprogressBar, 4, 0, 1, 2)


        uwidget = QWidget(self)
        uwidget.setLayout(ugrid)
        nwidget = QWidget(self)
        nwidget.setLayout(ngrid)
        hwidget = QWidget(self)
        hwidget.setLayout(hgrid)
        self.tabs.addTab(uwidget, "UniProt")
        self.tabs.addTab(nwidget, "NetMHCpan")
        self.tabs.addTab(hwidget, "Heatmap")

        self.setCentralWidget(self.tabs)

        self.nmpBtn.clicked.connect(self.on_nmpBtn_clicked)
        self.hmpBtn.clicked.connect(self.on_hmpBtn_clicked)
        self.subUP.clicked.connect(self.on_submitUP_clicked)
        self.subNMP.clicked.connect(self.on_submitNMP_clicked)
        self.subHMP.clicked.connect(self.on_submitHMP_clicked)
        self.udirBtn.clicked.connect(self.on_directoryButton_clicked)
        self.ndirBtn.clicked.connect(self.on_directoryButton_clicked)
        self.hdirBtn.clicked.connect(self.on_directoryButton_clicked)

        self.unp = unpWorker(self)
        self.nmp = nmpWorker(self)
        self.hmp = hmpWorker(self)
        self.unp.finished.connect(self.unp_finished)
        self.nmp.finished.connect(self.nmp_finished)
        self.hmp.finished.connect(self.hmp_finished)

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

        Path(self.dir.text()).mkdir(parents=True, exist_ok=True)

        start = time.time()
        self.unp.start()
        end = time.time()
        print("UniProt done in " + add.calc_time(end - start))

    def on_submitNMP_clicked(self):
        start = time.time()
        self.nmp.start()
        end = time.time()
        print("NetMHCpan done in " + add.calc_time(end-start))

    def on_directoryButton_clicked(self):
        webbrowser.open('file:///' + self.dir.text())

    def on_submitHMP_clicked(self):
        start = time.time()
        self.hmp.start()
        end = time.time()
        print("Heatmap produced in " + add.calc_time(end-start))

    def unp_finished(self):
        self.udirBtn.setEnabled(True)
        self.nmpBtn.setEnabled(True)
        self.subNMP.setEnabled(True)
        self.dir.setEnabled(False)
        print("UniProt finished")

    def nmp_finished(self):
        self.hmpBtn.setEnabled(True)
        self.ndirBtn.setEnabled(True)
        self.subHMP.setEnabled(True)
        print("NetMHCpan finished")

    def hmp_finished(self):
        self.hdirBtn.setEnabled(True)
        print("Heatmaps finished")

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
    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_UNPprogress)
        self.result = None
        self.org = parent.org
        self.loc = parent.loc
        self.rev = parent.rev
        self.fmt = parent.fmt
        self.del_redundant = parent.del_redundant
        self.dir = parent.dir
        self.parent = parent

    def run(self):
        self.parent.protein_data = Unp.exec_uniprot([re.sub(r"[\t\n]", "", self.org.text()), self.loc.text(),
                                                     self.rev.isChecked()],
                                                    self.fmt.text(), self.del_redundant, self.dir.text(),
                                                    self.signals.signal_int, self.parent)
        print(self.parent.protein_data)
        if self.parent.protein_data == -1:
            #e = QErrorMessage(self.parent.MainWindow4)
            #e.showMessage("UniProt error occured.\n Check batch output and try again.")
            self.signals.signal_int.emit(0)
        else:
            self.signals.signal_int.emit(100)


class nmpWorker(QThread):
    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_NMPprogress)
        self.unp_result = parent.protein_data
        self.allele = parent.allele
        self.len = parent.len
        self.thd = parent.thd
        self.dir = parent.dir
        self.parent = parent

    def run(self):
        if self.parent.protein_data == -1:
            #e = QErrorMessage(self.parent.MainWindow4)
            #e.showMessage("NetMHCpan error occured.\n Check batch output and try again.")
            self.signals.signal_int.emit(0)
        else:
            self.parent.prediction_data = Nmp.exec_netmhcpan(self.parent.protein_data, [self.allele.text(), self.len.text(),
                                                                               self.thd.text()],
                                                             self.dir.text(), self.signals.signal_int)
            self.signals.signal_int.emit(100)


class hmpWorker(QThread):
    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.signals = MySignals()
        self.signals.signal_int.connect(parent.update_HMPprogress)
        self.unp_result = parent.protein_data
        self.nmp_result = parent.prediction_data
        self.allele = parent.allele
        self.hco = parent.hco
        self.dir = parent.dir
        self.parent = parent

    def run(self):
        self.signals.signal_int.emit(1)
        self.parent.heatmap_data = Hmp.produce_heatmap(self.parent.protein_data, self.parent.prediction_data,
                                                       self.dir.text(), self.parent.hco.text(), self.signals.signal_int)
        self.signals.signal_int.emit(100)

def open_win():
    app = QApplication()

    window = MainWindow4()
    window.show()

    app.exec()


open_win()
