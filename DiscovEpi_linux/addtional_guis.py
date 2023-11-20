#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscovEpi
# Version 1.0
# A tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 20.11.2023: Product.
#################################################

import time
import platform
import ctypes.wintypes
from PySide6.QtWidgets import *


class DeleteRedundantDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.answer = None
        self.setWindowTitle("Redundant Sequences")

        self.deleteBox = QDialogButtonBox()
        self.deleteBox.addButton("Yes", QDialogButtonBox.AcceptRole)
        self.deleteBox.addButton("No", QDialogButtonBox.RejectRole)
        self.deleteBox.addButton("Back", QDialogButtonBox.HelpRole)

        self.deleteBox.accepted.connect(self.yes_clicked)
        self.deleteBox.rejected.connect(self.no_clicked)
        self.deleteBox.helpRequested.connect(self.back_clicked)

        self.layout = QVBoxLayout()
        self.message = QLabel("Do You Want To Discard Redundant Sequences?")
        self.layout.addWidget(self.message)
        self.layout.addWidget(self.deleteBox)
        self.setLayout(self.layout)

        self.show()

    def yes_clicked(self):
        self.answer = 1
        self.close()
        print("yes")

    def no_clicked(self):
        self.answer = 0
        self.close()
        print("no")

    def back_clicked(self):
        self.answer = -1
        self.close()
        print("back")


def get_folder():
    if platform.system() == "Windows":
        dir_docs = 5  # My Documents
        dir_current = 0  # Get current, not default value

        buf = ctypes.create_unicode_buffer(ctypes.wintypes.MAX_PATH)
        ctypes.windll.shell32.SHGetFolderPathW(None, dir_docs, None, dir_current, buf)

        return buf.value
    else:
        return ""


class wait_dialog(QDialog):
    def __init__(self, parent=None):
        super(wait_dialog, self).__init__(parent)
        self.setWindowTitle("Process ongoing")

        self.message = QLabel("Please wait...")
        self.cnclbtn = QPushButton("Cancel")

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.message)
        self.layout.addWidget(self.cnclbtn)

        self.cnclbtn.clicked.connect(self.close_event(parent))

        self.show()

    def close_event(self, window):
        window.close()

class progress_dialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Progress")
        self.progress_value = None

        self.layout = QVBoxLayout()
        self.layout_horiz1 = QHBoxLayout()
        self.layout_horiz2 = QHBoxLayout()
        self.message = QLabel("Parameter:")
        self.message1 = QLabel("Staphylococcus aureus\nCell Wall\nyes\nyes\nPlatzhalter")

        self.pbar = QProgressBar(self)
        self.pbar.setMinimum(0)
        self.pbar.setMaximum(100)
        self.pbar.setValue(0)

        self.startBtn = QPushButton("Start", self)
        #self.startBtn.clicked.connect(Unp.exec_uniprot())
        self.contBtn = QPushButton("Continue", self)
        #self.contBtn.clicked.connect(TAB WECHSELN / DIALOG SCHLIESSEN)

        self.layout_horiz1.addWidget(self.message)
        self.layout_horiz1.addWidget(self.message1)
        self.layout.addLayout(self.layout_horiz1)
        self.layout_horiz2.addWidget(self.startBtn)
        self.layout_horiz2.addWidget(self.contBtn)
        self.layout.addLayout(self.layout_horiz2)
        self.layout.addWidget(self.pbar)
        self.setLayout(self.layout)


        self.show()

    # Provide estimated time and percentage of completed tasks
    def pbar_update(self, count, count_all, *pred_start_time):

        # self.percentage_old = (count - 1) / count_all * 100

        if pred_start_time:
            print("est time")
            if count == 3:
                end = time.time()
                estimated_seconds = (end - pred_start_time) * count_all / 3
                self.message1.setText("\nEstimated time:" + str(calc_time(estimated_seconds)))

        percentage_new = count / count_all * 100
        print(percentage_new)
        self.pbar.setValue(percentage_new)

        if count == count_all:
            self.contBtn.setEnabled(True)


# Format seconds in hours:minutes:seconds.
def calc_time(passed_seconds):
    seconds = passed_seconds % 60
    minutes = passed_seconds / 60 % 60
    hours = passed_seconds / 3600
    runtime = "%02d:%02d:%02d" % (hours, minutes, seconds)
    return runtime

