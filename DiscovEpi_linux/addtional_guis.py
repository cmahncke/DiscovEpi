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

import time
import platform
import ctypes.wintypes
from PySide6.QtWidgets import *
import os


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
        return os.path.expanduser("~/Documents")


# Format seconds in hours:minutes:seconds.
def calc_time(passed_seconds):
    seconds = passed_seconds % 60
    minutes = passed_seconds / 60 % 60
    hours = passed_seconds / 3600
    runtime = "%02d:%02d:%02d" % (hours, minutes, seconds)
    return runtime
