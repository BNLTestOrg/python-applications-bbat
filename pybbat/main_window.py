# Please use relative imports!
# https://www.python.org/dev/peps/pep-0328/#rationale-for-relative-imports
from qtpy.QtWidgets import (
    QLabel,
    QWidget,
    QGridLayout,
    QLineEdit,
    QPushButton,
    QComboBox,
)
from qtpy.QtGui import QFont
from cad_ui.general import CADMainWindow
from cad_ui.plotting import CadPlot
import numpy as np


class MoreResults(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("bkt")
        layout = QGridLayout()

        title = QLabel(
            "Update these variables by clicking on either Bucket OK or Bunch OK"
        )
        layout.addWidget(title, 0, 0, 1, 2)

        bfield = QLabel("Bfield (T)")
        bline = QLineEdit("0.30662093")
        layout.addWidget(bfield, 1, 0, 1, 1)
        layout.addWidget(bline, 1, 1, 1, 1)

        frf = QLabel("frf (MHz)")
        frfline = QLineEdit("28.0")
        layout.addWidget(frf, 2, 0, 1, 1)
        layout.addWidget(frfline, 2, 1, 1, 1)

        gamma = QLabel("gamma")
        gline = QLineEdit("9.68")
        layout.addWidget(gamma, 3, 0, 1, 1)
        layout.addWidget(gline, 3, 1, 1, 1)

        beta = QLabel("beta")
        bline = QLineEdit("0.994649")
        layout.addWidget(beta, 4, 0, 1, 1)
        layout.addWidget(bline, 4, 1, 1, 1)

        eta = QLabel("eta")
        eline = QLineEdit("-0.0087482")
        layout.addWidget(eta, 5, 0, 1, 1)
        layout.addWidget(eline, 5, 1, 1, 1)

        kineng = QLabel("Kinetic energy (GeV/u)")
        keline = QLineEdit("8.083413")
        layout.addWidget(kineng, 6, 0, 1, 1)
        layout.addWidget(keline, 6, 1, 1, 1)

        moment = QLabel("Momentum (GeV/c/u)")
        momline = QLineEdit("8.966442")
        layout.addWidget(moment, 7, 0, 1, 1)
        layout.addWidget(momline, 7, 1, 1, 1)

        bunlen = QLabel("Bunch length (deg)")
        bunLine = QLineEdit("20.16")
        layout.addWidget(bunlen, 8, 0, 1, 1)
        layout.addWidget(bunLine, 8, 1, 1, 1)

        buclen = QLabel("Bucket length (deg)")
        bucLine = QLineEdit("360.0")
        layout.addWidget(buclen, 9, 0, 1, 1)
        layout.addWidget(bucLine, 9, 1, 1, 1)

        nslen = QLabel("Bucket length (ns)")
        nsLine = QLineEdit("35.714285")
        layout.addWidget(nslen, 10, 0, 1, 1)
        layout.addWidget(nsLine, 10, 1, 1, 1)

        self.setLayout(layout)


class ControlPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        layout = QGridLayout()

        headerFont = QFont()
        headerFont.setBold(True)

        mb_header = QLabel("Editing machine and beam parameters")
        mb_header.setFont(headerFont)
        layout.addWidget(mb_header, 0, 0, 1, 3)

        machine = QLabel("Machine")
        layout.addWidget(machine, 1, 0, 1, 2)
        machine_line = QComboBox()
        machine_line.addItems(["AGS", "Booster", "RHIC", "others"])
        layout.addWidget(machine_line, 1, 2, 1, 1)

        bd_header = QLabel("B Dot (T/s)")
        bd_header.setFont(headerFont)
        layout.addWidget(bd_header, 2, 0, 1, 2)
        bd_line = QLineEdit("0.0")
        layout.addWidget(bd_line, 2, 2, 1, 1)

        rff_header = QComboBox()
        rff_header.addItems(
            [
                "RF Frequency (MHz)",
                "B field (T)",
                "RF frequency (MHz)",
                "Gamma",
                "pc (GeV/c/u)",
                "K (MeV/u)",
                "K (GeV/u)",
            ]
        )
        layout.addWidget(rff_header, 3, 0, 1, 2)
        rff_line = QLineEdit("28.0")
        layout.addWidget(rff_line, 3, 2, 1, 1)

        rfv_header = QLabel("RF voltage/turn (kV)")
        layout.addWidget(rfv_header, 4, 0, 1, 2)
        rfv_line = QLineEdit("100")
        layout.addWidget(rfv_line, 4, 2, 1, 1)

        rfh_header = QLabel("RF harmonic number")
        layout.addWidget(rfh_header, 5, 0, 1, 2)
        rfh_line = QLineEdit("360")
        layout.addWidget(rfh_line, 5, 2, 1, 1)

        w_minus = QPushButton("--W")
        bucket_ok = QPushButton("Bucket OK")
        w_plus = QPushButton("++W")
        layout.addWidget(w_minus, 6, 0, 1, 1)
        layout.addWidget(bucket_ok, 6, 1, 1, 1)
        layout.addWidget(w_plus, 6, 2, 1, 1)

        species_label = QLabel("Species")
        layout.addWidget(species_label, 7, 0, 1, 2)
        species_line = QComboBox()
        species_line.addItems(
            [
                "Proton",
                "Deuteron",
                "Gold",
                "Iron",
                "Carbon",
                "Sulfur",
                "Silicon",
                "Copper",
                "Iodine",
                "Uranium",
                "Ruthenium",
                "Zirconium",
                "others",
            ]
        )
        layout.addWidget(species_line, 7, 2, 1, 1)

        charge_header = QLabel("Charge State")
        charge_header.setFont(headerFont)
        layout.addWidget(charge_header, 8, 0, 1, 2)
        charge_line = QLineEdit("79")
        layout.addWidget(charge_line, 8, 2, 1, 1)

        # bunch_header = QLabel("Bunch Length (ns)")
        bunch_header = QComboBox()
        bunch_header.addItems(
            [
                "Bunch Length (ns)",
                "Bunch Emittance (eVs/u)",
                "Bunch Length (rad)",
                "Bunch Length (deg)",
            ]
        )
        layout.addWidget(bunch_header, 9, 0, 1, 2)
        bunch_line = QLineEdit("2")
        layout.addWidget(bunch_line, 9, 2, 1, 1)

        bw_minus = QPushButton("--W")
        bunch_ok = QPushButton("Bunch OK")
        bw_plus = QPushButton("++W")
        layout.addWidget(bw_minus, 10, 0, 1, 1)
        layout.addWidget(bunch_ok, 10, 1, 1, 1)
        layout.addWidget(bw_plus, 10, 2, 1, 1)

        synchPhase = QLabel("Synchronous Phase (deg)")
        synchPhase.setFont(headerFont)
        layout.addWidget(synchPhase, 11, 0, 1, 2)
        synchLine = QLineEdit("0.0")
        layout.addWidget(synchLine, 11, 2, 1, 1)

        statBuck = QLabel("Stationary Bucket Area (eVs/u)")
        statBuck.setFont(headerFont)
        layout.addWidget(statBuck, 12, 0, 1, 2)
        statLine = QLineEdit("0.3866")
        layout.addWidget(statLine, 12, 2, 1, 1)

        movBuck = QLabel("Moving Bucket Area (eVs/u)")
        movBuck.setFont(headerFont)
        layout.addWidget(movBuck, 13, 0, 1, 2)
        movLine = QLineEdit("0.3866")
        layout.addWidget(movLine, 13, 2, 1, 1)

        synchFreq = QLabel("Synchronous Frequency (Hz)")
        synchFreq.setFont(headerFont)
        layout.addWidget(synchFreq, 14, 0, 1, 2)
        synchFreq = QLineEdit("116.76")
        layout.addWidget(synchFreq, 14, 2, 1, 1)

        bunchArea = QLabel("Bunch Area (eVs/u)")
        bunchArea.setFont(headerFont)
        layout.addWidget(bunchArea, 15, 0, 1, 2)
        bunchLine = QLineEdit("0.0023")
        layout.addWidget(bunchLine, 15, 2, 1, 1)

        bunchLength = QLabel("Bunch Length (ns)")
        bunchLength.setFont(headerFont)
        layout.addWidget(bunchLength, 16, 0, 1, 2)
        LengthLine = QLineEdit("2")
        layout.addWidget(LengthLine, 16, 2, 1, 1)

        self.setLayout(layout)


class ButtonPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        quit = QPushButton("Quit")
        refresh = QPushButton("Refresh")
        printButton = QPushButton("Print")
        config = QPushButton("Config")
        secondRf = QPushButton("Second RF")
        help = QPushButton("Help")

        layout = QGridLayout()
        layout.addWidget(quit, 0, 0, 1, 1)
        layout.addWidget(refresh, 0, 1, 1, 1)
        layout.addWidget(printButton, 0, 2, 1, 1)
        layout.addWidget(config, 0, 3, 1, 1)
        layout.addWidget(secondRf, 0, 4, 1, 1)
        layout.addWidget(help, 0, 5, 1, 1)
        self.setLayout(layout)


class Window(CADMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(title="pybbat", menubar=False, *args, **kwargs)
        self.w = None
        self.setWindowTitle("pybbat")
        control = ControlPanel()
        buttons = ButtonPanel()
        layout = QGridLayout()
        layout.addWidget(control, 0, 0, 15, 3)
        more = QPushButton("More Results")
        more.clicked.connect(self.show_new_window)
        layout.addWidget(more, 15, 1, 1, 1)
        layout.addWidget(buttons, 0, 3, 1, 1)

        self.plot = CadPlot()
        self.plot.addDataset(
            "Random", np.arange(0, 10, 1), np.arange(0, 10, 1), color="b"
        )
        layout.addWidget(self.plot, 1, 3, 10, 1)

        wid = QWidget()
        wid.setLayout(layout)
        self.setCentralWidget(wid)

    def show_new_window(self, checked):
        if self.w is None:
            self.w = MoreResults()
            self.w.show()
        else:
            self.w.close()  # Close window.
            self.w = None  # Discard reference.
