from qtpy.QtWidgets import (
    QLabel,
    QWidget,
    QGridLayout,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSlider,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from cad_ui.general import CADMainWindow
from cad_ui.plotting import CadPlot
import numpy as np
import blt
import os
import math


class Utilities:
    def bltvector(self, n):
        try:
            # assumes n is a global variable
            blt.vector(n)
        except:
            pass

    def environment(self):
        # I'm iffy about this, there's a different way to find the path
        env = os.environ
        try:
            path = env["PATH"]
        except KeyError:
            return 1
        paths = path.split(":")


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


class SecondRF(QWidget):
    def __init__(self) -> None:
        super().__init__()
        # translating global variables to class variables
        self.n = 2
        self.Vrf = 0
        self.Vn = 0
        self.vrf = 0
        self.vn = 0
        self.theta = 0
        self.phis_1 = 0
        self.vzerox = 0
        self.bdot = 0
        self.vzeroy = 0
        self.g2rfsep = 0
        self.g2rflist = 0
        self.v1max = 0
        self.e = 0
        self.h = 0
        self.A_1 = 0
        self.B_1 = 0
        self.B1 = 0
        self.sop = 0
        self.A2bun = 0
        self.fnu2rf = 0
        self.radeg = 180 / math.pi

        self.setWindowTitle("Dr. BBat")
        quit = QPushButton("Quit")
        redraw = QPushButton("Redraw")
        vv = QPushButton("V")
        printButton = QPushButton("Print")
        help = QPushButton("Help")

        layout = QGridLayout()
        layout.addWidget(quit, 0, 0, 1, 1)
        layout.addWidget(redraw, 0, 1, 1, 1)
        layout.addWidget(printButton, 0, 3, 1, 1)
        layout.addWidget(vv, 0, 2, 1, 1)
        layout.addWidget(help, 0, 4, 1, 1)

        title = QLabel("Editing RF Parameters")
        font = QFont()
        font.setBold(True)
        title.setFont(font)
        layout.addWidget(title, 1, 2, 1, 3)

        vrf = QLabel("Vrf")
        vrf_line = QLineEdit("0")
        layout.addWidget(vrf, 2, 0, 1, 2)
        layout.addWidget(vrf_line, 2, 2, 1, 3)

        vn = QLabel("Vn")
        vn_line = QLineEdit("100.0")
        layout.addWidget(vn, 3, 0, 1, 2)
        layout.addWidget(vn_line, 3, 2, 1, 3)

        n = QLabel("n")
        n_line = QLineEdit("2")
        layout.addWidget(n, 4, 0, 1, 2)
        layout.addWidget(n_line, 4, 2, 1, 3)

        t = QLabel("Theta (deg)")
        t_line = QLineEdit("86")
        layout.addWidget(t, 5, 0, 1, 2)
        layout.addWidget(t_line, 5, 2, 1, 3)

        vn_vrf = QLabel("Vn/Vrf")
        vnvrf_line = QLineEdit("1.0")
        layout.addWidget(vn_vrf, 6, 0, 1, 2)
        layout.addWidget(vnvrf_line, 6, 2, 1, 3)

        phis = QLabel("Phis (deg)")
        phis_line = QLineEdit("0")
        layout.addWidget(phis, 7, 0, 1, 2)
        layout.addWidget(phis_line, 7, 2, 1, 3)

        phi2s = QLabel("Phi2s = (phis+theta)(deg)")
        phi2s_line = QLineEdit("86")
        layout.addWidget(phi2s, 8, 0, 1, 2)
        layout.addWidget(phi2s_line, 8, 2, 1, 3)

        phis_1 = QSlider(Qt.Horizontal)
        phis_1.setMinimum(-180)
        phis_1.setMaximum(180)
        phis_1.setValue(0)
        phis_1.setTickPosition(QSlider.TicksBelow)
        phis1_l = QLabel("phis_1")
        layout.addWidget(phis1_l, 9, 0, 1, 1)
        layout.addWidget(phis_1, 9, 1, 1, 4)

        V_1 = QSlider(Qt.Horizontal)
        V_1.setMinimum(0)
        V_1.setMaximum(100)
        V_1.setValue(100)
        V_1.setTickPosition(QSlider.TicksBelow)
        V1_l = QLabel("V_1")
        layout.addWidget(V1_l, 10, 0, 1, 1)
        layout.addWidget(V_1, 10, 1, 1, 4)

        Vn = QSlider(Qt.Horizontal)
        Vn.setMinimum(0)
        Vn.setMaximum(100)
        Vn.setValue(100)
        Vn.setTickPosition(QSlider.TicksBelow)
        vn_l = QLabel("V_n")
        layout.addWidget(vn_l, 11, 0, 1, 1)
        layout.addWidget(Vn, 11, 1, 1, 4)

        theta = QSlider(Qt.Horizontal)
        theta.setMinimum(0)
        theta.setMaximum(100)
        theta.setValue(86)
        theta.setTickPosition(QSlider.TicksBelow)
        theta_l = QLabel("theta")
        layout.addWidget(theta_l, 12, 0, 1, 1)
        layout.addWidget(theta, 12, 1, 1, 4)

        k0 = QSlider(Qt.Horizontal)
        k0.setMinimum(0)
        k0.setMaximum(100)
        k0.setValue(0)
        k0.setTickPosition(QSlider.TicksBelow)
        k0_l = QLabel("k_0")
        layout.addWidget(k0_l, 13, 0, 1, 1)
        layout.addWidget(k0, 13, 1, 1, 4)

        W_pp = QPushButton("++W")
        W_mm = QPushButton("--W")
        W_line = QLineEdit("1")
        U_pp = QPushButton("++U")
        U_mm = QPushButton("--U")
        U_line = QLineEdit("10")
        layout.addWidget(W_pp, 14, 0, 1, 1)
        layout.addWidget(W_mm, 14, 1, 1, 1)
        layout.addWidget(W_line, 14, 2, 1, 1)
        layout.addWidget(U_pp, 15, 0, 1, 1)
        layout.addWidget(U_mm, 15, 1, 1, 1)
        layout.addWidget(U_line, 15, 2, 1, 1)

        refresh = QPushButton("refresh")
        layout.addWidget(refresh, 16, 1, 1, 3)

        srfPlot = CadPlot()
        srfPlot.addDataset(
            "Random", np.arange(0, 10, 1), np.arange(0, 10, 1), color="b"
        )
        layout.addWidget(srfPlot, 0, 5, 20, 4)

        textWid = QWidget()
        textLayout = QGridLayout()

        phase = QLabel("phase(deg):")
        phase_value = QLabel("54.1689284")
        phasens = QLabel("phase(ns):")
        phasens_value = QLabel("2.45909358")
        textLayout.addWidget(phase, 0, 0, 1, 1)
        textLayout.addWidget(phase_value, 0, 1, 1, 2)
        textLayout.addWidget(phasens, 0, 3, 1, 1)
        textLayout.addWidget(phasens_value, 0, 4, 1, 2)
        de_mev = QLabel("dE(MeV):")
        demev_line = QLabel("something")
        de_es = QLabel("dE/Es(10^-3):")
        dees_line = QLabel("something")
        textLayout.addWidget(de_mev, 1, 0, 1, 1)
        textLayout.addWidget(demev_line, 1, 1, 1, 1)
        textLayout.addWidget(de_es, 1, 3, 1, 1)
        textLayout.addWidget(dees_line, 1, 4, 1, 1)

        textWid.setLayout(textLayout)
        layout.addWidget(textWid, 20, 5, 4, 4)

        self.setLayout(layout)


class ButtonPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        quit = QPushButton("Quit")
        refresh = QPushButton("Refresh")
        printButton = QPushButton("Print")
        config = QPushButton("Config")
        self.secondRf = QPushButton("Second RF")
        self.srf = None
        self.secondRf.clicked.connect(self.show_SRF)
        help = QPushButton("Help")
        layout = QGridLayout()
        layout.addWidget(quit, 0, 0, 1, 1)
        layout.addWidget(refresh, 0, 1, 1, 1)
        layout.addWidget(printButton, 0, 2, 1, 1)
        layout.addWidget(config, 0, 3, 1, 1)
        layout.addWidget(self.secondRf, 0, 4, 1, 1)
        layout.addWidget(help, 0, 5, 1, 1)
        self.setLayout(layout)

    def show_SRF(self, checked):
        self.srf = SecondRF()
        self.srf.show()


class Window(CADMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(title="bbat", menubar=False, *args, **kwargs)
        self.w = None
        self.setWindowTitle("bbat")
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
            self.w.close()
            self.w = None
