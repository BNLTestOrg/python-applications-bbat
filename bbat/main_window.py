from qtpy.QtWidgets import (
    QLabel,
    QWidget,
    QGridLayout,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSpinBox,
    QTabWidget,
    QMessageBox,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from cad_ui.general import CADMainWindow, PrintMenu
from cad_ui.plotting import CadPlot
import math
import pyqtgraph as pg

from bbat.btools import bTools
from .bmath import bmath


class BBat(CADMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(title="bbat", menubar=False, *args, **kwargs)
        self.machine = "RHIC"
        self.n = 2
        self.Vrf = 100
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vrf_i = 0
        self.Vn = 100
        self.Vn_k = self.Vn * bmath.kilo
        self.theta = 86
        self.phis1 = 10
        self.phis_1 = 0
        self.df_f = 0
        self.vzerox = [-190, 190]
        self.vzeroy = [0, 0]
        self.bdot = (8.7 * 201 * 13.7) / 1000
        self._g2rfsep = []
        self.sep_x = {}
        self.sep_y = {}
        self.g2rflist = []
        self.gammas_1 = 0
        self.betas_1 = 1
        # ****** these are undefined IDK ABOUT THIS ******
        # self._g2rflist = [self.RFV2, self.RFU2, self.vzero, self.vinf, self.phis1, self.hline]
        # v1max is in para.py
        self.fnu = 116.764913687
        self.Tnu = 0
        self.st_bkt = 0.3866
        self.dfnu = 0
        self.e = 79
        self.h = 360
        self.A_1 = 1
        self.B_1 = 1
        self.B1 = self.e / (2 * bmath.pi * self.h)  # B1 is not B_1, B1= Vrf*B_1;
        self.sop = 0
        self.sop_x = 0
        self.A2bkt = 0
        self.A2Bun = 0
        self.fnu2rf = 0
        self.dR_R = 0
        self.W2Max = 3
        self.W2Min = -3
        self.U2Max = 200
        self.U2Min = 200
        self.RFU2_x = []
        self.RFU2_y = []
        self.RFV2_x = []
        self.RFV2_y = []

        self.phis_1_x = []
        self.phis_1_y = []
        self.vinf_x = []
        self.vinf_y = []
        self.vzero_x = []
        self.vzero_y = []
        self.hline_x = []
        self.hline_y = []

        self.dE = 0
        # self.frf_1 = 28.0
        self.frf_1 = 0
        self.dE_Es = 0
        self.Es = 1
        # this should be a list
        self.phis = 0
        # self.fnu1 = 0
        self.dP = 0
        self.dP_Ps = 0
        self.betas1 = 1
        self.lotheta = -180 / self.n
        self.uptheta = 180 / self.n
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.phi2s = self.phis_1 + self.theta
        self.phis_2 = 0
        self.phi12_1 = 0
        self.phi12t_1 = 0
        self.BUNe_1 = 0

        # pi and radeg are in the para
        self.BUNlt = 2
        self.BUNlr = 0
        self.BUNe = 0
        self.BUNld = 0
        self.species = "Gold (au)"
        self.fnu2rf = 0
        self.A = 0
        self.ww = 1
        self.uu = 10

        self.w = None
        self.Bf = 0.30662093
        self.frf = 28
        self.gammas = 5
        self.pc = 8.96644213239
        self.km = 8083.41302244
        self.kg = 8.08341302244
        self.Bdot = 0.0
        self.h = 360
        self.Wmax = 2
        self.Wmin = -2
        self.phase_t1 = 0
        self.etas_1 = 1

        self.BKT_x = []
        self.BKT_y = []
        self.BUN_x = []
        self.BUN_y = []

        self.BUN = "lt"
        self.x = 1
        self.srfWid = None

        self.machineValues = bmath.machineParameters[self.machine]

        self.rf_dict = {
            "RF Frequency (MHz)": self.frf,
            "B field (T)": self.Bf,
            "Gamma": self.gammas,
            "pc (GeV/c/u)": self.pc,
            "K (MeV/u)": self.km,
            "K (GeV/u)": self.kg,
        }
        self.gbfpk_name = {
            "RF Frequency (MHz)": "frf",
            "B field (T)": "bf",
            "Gamma": "gamma",
            "pc (GeV/c/u)": "pc",
            "K (MeV/u)": "km",
            "K (GeV/u)": "kg",
        }
        self.charge_dict = {
            "Proton (p)": bmath.Q_p,
            "Deuteron (d)": bmath.Q_d,
            "Gold (au)": bmath.Q_au,
            "Iron (fe)": bmath.Q_fe,
            "Carbon (c)": bmath.Q_c,
            "Sulfur (s)": bmath.Q_s,
            "Silicon (si)": bmath.Q_si,
            "Copper (cu)": bmath.Q_cu,
            "Iodine (i)": bmath.Q_i,
            "Uranium (u)": bmath.Q_u,
        }

        self.v_wid = QWidget()
        self.v_wid.setWindowTitle("Vrf form")
        self.v_layout = QGridLayout()
        self.v_plot = CadPlot()
        self.vplot_vb = self.v_plot.getViewBox()
        self.vplot_vb.setXRange(-370, 370)
        # set the tick marks at 45 if possible
        self.v_layout.addWidget(self.v_plot, 0, 0, 6, 10)
        self.v_ok = QPushButton("OK")
        self.v_ok.clicked.connect(self.v_wid.close)
        self.v_layout.addWidget(self.v_ok, 6, 0, 1, 10)
        self.v_wid.setLayout(self.v_layout)

        self.vzerox = [-370, 370]
        self.vzeroy = [0, 0]

        # UI components
        self.setWindowTitle("bbat")
        control = self.ControlPanel()
        buttons = self.ButtonPanel()
        layout = QGridLayout()
        layout.addWidget(control, 0, 0, 17, 3)
        more = QPushButton("More Results")
        more.clicked.connect(self.showMoreRes)
        layout.addWidget(more, 17, 1, 1, 1)
        layout.addWidget(buttons, 0, 3, 1, 1)

        self.plot = CadPlot()
        self.plot.plotItem.setLabel("bottom", "RF Phase (deg)")
        self.plot.plotItem.setLabel("left", "W")
        self.plot.plot_vb.setXRange(-190, 370)
        self.plot.plot_vb.setYRange(self.Wmin, self.Wmax)
        self.plot.xmin_autorange = True
        self.plot.ymin_autorange = True
        self.plot.xmax_autorange = True
        self.plot.ymax_autorange = True
        self.hline = pg.InfiniteLine(angle=0, movable=False)
        self.vline = pg.InfiniteLine(angle=90, movable=False)
        self.plot.plotItem.addItem(self.hline)
        self.plot.plotItem.addItem(self.vline)
        self.plot.scene().sigMouseMoved.connect(self.mainMouseMoved)
        layout.addWidget(self.plot, 1, 3, 10, 1)
        textbox = self.textBox()
        layout.addWidget(textbox, 11, 3, 8, 1)

        wid = QWidget()
        wid.setLayout(layout)
        self.setCentralWidget(wid)

        self.bfield_label = QLabel("B field (T)")
        self.bline = QLineEdit("0.306")
        self.frf_label = QLabel("frf (MHz)")
        self.frfline = QLineEdit("28.0")
        self.gamma_label = QLabel("gamma")
        self.gline = QLineEdit("9.68")
        self.beta_label = QLabel("beta")
        self.betaline = QLineEdit("0.994")
        self.eta_label = QLabel("eta")
        self.eline = QLineEdit("-0.008")
        self.kineng_label = QLabel("Kinetic energy (GeV/u)")
        self.keline = QLineEdit("8.083")
        self.moment_label = QLabel("Momentum (GeV/c/u)")
        self.momline = QLineEdit("8.966")
        self.bunlen_label = QLabel("Bunch length (deg)")
        self.bunLine = QLineEdit("20.16")
        self.buclen_label = QLabel("Bucket length (deg)")
        self.bucLine = QLineEdit("360.0")
        self.nslen_label = QLabel("Bucket length (ns)")
        self.nsLine = QLineEdit("35.714")

        self.morRes = QWidget()
        self.MoreResults()

        self.calcBucket()
        self.CBunch()

        self.second_g2rf = CadPlot(show_legend=False)
        self.second_g2rf.plotItem.vb.setXRange(-190, 370)
        self.second_g2rf.plotItem.vb.setYRange(-3, 3)
        self.second_g2rf.xmin_autorange = True
        self.second_g2rf.ymin_autorange = True
        self.second_g2rf.xmax_autorange = True
        self.second_g2rf.ymax_autorange = True
        self.second_g2rf.plotItem.setLabel("bottom", "RF Phase (deg)")
        self.second_g2rf.plotItem.setLabel("left", "W")

        self.g2rf_hline = pg.InfiniteLine(angle=0, movable=False)
        self.g2rf_vline = pg.InfiniteLine(angle=90, movable=False)
        self.second_g2rf.plotItem.addItem(self.g2rf_hline)
        self.second_g2rf.plotItem.addItem(self.g2rf_vline)
        self.second_g2rf.scene().sigMouseMoved.connect(self.mouseMoved)

        for i in range(0, 100):
            self.sep_x["sep_" + str(i)] = []
            self.sep_y["sep_" + str(i)] = []

    def mainMouseMoved(self, event):
        x = event.x()
        y = event.y()
        mousePt = self.plot.plotItem.vb.mapSceneToView(event)
        self.vline.setPos(mousePt.x())
        self.hline.setPos(mousePt.y())
        self.bltCoor(mousePt.x(), mousePt.y())
        self.bun(mousePt.x(), mousePt.y())

    def closeEvent(self, event):
        super().closeEvent(event)
        if self.morRes != None and self.morRes.isVisible():
            self.morRes.close()
        if self.srfWid != None and self.srfWid.isVisible():
            self.srfWid.close()
        if self.v_wid != None and self.v_wid.isVisible():
            self.v_wid.close()

    def setSTbuck(self):
        self.st_bkt = float(self.statLine.text())

    def secondRFWindow(self):
        self.Vrf_V = 0
        self.B1 = self.e / (2 * bmath.pi * self.h)
        self.n = 2
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.bdot = (8.7 * 201 * 13.7) / 1000
        self.srfWid = QWidget()
        self.srfWid.setWindowTitle("Dr. BBat")
        quit = QPushButton("Close")
        quit.clicked.connect(self.srfWid.close)
        redraw = QPushButton("Redraw")
        redraw.clicked.connect(self.draw2rfButtonCommand)
        vv = QPushButton("V")
        vv.clicked.connect(self.vplot)
        printButton = QPushButton("Print")
        pm = PrintMenu(self)
        printButton.setMenu(pm)
        help = QPushButton("Help")
        help.clicked.connect(self.rf2Help)

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
        self.vrf_line = QLineEdit("0")
        self.vrf_line.returnPressed.connect(self.setVrfValue)
        layout.addWidget(vrf, 2, 0, 1, 2)
        layout.addWidget(self.vrf_line, 2, 2, 1, 3)

        vn = QLabel("Vn")
        self.vn_line = QLineEdit("100.0")
        self.vn_line.returnPressed.connect(self.setVnalue)
        layout.addWidget(vn, 3, 0, 1, 2)
        layout.addWidget(self.vn_line, 3, 2, 1, 3)

        n = QLabel("n")
        self.n_line = QLineEdit("2")
        self.n_line.returnPressed.connect(self.set_n)
        layout.addWidget(n, 4, 0, 1, 2)
        layout.addWidget(self.n_line, 4, 2, 1, 3)

        t = QLabel("Theta (deg)")
        self.t_line = QLineEdit("86")
        self.t_line.returnPressed.connect(self.setTValue)
        layout.addWidget(t, 5, 0, 1, 2)
        layout.addWidget(self.t_line, 5, 2, 1, 3)

        vn_vrf = QLabel("Vn/Vrf")
        self.vnvrf_line = QLineEdit("1.0")
        layout.addWidget(vn_vrf, 6, 0, 1, 2)
        layout.addWidget(self.vnvrf_line, 6, 2, 1, 3)

        phis = QLabel("Phis (deg)")
        self.phis_line = QLineEdit("0")
        self.phis_line.returnPressed.connect(self.setValue)
        layout.addWidget(phis, 7, 0, 1, 2)
        layout.addWidget(self.phis_line, 7, 2, 1, 3)

        phi2s = QLabel("Phi2s = (phis+theta)(deg)")
        self.phi2s_line = QLineEdit("86")
        layout.addWidget(phi2s, 8, 0, 1, 2)
        layout.addWidget(self.phi2s_line, 8, 2, 1, 3)

        self.phis_1spinbox = QSpinBox()
        self.phis_1spinbox.setMinimum(-180)
        self.phis_1spinbox.setMaximum(180)
        self.phis_1spinbox.setValue(0)
        self.phis_1spinbox.valueChanged.connect(self.setValueSB)
        phis1_l = QLabel("phis_1")
        layout.addWidget(phis1_l, 9, 0, 1, 1)
        layout.addWidget(self.phis_1spinbox, 9, 1, 1, 4)

        self.V_1spinbox = QSpinBox()
        self.V_1spinbox.setMinimum(0)
        self.V_1spinbox.setMaximum(100)
        self.V_1spinbox.setValue(0)
        self.V_1spinbox.valueChanged.connect(self.setVrfValueSB)
        V1_l = QLabel("V_1")
        layout.addWidget(V1_l, 10, 0, 1, 1)
        layout.addWidget(self.V_1spinbox, 10, 1, 1, 4)

        self.Vn_spinbox = QSpinBox()
        self.Vn_spinbox.valueChanged.connect(self.setVnalueSB)
        self.Vn_spinbox.setMinimum(0)
        self.Vn_spinbox.setMaximum(100)
        self.Vn_spinbox.setValue(100)
        vn_l = QLabel("V_n")
        layout.addWidget(vn_l, 11, 0, 1, 1)
        layout.addWidget(self.Vn_spinbox, 11, 1, 1, 4)

        self.thetaSB = QSpinBox()
        self.thetaSB.valueChanged.connect(self.setTValueSB)
        self.thetaSB.setMinimum(-180)
        self.thetaSB.setMaximum(180)
        self.thetaSB.setValue(86)
        theta_l = QLabel("theta")
        layout.addWidget(theta_l, 12, 0, 1, 1)
        layout.addWidget(self.thetaSB, 12, 1, 1, 4)

        self.k0 = QSpinBox()
        self.k0.valueChanged.connect(self.setKValueSB)
        self.k0.setMinimum(0)
        self.k0.setMaximum(100)
        self.k0.setValue(0)
        k0_l = QLabel("k_0")
        layout.addWidget(k0_l, 13, 0, 1, 1)
        layout.addWidget(self.k0, 13, 1, 1, 4)

        W_pp = QPushButton("++W")
        W_pp.clicked.connect(self.W2_plus)
        W_mm = QPushButton("--W")
        W_mm.clicked.connect(self.W2_minus)
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
        refresh.clicked.connect(self.init_sep)
        refresh.clicked.connect(self.draw2rfButtonCommand)
        layout.addWidget(refresh, 16, 1, 1, 3)

        layout.addWidget(self.second_g2rf, 0, 5, 20, 4)

        textWid = QWidget()
        textLayout = QGridLayout()

        phase = QLabel("phase(deg):")
        self.phase_value2 = QLabel(str(self.sop_x))
        phasens = QLabel("phase(ns):")
        self.phasens_value2 = QLabel(str(self.phase_t1))
        textLayout.addWidget(phase, 0, 0, 1, 1)
        textLayout.addWidget(self.phase_value2, 0, 1, 1, 1)
        textLayout.addWidget(phasens, 0, 2, 1, 1)
        textLayout.addWidget(self.phasens_value2, 0, 3, 1, 1)
        de_mev = QLabel("dE(MeV):")
        self.demev_line2 = QLabel(str(self.dE))
        de_es = QLabel("dE/Es(10^-3):")
        self.dees_line2 = QLabel(str(self.dE_Es))
        textLayout.addWidget(de_mev, 1, 0, 1, 1)
        textLayout.addWidget(self.demev_line2, 1, 1, 1, 1)
        textLayout.addWidget(de_es, 1, 2, 1, 1)
        textLayout.addWidget(self.dees_line2, 1, 3, 1, 1)

        dp_mevc = QLabel("dp(MeV/c):")
        self.dpmevc_val2 = QLabel(str(self.dP))
        dp_ps = QLabel("dP/Ps(10^-3)")
        self.dpps_val2 = QLabel(str(self.dP_Ps))
        textLayout.addWidget(dp_mevc, 2, 0, 1, 1)
        textLayout.addWidget(self.dpmevc_val2, 2, 1, 1, 1)
        textLayout.addWidget(dp_ps, 2, 2, 1, 1)
        textLayout.addWidget(self.dpps_val2, 2, 3, 1, 1)
        self.betas_1 = 1

        dr_r = QLabel("dR/R(10^-3):")
        self.drr_val2 = QLabel(str(self.dR_R))
        df_f = QLabel("df/f(10^3):")
        self.dff_val2 = QLabel(str(self.df_f))
        textLayout.addWidget(dr_r, 3, 0, 1, 1)
        textLayout.addWidget(self.drr_val2, 3, 1, 1, 1)
        textLayout.addWidget(df_f, 3, 2, 1, 1)
        textLayout.addWidget(self.dff_val2, 3, 3, 1, 1)
        self.dR_Ps = 0
        self.betas_1 = 1

        fnu = QLabel("fnu(Hz):")
        # self.fnu_val2 = QLabel(str(self.fnu2rf))
        textLayout.addWidget(fnu, 4, 0, 1, 1)
        # textLayout.addWidget(self.fnu_val2, 4, 1, 1, 1)

        abkt = QLabel("Abkt(eVs/u):")
        # self.abkt_val2 = QLabel(str(self.A2bkt))
        abun = QLabel("Abun(eVs/u):")
        # self.abun_val2 = QLabel(str(self.A2Bun))
        textLayout.addWidget(abkt, 5, 0, 1, 1)
        # textLayout.addWidget(self.abkt_val2, 5, 1, 1, 1)
        textLayout.addWidget(abun, 6, 0, 1, 1)
        # textLayout.addWidget(self.abun_val2, 6, 1, 1, 1)

        textWid.setLayout(textLayout)
        layout.addWidget(textWid, 20, 5, 4, 4)

        self.g2rflist = [
            self.RFV2_x,
            self.RFV2_y,
            self.vzerox,
            self.vzeroy,
            self.vinf_x,
            self.vinf_y,
            self.phis_1_x,
            self.phis_1_y,
            self.hline_x,
            self.hline_y,
        ]
        self.draw2rfButtonCommand()

        self.srfWid.setLayout(layout)
        return self.srfWid

    def mouseMoved(self, event):
        x = event.x()
        y = event.y()
        mousePoint = self.second_g2rf.plotItem.vb.mapSceneToView(event)
        self.g2rf_vline.setPos((mousePoint.x(), 0))
        self.g2rf_hline.setPos((0, mousePoint.y()))
        self.blt2rfCoor(mousePoint.x(), mousePoint.y())
        self.blt2RfBun(mousePoint.x(), mousePoint.y())

    def rf2Help(self):
        self.srf_wid = QWidget()
        srf_layout = QGridLayout()
        title = QLabel("Dr. BBat is a Double RF Bunch and Bucket Analysis Tool")
        titleFont = QFont("Helvetica", 11, QFont.Bold)
        title.setFont(titleFont)
        srf_layout.addWidget(title)

        menu_curves = QLabel(
            "Menu Curves\nGreen\nIt's the total RF waveform V(phi) = vrf*sin(phi)+Vn*Sin(n*phi+n*theta)\n\nBlue\nIt's the RF potential U(phi)=integral{V(phi)}\n\n"
        )
        srf_layout.addWidget(menu_curves)
        curves = QLabel(
            "Orange\nThe buckets.\n\nScales\n To update the green and blue curve scales press --U and ++U\nYou can change the bucket curves by using the --W and ++W buttons."
        )
        srf_layout.addWidget(curves)
        input_label = QLabel(
            "Vrf\nContains the same value as BBat. Any change made in either window will effect each other. This also works for the phis parameter (in degrees).\n\nVn\nVoltage (kV) on the second Rf system.\n\nn\nThe n-th harmonic.\n\nTheta\nThe phase difference between these two rf voltages (-pi/n < theta < pi/n).\n\nPhis\nSynchronous phase angle defined by Vrf*Sin(phis)+Vn*Sin(n*phis+n*theta)=C*rho*Bdot\nIf Vn+0, it has a simple physical meaning: Vrf*sin(phis) is the energy gain per turn. When Vn !=0, the physical meaning is not so easy... but if phis+theta=k*pi (where k is an integer)\nphis regains its simple physical meaning for the second rf system doesn't supply energy to the beam."
        )
        srf_layout.addWidget(input_label)
        self.srf_wid.setLayout(srf_layout)
        self.srf_wid.show()

    def blt2rfCoor(self, x, y):
        self.phase_value2.setText("{:.3f}".format(x))
        self.phase_t1 = (x / 360) / (self.frf_1 + 1.0e-20) * bmath.kilo
        self.phasens_value2.setText("{:.3f}".format(self.phase_t1))
        self.dE = y * 2.0 * bmath.pi * self.frf_1  # frf_1 in MHz, dE in MeV
        self.demev_line2.setText("{:.3f}".format(self.dE))
        self.dE_Es = self.dE / self.Es * bmath.kilo * bmath.mega  # Es in eV
        self.dees_line2.setText("{:.3f}".format(self.dE_Es))
        self.dP = self.dE / self.betas1
        self.dP_Ps = self.betas1 * self.betas1 * self.dE_Es
        self.dpmevc_val2.setText("{:.3f}".format(self.dP))
        self.dpps_val2.setText("{:.3f}".format(self.dP_Ps))
        self.gamma_tr = self.machineValues["gamma_tr"]
        self.dR_R = self.dP_Ps / self.gamma_tr**2
        self.df_f = -self.dP_Ps * self.etas_1
        self.drr_val2.setText("{:.3f}".format(self.dR_R))
        self.dff_val2.setText("{:.3f}".format(self.df_f))

    def set_n(self):
        self.n = float(self.n_line.text())

    def init_sep(self):
        i = int(1 + math.ceil(self.n))
        data_names = list(self.second_g2rf._entries.keys())
        if "2RFU" in data_names:
            for j in self._g2rfsep:
                if j in data_names:
                    self.second_g2rf.removeDataset(j)
                    data_names.remove(j)
            while i >= 0:
                self._g2rfsep.append("sep" + str(i))
                self.second_g2rf.addOrUpdateDataset(
                    "sep" + str(i),
                    self.sep_x["sep_" + str(i)],
                    self.sep_y["sep_" + str(i)],
                    color="orange",
                    width=0.8,
                )
                i = i - 1

    def refresh_command(self):
        data_names = list(self.second_g2rf._entries.keys())
        if "2RFV" in data_names:
            to_delete = self._g2rfsep + self.g2rflist
            for d in to_delete:
                if d in data_names:
                    self.second_g2rf.removeDataset(d)
        self.second_g2rf.plotItem.vb.setXRange(-370, 370)
        self.second_g2rf.plotItem.vb.setYRange(-3, 3)
        self.second_g2rf.addOrUpdateDataset(
            "2RFU", self.RFU2_x, self.RFU2_y, color="blue", width=0.8
        )
        self.second_g2rf.addOrUpdateDataset(
            "2RFV", self.RFV2_x, self.RFV2_y, color="green", width=0.8
        )
        ##############################phis_1
        self.second_g2rf.addOrUpdateDataset(
            "phis_1", self.phis_1_x, self.phis_1_y, color="brown", width=0.8
        )
        ##############################coordinates
        self.second_g2rf.addOrUpdateDataset(
            "vzero", self.vzero_x, self.vzero_y, color="red", width=0.8
        )
        self.second_g2rf.addOrUpdateDataset(
            "vinf", self.vinf_x, self.vinf_y, color="red", width=0.8
        )
        self.second_g2rf.addOrUpdateDataset(
            "hline", self.hline_x, self.hline_y, color="orange", width=0.8
        )
        self.init_sep()

    def draw2rfButtonCommand(self):
        # """ this is the action the redraw button completes
        #    i believe it resets the values in the rf parameter list
        #    and redraws the main lines in the graph. (idk about hline, sep3,2,1,0)
        # """
        self.init_sep()
        vec_2rfv = bTools.Draw2rf(
            self.Vrf_V, self.Vn, self.n, self.theta, self.phase_t1
        )
        self.RFV2_x = vec_2rfv[::2]
        self.RFV2_y = vec_2rfv[1::2]
        # self.second_g2rf.addOrUpdateDataset("2RFV", self.RFV2_x, self.RFV2_y, color="green", width=0.8)

        vec_2rfu = bTools.Draw2rfU(
            self.A_1,
            self.B_1,
            self.Vrf,
            self.Vn,
            self.n,
            self.theta,
            self.phis_1,
            self.phase_t1,
        )
        self.RFU2_x = vec_2rfu[::2]
        self.RFU2_y = vec_2rfu[1::2]
        # self.second_g2rf.addOrUpdateDataset("2RFU", self.RFU2_x, self.RFU2_y, color="blue", width=0.8)

        vec_phis = bTools.DrawPhis(self.phis_1, 0.5)
        vec_phisx = vec_phis[::2]
        vec_phisy = vec_phis[1::2]
        self.second_g2rf.addOrUpdateDataset(
            "phis_1", vec_phisx, vec_phisy, color="brown", width=0.8
        )

        self.second_g2rf.addOrUpdateDataset(
            "vzero", self.vzero_x, self.vzero_y, color="red", width=0.8
        )
        self.second_g2rf.addOrUpdateDataset(
            "vinf", self.vinf_x, self.vinf_y, color="red", width=0.8
        )

        sep_data = bTools.Draw2rfSep(
            self.n,
            self.phase_t1,
            self.theta,
            self.phis_1,
            self.A_1,
            self.Vrf_k,
            self.Vn_k,
            self.B1,
        )
        xsep = sep_data[::2]
        ysep = sep_data[1::2]
        self.second_g2rf.addOrUpdateDataset(
            "sep", xsep, ysep, color="orange", width=0.8
        )
        # print(sep_data)

        # for i in range(0, len(sep_data)):
        #    data = sep_data[i]
        #    self.second_g2rf.addOrUpdateDataset("sep" + str(i), data[0], data[1])
        # Draw2rfSep  ".second.g2rf" "sep" $A_1 $B1 $vrf $vn $n $theta $phis_1;
        # self.A2bkt = bTools.BKT2rf(self.phase_t1, self.n, self.theta, self.phis, self.A_1, self.Vrf, self.Vn, self.B1, self.phis_1)
        contour_data = bTools.Draw2rfHcontour(
            self.theta,
            self.phis_1,
            self.phi2,
            self.n,
            self.phase_t1,
            self.A_1,
            self.Vrf_V,
            self.Vn,
            self.B1,
        )
        self.hline_x = contour_data[::2]
        self.hline_y = contour_data[1::2]
        self.second_g2rf.addOrUpdateDataset(
            "hline", self.hline_x, self.hline_y, color="orange", width=0.8
        )
        self.updateVPlot()

    def helpText(self):
        """This help text will replace the help menu in the widget
        HLP_load .second.menu.help.m "$DIR/cSecondHelp.txt"
        """
        txt = "Dr. BBat is a Double RF Bunch and Bucket Analysis Tool\nCurves\n Green: Total RF Waveform\n"

    # these are for the line edits
    def setVrfValue(self):
        """This function is the same as Vnvalue (line 250)

        Args:
            value (_type_): _description_
        """
        value = float(self.vrf_line.text())
        self.V_1spinbox.setValue(value)
        self.A = value
        self.Vrf = value
        self.Vrf_V = value
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vn_k = self.Vn * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1.0e-20)
        self.lotheta = -180 / self.n
        self.uptheta = 180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.init_sep()
        self.draw2rfButtonCommand()

    def setVnalue(self):
        value = float(self.vn_line.text())
        self.Vn_spinbox.setValue(value)
        self.A = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.lotheta = -180 / self.n
        self.uptheta = 180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.init_sep()
        self.draw2rfButtonCommand()

    def setValue(self):
        value = float(self.phis_line.text())
        self.A = value
        self.phis_1 = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.lotheta = -180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.init_sep()
        self.draw2rfButtonCommand()

    def setTValue(self):
        value = float(self.t_line.text())
        self.thetaSB.setValue(value)
        self.A = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.init_sep()
        self.draw2rfButtonCommand()

    # these are for the spinbox
    def setVrfValueSB(self):
        """This function is the same as Vnvalue (line 250)

        Args:
            value (_type_): _description_
        """
        value = self.V_1spinbox.value()
        self.vrf_line.setText(str(value))
        self.A = value
        self.Vrf = value
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vn_k = self.Vn * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1.0e-20)
        self.lotheta = -180 / self.n
        self.uptheta = 180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.vrf_line.setText(str(self.Vrf))
        self.vnvrf_line.setText(str(self.deltav))
        self.init_sep()
        self.draw2rfButtonCommand()

    def setVnalueSB(self):
        value = self.Vn_spinbox.value()
        self.vn_line.setText(str(value))
        self.Vn = value
        self.A = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.lotheta = -180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.vnvrf_line.setText(str(self.deltav))
        self.init_sep()
        self.draw2rfButtonCommand()

    def setValueSB(self):
        value = self.phis_1spinbox.value()
        self.phis_line.setText(str(value))
        self.phis_1 = value
        self.A = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf + 1e-20)
        self.lotheta = -180 / self.n
        self.phi2s = self.phis_1 + self.theta
        self.phis_line.setText(str(self.phis_1))
        self.phi2s_line.setText(str(self.phi2s))
        self.init_sep()
        self.draw2rfButtonCommand()

    def setKValueSB(self):
        value = self.k0.value()
        self.A = value
        self.k = value
        self.Vrf = self.A * bmath.v1max / 100
        self.vrf_line.setText(str(self.Vrf))
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vn = (100.0 - self.A) * bmath.v1max / 100
        self.vn_line.setText(str(self.Vn))
        self.Vn_k = self.Vn * bmath.kilo
        self.init_sep()
        self.draw2rfButtonCommand()

    def setTValueSB(self):
        value = self.thetaSB.value()
        self.A = value
        self.Vn_k = self.Vn * bmath.kilo
        self.Vrf_k = self.Vrf * bmath.kilo
        self.theta = value
        self.t_line.setText(str(value))
        self.init_sep()
        self.draw2rfButtonCommand()

    def W2_plus(self):
        # connects to W++
        self.W2Max = self.W2Max + self.ww
        self.W2Min = self.W2Min - self.ww
        self.second_g2rf.plotItem.vb.setYRange(self.W2Min, self.W2Max)

    def W2_minus(self):
        # connects to W--
        self.W2Max = self.W2Max - self.ww
        self.W2Min = self.W2Min + self.ww
        self.second_g2rf.plotItem.vb.setYRange(self.W2Min, self.W2Max)

    def U2_plus(self):
        # connects to U++
        self.U2Max = self.U2Max + self.uu
        self.U2Min = self.U2Min - self.uu
        # there is no second y axis
        # configure the y2 axis to have the new max and min values

    def U2_minus(self):
        # connects to U--
        self.U2Max = self.U2Max - self.uu
        self.U2Min = self.U2Min + self.uu
        # there is no second y axis
        # configure the y2 axis to have the new max and min values

    def vplot(self):
        self.updateVPlot()
        self.v_wid.show()

    def updateVPlot(self):
        X_data = []
        Y_data = []
        Yn_data = []
        Yt_data = []
        for i in range(-360, 360):
            X_data.append(i)
            Y_data.append(self.Vrf_V * math.sin(i / bmath.radeg))
            Yn_data.append(self.Vn * math.sin(self.n * (i + self.theta) / bmath.radeg))
            Yt_data.append(
                self.Vrf_V * math.sin(i / bmath.radeg)
                + self.Vn * math.sin(self.n * (i + self.theta) / bmath.radeg)
            )

        self.v_plot.addOrUpdateDataset("Vrf", X_data, Y_data, color="blue", width=0.8)
        self.v_plot.addOrUpdateDataset("Vn", X_data, Yn_data, color="red", width=0.8)
        self.v_plot.addOrUpdateDataset(
            "Vrfn", X_data, Yt_data, color="green", width=0.8
        )
        self.v_plot.addOrUpdateDataset(
            "phis_1",
            [self.phis_1, self.phis1],
            [0, self.Vrf_V],
            color="brown",
            width=0.8,
        )
        self.v_plot.addOrUpdateDataset(
            "vzero", self.vzerox, self.vzeroy, color="red", width=0.8
        )
        self.v_plot.addOrUpdateDataset(
            "vinf", self.vzeroy, self.vzerox, color="red", width=0.8
        )

    def blt2RfBun(self, x, y):
        # scan [$graph invtransform $winX $winY] "%s %s" sop(x) sop(y)
        self.phi2 = x
        contour_data = bTools.Draw2rfHcontour(
            self.theta,
            self.phis_1,
            self.phi2,
            self.n,
            self.phase_t1,
            self.A_1,
            self.Vrf,
            self.Vn,
            self.B1,
        )
        self.hline_x = contour_data[::2]
        self.hline_y = contour_data[1::2]

        self.A2Bun = bTools.BUN2rf(
            self.theta,
            self.phis,
            self.phi2,
            self.n,
            self.phase_t1,
            self.A_1,
            self.Vrf,
            self.Vn,
            self.B_1,
            bmath.A[self.species],
        )
        self.fnu2rf = bTools.fnu2rf(
            self.theta,
            self.phis,
            self.phi2,
            self.n,
            self.phase_t1,
            self.A_1,
            self.Vrf,
            self.Vn,
            self.B_1,
        )

    def MoreResults(self):
        self.morRes.setWindowTitle("bkt")
        layout = QGridLayout()
        title = QLabel(
            "Update these variables by clicking on either Bucket OK or Bunch OK"
        )
        layout.addWidget(title, 0, 0, 1, 2)

        layout.addWidget(self.bfield_label, 1, 0, 1, 1)
        layout.addWidget(self.bline, 1, 1, 1, 1)

        layout.addWidget(self.frf_label, 2, 0, 1, 1)
        layout.addWidget(self.frfline, 2, 1, 1, 1)

        layout.addWidget(self.gamma_label, 3, 0, 1, 1)
        layout.addWidget(self.gline, 3, 1, 1, 1)

        layout.addWidget(self.beta_label, 4, 0, 1, 1)
        layout.addWidget(self.betaline, 4, 1, 1, 1)

        layout.addWidget(self.eta_label, 5, 0, 1, 1)
        layout.addWidget(self.eline, 5, 1, 1, 1)

        layout.addWidget(self.kineng_label, 6, 0, 1, 1)
        layout.addWidget(self.keline, 6, 1, 1, 1)

        layout.addWidget(self.moment_label, 7, 0, 1, 1)
        layout.addWidget(self.momline, 7, 1, 1, 1)

        layout.addWidget(self.bunlen_label, 8, 0, 1, 1)
        layout.addWidget(self.bunLine, 8, 1, 1, 1)

        layout.addWidget(self.buclen_label, 9, 0, 1, 1)
        layout.addWidget(self.bucLine, 9, 1, 1, 1)

        layout.addWidget(self.nslen_label, 10, 0, 1, 1)
        layout.addWidget(self.nsLine, 10, 1, 1, 1)

        ok = QPushButton("OK")
        ok.clicked.connect(self.morRes.close)
        layout.addWidget(ok, 11, 0, 1, 2)

        self.morRes.setLayout(layout)

    def showMoreRes(self):
        self.morRes.show()

    def ButtonPanel(self):
        self.buttonWid = QWidget()
        self.title_font = QFont()
        self.title_font.setPointSize(10.5)
        self.title_font.setBold(True)

        quit = QPushButton("Close")
        quit.clicked.connect(self.close)
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self.refreshPlot)
        printButton = QPushButton("Print")
        pm = PrintMenu(self)
        printButton.setMenu(pm)
        config = QPushButton("Config")
        self.secondRf = QPushButton("Second RF")
        self.srf = None
        self.secondRf.clicked.connect(self.show_SRF)
        help = QPushButton("Help")
        help.clicked.connect(self.helpWidget)
        layout = QGridLayout()
        layout.addWidget(quit, 0, 0, 1, 1)
        layout.addWidget(refresh, 0, 1, 1, 1)
        layout.addWidget(printButton, 0, 2, 1, 1)
        layout.addWidget(config, 0, 3, 1, 1)
        layout.addWidget(self.secondRf, 0, 4, 1, 1)
        layout.addWidget(help, 0, 5, 1, 1)
        self.buttonWid.setLayout(layout)
        return self.buttonWid

    def refreshPlot(self):
        curr_items = list(self.plot._entries.keys())
        if "BKT" in curr_items:
            self.plot.removeDataset("BKT")
        if "BUN" in curr_items:
            self.plot.removeDataset("BUN")
        self.plot.plotItem.vb.setXRange(-190, 370)
        self.BKT_x = []
        self.BKT_y = []
        self.plot.addOrUpdateDataset(
            "BKT", self.BKT_x, self.BKT_y, color="red", width=0.8
        )
        self.BUN_x = []
        self.BUN_y = []
        self.plot.addOrUpdateDataset(
            "BUN", self.BUN_x, self.BUN_y, color="blue", width=0.8
        )

    def show_SRF(self, checked):
        self.srf = self.secondRFWindow()
        self.srf.show()

    def helpWidget(self):
        self.helpWid = QTabWidget()
        self.helpWid.setFixedSize(600, 300)
        self.helpWid.setWindowTitle("Help")
        overview = QWidget()
        overLayout = QGridLayout()
        # title = QLabel("BBAT Version 2 by Jennefer Maldonado (jmaldonad@bnl.gov)")
        # title.setFont(self.title_font)
        body_text = QLabel(
            "BBat pronunced b-bat is a Bunch and Bucket Analysis Tool for a single RF system."
        )
        # overLayout.addWidget(title, 0, 0, 1, 3)
        overLayout.addWidget(body_text, 1, 0, 1, 3)
        overview.setLayout(overLayout)
        self.helpWid.addTab(overview, "Overview")

        machine = QWidget()
        machText = QLabel(
            "Press the machine button to list the built in machines. AGS is the default machine. Choosing others, you'll need to key in some machine parameters. Click ok when you finish typing or chaning some other parameters."
        )
        machText.setWordWrap(True)
        machLO = QGridLayout()
        machLO.addWidget(machText, 0, 0, 1, 1)
        machine.setLayout(machLO)
        self.helpWid.addTab(machine, "Machine")

        gamWid = QWidget()
        gamText = QLabel(
            "This button is for beam energy. The choices to enter are B field, RF frequency, Momentum (per nucleon), Kinetic Energy (per nucleon)."
        )
        gamText.setWordWrap(True)
        gamLO = QGridLayout()
        gamLO.addWidget(gamText)
        gamWid.setLayout(gamLO)
        self.helpWid.addTab(gamWid, "Gamma")

        buckWid = QWidget()
        buckText = QLabel(
            "Clicking the 'Bucket OK' or 'Bunch OK' button will calculate the bucket parameters and draw the bucket in phase space. Clicking on the ++W and --W buttons will change the vertical scale of the phase space."
        )
        buckText.setWordWrap(True)
        buckLO = QGridLayout()
        buckLO.addWidget(buckText)
        buckWid.setLayout(buckLO)
        self.helpWid.addTab(buckWid, "Bucket OK/Bunch OK")

        bunWid = QWidget()
        bunText = QLabel(
            "Determine how to specify the bunch parameters, either the bunch length or the bunch area."
        )
        bunText.setWordWrap(True)
        bunLO = QGridLayout()
        bunLO.addWidget(bunText)
        bunWid.setLayout(bunLO)
        self.helpWid.addTab(bunWid, "Bunch Length")

        phaseWid = QWidget()
        phaseTitle = QLabel("How to explore phase space?")
        phaseTitle.setFont(self.title_font)
        phaseText = QLabel(
            "The phase space coordinates are RF Phase vs. W. RF phase is measured in degrees (the vertical cursor) and W is the dE/2pi/frf (horizontal cursor) and thus measured in eVs (electron-volt second).\nWhat you initally see is a bucket and a bunch in AGS with some random parameters. The bunch length is also displayed in ns.\nWhen you move the cursor line the vertical cursor will record the RF phase. The horizontal cursor will record the dE, dP, dE/Es, dP/Ps, df/f, and dR/R.\nWhen you click the right mouse button and the vertical cursor is in the right region, the bunch will be redraw with the new phase. The bunch area is also recalculated.\nThe last item show is the voltage, which corresponds to the position of the cursor line. If you want a bunch with that phase and that bunch height, you need to supply the voltage as shown."
        )
        phaseText.setWordWrap(True)
        phaseLO = QGridLayout()
        phaseLO.addWidget(phaseText)
        phaseWid.setLayout(phaseLO)
        self.helpWid.addTab(phaseWid, "Phase Space")

        self.helpWid.show()

    def show_new_window(self, checked):
        if self.w is None:
            self.w = self.MoreResults()
            self.w.show()
        else:
            self.w.close()
            self.w = None

    def textBox(self):
        self.mainTextWid = QWidget()
        self.mainTextLayout = QGridLayout()
        phase = QLabel("phase(deg):")
        self.phase_value = QLabel("54.168")
        phasens = QLabel("phase(ns):")
        self.phasens_value = QLabel("2.459")
        self.mainTextLayout.addWidget(phase, 0, 0, 1, 1)
        self.mainTextLayout.addWidget(self.phase_value, 0, 1, 1, 1)
        self.mainTextLayout.addWidget(phasens, 0, 2, 1, 1)
        self.mainTextLayout.addWidget(self.phasens_value, 0, 3, 1, 1)

        bunch = QLabel("bunch length (ns): ")
        bunch_value = QLabel("2")
        self.mainTextLayout.addWidget(bunch, 1, 0, 1, 2)
        self.mainTextLayout.addWidget(bunch_value, 1, 1, 1, 2)

        de_mev = QLabel("dE(MeV):")
        self.demev_line = QLabel(str(self.dE))
        de_es = QLabel("dE/Es(10^-3):")
        self.dees_line = QLabel(str(self.dE_Es))
        self.mainTextLayout.addWidget(de_mev, 2, 0, 1, 1)
        self.mainTextLayout.addWidget(self.demev_line, 2, 1, 1, 1)
        self.mainTextLayout.addWidget(de_es, 2, 2, 1, 1)
        self.mainTextLayout.addWidget(self.dees_line, 2, 3, 1, 1)

        dp_mevc = QLabel("dp(MeV/c):")
        self.dpmevc_val = QLabel(str(self.dP))
        dp_ps = QLabel("dP/Ps(10^-3)")
        self.dpps_val = QLabel(str(self.dP_Ps))
        self.mainTextLayout.addWidget(dp_mevc, 3, 0, 1, 1)
        self.mainTextLayout.addWidget(self.dpmevc_val, 3, 1, 1, 1)
        self.mainTextLayout.addWidget(dp_ps, 3, 2, 1, 1)
        self.mainTextLayout.addWidget(self.dpps_val, 3, 3, 1, 1)

        dr_r = QLabel("dR/R(10^-3):")
        self.drr_val = QLabel(str(self.dR_R))
        self.dff_label = QLabel("df/f(10^3):")
        self.dff_val = QLabel(str(self.df_f))
        self.mainTextLayout.addWidget(dr_r, 4, 0, 1, 1)
        self.mainTextLayout.addWidget(self.drr_val, 4, 1, 1, 1)
        self.mainTextLayout.addWidget(self.dff_label, 4, 2, 1, 1)
        self.mainTextLayout.addWidget(self.dff_val, 4, 3, 1, 1)

        fnu = QLabel("fnu(Hz):")
        self.fnu_val = QLabel(str(self.fnu))
        self.mainTextLayout.addWidget(fnu, 5, 0, 1, 1)
        self.mainTextLayout.addWidget(self.fnu_val, 5, 1, 1, 1)
        dfnu_label = QLabel("dfnu(Hz):")
        self.dfnu_val = QLabel(str(self.dfnu))
        self.mainTextLayout.addWidget(dfnu_label, 5, 2, 1, 1)
        self.mainTextLayout.addWidget(self.dfnu_val, 5, 3, 1, 1)

        vrfkv = QLabel("Vrf(kV)")
        self.vrfkv_value = QLabel(str(self.Vrf_i))
        self.mainTextLayout.addWidget(vrfkv, 6, 0, 1, 2)
        self.mainTextLayout.addWidget(self.vrfkv_value, 6, 2, 1, 2)

        self.mainTextWid.setLayout(self.mainTextLayout)

        return self.mainTextWid

    def ControlPanel(self):
        cp_widget = QWidget()
        layout = QGridLayout()

        headerFont = QFont()
        headerFont.setBold(True)

        mb_header = QLabel("Editing machine and beam parameters")
        mb_header.setFont(headerFont)
        layout.addWidget(mb_header, 0, 0, 1, 3)

        machine_label = QLabel("Machine")
        layout.addWidget(machine_label, 1, 0, 1, 2)
        self.machine_line = QComboBox()
        self.machine_line.addItems(["AGS", "Booster", "RHIC"])
        self.machine_line.setCurrentIndex(2)
        self.machine_line.currentTextChanged.connect(self.updateMachine)
        layout.addWidget(self.machine_line, 1, 2, 1, 1)

        bd_header = QLabel("B Dot (T/s)")
        bd_header.setFont(headerFont)
        layout.addWidget(bd_header, 2, 0, 1, 2)
        self.bd_line = QLineEdit("0.0")
        self.bd_line.returnPressed.connect(self.setBdot)
        layout.addWidget(self.bd_line, 2, 2, 1, 1)

        self.rff_line = QLineEdit("28.0")
        # drop down for rf frequency
        self.rf_items = [
            "RF Frequency (MHz)",
            "B field (T)",
            "Gamma",
            "pc (GeV/c/u)",
            "K (MeV/u)",
            "K (GeV/u)",
        ]
        self.rff_header = QComboBox()
        self.rff_header.currentTextChanged.connect(self.RFFreq_action)
        self.rff_header.addItems(self.rf_items)

        layout.addWidget(self.rff_header, 3, 0, 1, 2)
        layout.addWidget(self.rff_line, 3, 2, 1, 1)

        self.rfv_header = QLabel("RF voltage/turn (kV)")
        layout.addWidget(self.rfv_header, 4, 0, 1, 2)
        self.rfv_line = QLineEdit("100")
        self.rfv_line.returnPressed.connect(self.setRFV)
        layout.addWidget(self.rfv_line, 4, 2, 1, 1)

        self.rfh_header = QLabel("RF harmonic number")
        layout.addWidget(self.rfh_header, 5, 0, 1, 2)
        self.rfh_line = QLineEdit("360")
        self.rfh_line.returnPressed.connect(self.setRFH)
        layout.addWidget(self.rfh_line, 5, 2, 1, 1)

        w_minus = QPushButton("--W")
        w_minus.clicked.connect(self.W_minus)
        bucket_ok = QPushButton("Bucket OK")
        bucket_ok.clicked.connect(self.calcBucket)
        w_plus = QPushButton("++W")
        w_plus.clicked.connect(self.W_plus)
        layout.addWidget(w_minus, 6, 0, 1, 1)
        layout.addWidget(bucket_ok, 6, 1, 1, 1)
        layout.addWidget(w_plus, 6, 2, 1, 1)

        species_label = QLabel("Species")
        layout.addWidget(species_label, 7, 0, 1, 2)
        self.species_line = QComboBox()
        self.species_line.addItems(
            [
                "Proton (p)",
                "Deuteron (d)",
                "Gold (au)",
                "Iron (fe)",
                "Carbon (c)",
                "Sulfur (s)",
                "Silicon (si)",
                "Copper (cu)",
                "Iodine (i)",
                "Uranium (u)",
            ]
        )

        self.species_line.setCurrentIndex(2)
        self.species_line.currentTextChanged.connect(self.updateCharge)
        layout.addWidget(self.species_line, 7, 2, 1, 1)
        self.charge_header = QLabel("Charge State")
        self.charge_header.setFont(headerFont)
        layout.addWidget(self.charge_header, 8, 0, 1, 2)
        self.charge_line = QLineEdit("79")
        layout.addWidget(self.charge_line, 8, 2, 1, 1)

        # bunch_header = QLabel("Bunch Length (ns)")
        self.bunch_header = QComboBox()
        self.bunch_header.addItems(
            [
                "Bunch Length (ns)",
                "Bunch Emittance (eVs/u)",
                "Bunch Length (rad)",
                "Bunch Length (deg)",
            ]
        )
        self.bunch_header.currentTextChanged.connect(self.updateBunchValue)
        layout.addWidget(self.bunch_header, 9, 0, 1, 2)
        self.bunch_line = QLineEdit("2")
        layout.addWidget(self.bunch_line, 9, 2, 1, 1)

        bw_minus = QPushButton("--W")
        bw_minus.clicked.connect(self.W_minus)
        bunch_ok = QPushButton("Bunch OK")
        bunch_ok.clicked.connect(self.CBunch)
        bw_plus = QPushButton("++W")
        bw_plus.clicked.connect(self.W_plus)
        layout.addWidget(bw_minus, 10, 0, 1, 1)
        layout.addWidget(bunch_ok, 10, 1, 1, 1)
        layout.addWidget(bw_plus, 10, 2, 1, 1)

        self.synchPhase = QLabel("Synchronous Phase (deg)")
        self.synchPhase.setFont(headerFont)
        layout.addWidget(self.synchPhase, 11, 0, 1, 2)
        self.synchLine = QLineEdit("0.0")
        layout.addWidget(self.synchLine, 11, 2, 1, 1)

        self.statBuck = QLabel("Stationary Bucket Area (eVs/u)")
        self.statBuck.setFont(headerFont)
        layout.addWidget(self.statBuck, 12, 0, 1, 2)
        self.statLine = QLineEdit("0.386")
        self.statLine.returnPressed.connect(self.setSTbuck)
        layout.addWidget(self.statLine, 12, 2, 1, 1)

        self.movBuck = QLabel("Moving Bucket Area (eVs/u)")
        self.movBuck.setFont(headerFont)
        layout.addWidget(self.movBuck, 13, 0, 1, 2)
        self.movLine = QLineEdit("0.386")
        layout.addWidget(self.movLine, 13, 2, 1, 1)

        self.synchFreq = QLabel("Synchronous Frequency (Hz)")
        self.synchFreq.setFont(headerFont)
        layout.addWidget(self.synchFreq, 14, 0, 1, 2)
        self.synchFreq = QLineEdit("116.76")
        layout.addWidget(self.synchFreq, 14, 2, 1, 1)

        self.bunchArea = QLabel("Bunch Area (eVs/u)")
        self.bunchArea.setFont(headerFont)
        layout.addWidget(self.bunchArea, 15, 0, 1, 2)
        self.bunchLine = QLineEdit("0.002")
        layout.addWidget(self.bunchLine, 15, 2, 1, 1)

        self.bunchLength = QLabel("Bunch Length (ns)")
        self.bunchLength.setFont(headerFont)
        layout.addWidget(self.bunchLength, 16, 0, 1, 2)
        self.LengthLine = QLineEdit("2")
        layout.addWidget(self.LengthLine, 16, 2, 1, 1)

        cp_widget.setLayout(layout)
        return cp_widget

    def setBdot(self):
        self.Bdot = float(self.bd_line.text())

    def setRFV(self):
        self.Vrf = float(self.rfv_line.text())

    def setRFH(self):
        self.h = float(self.rfh_line.text())

    def updateMachine(self):
        self.machine = self.machine_line.currentText()

    def updateBunchValue(self):
        b_value = self.bunch_header.currentText()
        b_dict = {
            "Bunch Length (ns)": "lt",
            "Bunch Emittance (eVs/u)": "em",
            "Bunch Length (rad)": "lr",
            "Bunch Length (deg)": "ld",
        }
        self.BUN = b_dict[b_value]

    def updateBunchVar(self):
        if self.BUN == "lt":
            self.BUNlt = float(self.bunch_line.text())
            self.LengthLine.setText(str(self.BUNlt))
        elif self.BUN == "em":
            self.BUNe = float(self.bunch_line.text())
            self.LengthLine.setText(str(self.BUNe))
        elif self.BUN == "lr":
            self.BUNlr = float(self.bunch_line.text())
            self.LengthLine.setText(str(self.BUNlr))
        else:
            self.BUNld = float(self.bunch_line.text())
            self.LengthLine.setText(str(self.BUNld))

    def RFFreq_action(self):
        text = self.rff_header.currentText()
        value = self.rf_dict[text]
        # calculate something here?
        self.rff_line.setText("{:.3f}".format(value))
        self.gbfpk = self.gbfpk_name[text]
        # figure out what the actual value should be
        # self.rf_dict[text] = float(self.rff_line.text())
        # for key, value in self.rf_dict.items():
        #    print(key, value)

    def W_plus(self):
        # connects to W++
        self.Wmax = self.Wmax + 2
        self.Wmin = self.Wmin - 2
        self.plot.plot_vb.setYRange(self.Wmin, self.Wmax)
        # configure the y axis to have the new max and min values

    def W_minus(self):
        # connects to W--
        self.Wmax = self.Wmax - 2
        self.Wmin = self.Wmin + 2
        self.plot.plot_vb.setYRange(self.Wmin, self.Wmax)
        # configure the y axis to have the new max and min values

    def updateCharge(self):
        self.species = self.species_line.currentText()
        self.e = self.charge_dict[self.species]
        self.charge_line.setText("{:.2f}".format(self.e))

    def checkSpeciesCharge(self):
        if float(self.charge_line.text()) < 0:
            bTools.errorBox("atomic number is always > 0")

    def calcBucket(self):
        """From cBKT.tcl"""
        if self.e > bmath.Q[self.species]:
            bTools.errorBox("the charge state cannot be greater than the atomic number")

        machineValues = bmath.machineParameters[self.machine]
        self.Eo = bmath.Eo[self.species] * bmath.mega * bmath.A[self.species]
        self.gamma_tr = machineValues["gamma_tr"]
        self.C = machineValues["C"]
        self.Ro = self.C / (2 * bmath.pi)
        self.rho = machineValues["rho"]
        self.Vrf_k = self.Vrf * bmath.kV

        if self.gbfpk == "gamma":
            if self.gammas <= 1.0:
                bTools.errorBox("gammas is always > 1")
            else:
                self.gammas_1 = self.gammas
        elif self.gbfpk == "bf":
            if self.Bf <= 0:
                bTools.errorBox("BF is always > 0")
            else:
                self.Bf_1 = self.Bf
                self.gammas_1 = math.hypot(
                    1.0, self.e * bmath.c * self.rho / self.Eo * self.Bf_1
                )
        elif self.gbfpk == "frf":
            if self.frf <= 0:
                bTools.errorBox("frf is always > 0")
            else:
                self.frfmax = self.h * bmath.c / self.Ro / (2.0 * bmath.pi) / bmath.mega
                if self.frfmax < self.frf:
                    bTools.errorBox(
                        "the rf frequency can't be greater than the frfmax MHz"
                    )
                self.frf_1 = self.frf * bmath.mega
                self.betas = 2 * bmath.pi * self.frf_1 * self.Ro / (self.h * bmath.c)
                self.gammas_1 = 1.0 / math.sqrt(1 - self.betas * self.betas)
        elif self.gbfpk == "pc":
            if self.pc <= 0.0:
                bTools.errorBox("the pc is always > 0")
            self.pc_1 = self.pc * bmath.giga * bmath.A[self.species]
            self.gammas_1 - math.hypot(1.0, self.pc_1 / self.Eo)
        elif self.gbfpk == "km":
            if self.Ekm <= 0:
                bTools.errorBox("Ek is always > 0")
            self.Ek_1 = self.Ekm * bmath.mega * bmath.A[self.species]
            self.gammas_1 = 1.0 + self.Ek_1 / self.Eo
        elif self.gbfpk == "kg":
            if self.Ek <= 0:
                bTools.errorBox("Ek is always > 0")
            self.Ek_1 = self.Ekg * bmath.giga * bmath.A[self.species]
            self.gammas_1 = 1.0 + self.Ek_1 / self.Eo

        self.Es = self.Eo * self.gammas_1
        self.Ek = (self.gammas_1 - 1.0) * self.Eo
        self.pc = math.sqrt(1.0 - self.Eo * self.Eo / self.Es / self.Es) * self.Es
        self.betas = math.sqrt(1.0 - 1.0 / self.gammas_1 / self.gammas_1)
        self.etas = (1 / (self.gamma_tr**2)) - (1 / (self.gammas_1**2))
        self.A = self.etas * ((self.h * bmath.c / self.Ro) ** 2) / self.Es
        self.B = self.e * self.Vrf_k / (2 * bmath.pi * self.h)
        self.phis = bTools.Phi_s(self.etas, self.C, self.rho, self.Bdot, self.Vrf_k)
        if self.phis == float("inf"):
            bTools.errorBox("phis is infinity, ignoring the value")
            self.phis = 0
        self.fnu = math.sqrt(abs(self.A * self.B * math.cos(self.phis))) / (
            2.0 * bmath.pi
        )

        # stationary bucket area
        self.st_bkt = bmath.abkt_s(self.A, self.B)
        # stationary bucket area/u
        self.st_bkt = self.st_bkt / bmath.A[self.species]

        self.alpha = bmath.alpha_bkt(self.phis)

        # moving bucket area
        self.m_bkt = self.alpha * self.st_bkt
        self.m_bkt = self.m_bkt

        self.aphi2 = bTools.phi_2_bkt(self.phis)
        self.aphi1 = bTools.phi_1_bkt(self.phis)
        self.aphi12 = abs(self.aphi2 - self.aphi1)

        # DrawBkt $phis $A $B
        bkt_data = bTools.DrawBkt(self.phis, self.A, self.B)
        self.BKT_x = bkt_data[::2]
        self.BKT_x = self.BKT_x[12:]
        self.BKT_x = self.BKT_x[:-12]
        self.BKT_y = bkt_data[1::2]
        self.BKT_y = self.BKT_y[12:]
        self.BKT_y = self.BKT_y[:-12]

        self.phis_1 = self.phis * 180 / bmath.pi

        if self.gbfpk == "gamma":
            self.Bf_1 = math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "bf":
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "frf":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * self.Eo
            )
        elif self.gbfpk == "pc":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "km":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "kg":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C

        self.frf_1 = self.frf_1 / bmath.mega
        self.frf = self.frf_1
        self.gammas = self.gammas_1
        self.Bf = self.Bf_1
        self.A_1 = self.A
        self.B_1 = self.B
        self.betas_1 = self.betas
        self.fnu_1 = self.fnu
        self.etas_1 = self.etas
        self.Ek_1 = self.Ek / bmath.giga / bmath.A[self.species]
        self.Ek = self.Ek / bmath.giga / bmath.A[self.species]
        self.Ekm = self.Ek * bmath.kilo
        self.Ekg = self.Ek
        self.pc_1 = self.pc / bmath.giga / bmath.A[self.species]
        self.pc = self.pc / bmath.giga / bmath.A[self.species]
        self.BKTld = abs(self.aphi12) * bmath.radeg
        self.BKTlt = abs(self.aphi12) / (2 * bmath.pi * self.frf) * bmath.kilo
        self.BKTdW = math.sqrt(
            abs(
                2
                * self.B_1
                / self.A_1
                * (
                    (bmath.pi - 2 * self.phis) * math.sin(self.phis)
                    - 2 * math.cos(self.phis)
                )
            )
        )

        self.bline.setText("{:.3f}".format(self.Bf))
        self.frfline.setText("{:.3f}".format(self.frf))
        self.gline.setText("{:.3f}".format(self.gammas))
        self.betaline.setText("{:.3f}".format(self.betas_1))
        self.eline.setText("{:.3f}".format(self.etas_1))
        self.keline.setText("{:.3f}".format(self.Ek))
        self.bucLine.setText("{:.3f}".format(self.BKTld))
        self.nsLine.setText("{:.3f}".format(self.BKTlt))
        self.momline.setText("{:.3f}".format(self.pc))
        self.bunLine.setText("{:.3f}".format(self.BUNld))
        self.rf_dict = {
            "RF Frequency (MHz)": self.frf,
            "B field (T)": self.Bf,
            "Gamma": self.gammas,
            "pc (GeV/c/u)": self.pc,
            "K (MeV/u)": self.km,
            "K (GeV/u)": self.kg,
        }
        self.drawBKT_BUN()

    def drawBKT_BUN(self):

        self.plot.addOrUpdateDataset(
            "BKT", self.BKT_x, self.BKT_y, color="red", width=0.8
        )
        self.plot.addOrUpdateDataset(
            "BUN", self.BUN_x, self.BUN_y, color="blue", width=0.8
        )

    def CBunch(self):
        self.checkSpeciesCharge()
        self.updateMachine()
        machineValues = bmath.machineParameters[self.machine]
        self.Eo = bmath.Eo[self.species] * bmath.mega * bmath.A[self.species]
        self.gamma_tr = machineValues["gamma_tr"]
        self.C = machineValues["C"]
        self.Ro = self.C / (2 * bmath.pi)
        self.rho = machineValues["rho"]
        self.Vrf_k = self.Vrf * bmath.kV

        if self.e > bmath.Q[self.species]:
            bTools.errorBox("the charge state cannot be greater than the atomic number")

        if self.gbfpk == "gamma":
            if self.gammas <= 1.0:
                bTools.errorBox("gammas is alway > 1")
            self.gammas_1 = self.gammas
        elif self.gbfpk == "bf":
            if self.Bf <= 0:
                bTools.errorBox("Bf is always > 0")
            self.Bf_1 = self.Bf
            self.gammas_1 = math.hypot(
                1.0, (self.e * bmath.c * self.rho / self.Eo * self.Bf_1)
            )
        elif self.gbfpk == "frf":
            if self.frf <= 0.0:
                bTools.errorBox("frf is always > 0")
            self.frfmax = self.h * bmath.c / self.Ro / (2 * bmath.pi) / bmath.mega
            if self.frfmax < self.frf:
                bTools.errorBox("error rf frequency can't be greater than the maxmimum")
            self.frf_1 = self.frf * bmath.mega
            self.betas = 2 * bmath.pi * self.frf_1 * self.Ro / (self.h * bmath.c)
            self.gammas_1 = 1.0 / math.sqrt(1.0 - self.betas * self.betas)
        elif self.gbfpk == "pc":
            if self.pc <= 0.0:
                bTools.errorBox("the pc is always > 0")
            self.pc_1 = self.pc * bmath.giga * bmath.A[self.species]
            self.gammas_1 = math.hypot(1.0, self.pc_1 / self.Eo)
        elif self.gbfpk == "km":
            if self.Ekm <= 0:
                bTools.errorBox("Ek is always > 0")
            self.Ek_1 = self.Ekm * bmath.mega * bmath.A[self.species]
            self.gammas_1 = 1.0 + self.Ek_1 / self.Eo
        elif self.gbfpk == "kg":
            if self.Ek <= 0:
                bTools.errorBox("Ek is always > 0")
            self.Ek_1 = self.Ekg * bmath.giga * bmath.A[self.species]
            self.gammas_1 = 1.0 + self.Ek_1 / self.Eo

        self.Es = self.Eo * self.gammas_1
        self.Ek = (self.gammas_1 - 1.0) * self.Eo
        self.pc = math.sqrt(1.0 - self.Eo * self.Eo / self.Es / self.Es) * self.Es
        self.betas = math.sqrt(1.0 - 1.0 / self.gammas_1 / self.gammas_1)
        self.etas = (1 / (self.gamma_tr**2)) - (1 / (self.gammas_1**2))
        self.A = self.etas * (self.h * bmath.c / self.Ro) ** 2.0 / self.Es
        self.B = self.e * self.Vrf_k / (2 * bmath.pi * self.h)
        self.phis_2 = bTools.Phi_s(self.etas, self.C, self.rho, self.Bdot, self.Vrf_k)
        self.phis = self.phis_2
        self.phis_1 = self.phis_2 * 180 / bmath.pi
        self.st_bkt = 16.0 * math.sqrt(abs(self.B / self.A))
        self.alpha = bTools.Alpha_bkt(self.phis_2)
        self.m_bkt = self.st_bkt * self.alpha

        if self.gbfpk == "gamma":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * self.Eo
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "bf":
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "frf":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * self.Eo
            )
        elif self.gbfpk == "pc":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "km":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        elif self.gbfpk == "kg":
            self.Bf_1 = (
                math.sqrt(self.gammas_1 * self.gammas_1 - 1.0)
                / self.e
                / bmath.c
                / self.rho
                * (self.Eo)
            )
            self.frf_1 = self.h * self.betas * bmath.c / self.C
        self.updateBunchVar()
        if self.BUN == "em":
            if self.BUNe <= 0.0:
                bTools.errorBox(" Bunch area is always > 0")
            self.BUNe_1 = self.BUNe * bmath.A[self.species]
            self.abun = self.BUNe_1 / self.st_bkt
            self.phibun2 = bTools.i_Alpha_bun(self.abun, self.phis_2)
            self.phibun1 = bTools.Phi_1_bun(self.phis_2, self.phibun2)
            self.phibun1 = bTools.Proper_phi(self.phis_2, self.phibun1)
            self.phibun2 = bTools.Proper_phi(self.phis_2, self.phibun2)
            self.phi12 = abs(self.phibun2 - self.phibun1)
            self.BUNld = self.phi12 * bmath.radeg
        if self.BUN == "lt":
            if self.BUNlt <= 0:
                bTools.errorBox("Bunch length is always > 0")
            self.phi12 = bmath.ns * 2.0 * bmath.pi * self.frf_1 * self.BUNlt
            self.BUNld = self.phi12 * bmath.radeg
            self.BUNlr = self.phi12
            self.abun = bTools.Alpha_bun_phi12(self.phi12, self.phis_2)
            self.BUNe_1 = self.abun * self.st_bkt
            self.phibun2 = bTools.i_Alpha_bun(self.abun, self.phis_2)
            self.phibun1 = bTools.Phi_1_bun(self.phis_2, self.phibun2)
            self.phibun1 = bTools.Proper_phi(self.phis_2, self.phibun1)
            self.phibun2 = bTools.Proper_phi(self.phis_2, self.phibun2)
        if self.BUN == "lr":
            if self.BUNlr <= 0.0:
                bTools.errorBox(" Bunch length is always > 0")
            self.phi12 = self.BUNlr
            self.BUNlt = self.BUNlr / (2.0 * bmath.pi * self.frf_1) / bmath.ns
            self.BUNld = self.phi12 * bmath.radeg
            self.abun = bTools.Alpha_bun_phi12(self.phi12, self.phis_2)
            self.BUNe_1 = self.abun * self.st_bkt
            self.phibun2 = bTools.i_Alpha_bun(self.abun, self.phis_2)
            self.phibun1 = bTools.Phi_1_bun(self.phis_2, self.phibun2)
            self.phibun1 = bTools.Proper_phi(self.phis_2, self.phibun1)
            self.phibun2 = bTools.Proper_phi(self.phis_2, self.phibun2)
        if self.BUN == "ld":
            if self.BUNld <= 0.0:
                bTools.errorBox(" Bunch length is always > 0")
            self.phi12 = self.BUNld / bmath.radeg
            self.BUNlr = self.phi12
            self.BUNld = self.phi12 * bmath.radeg
            self.BUNlt = self.phi12 / (2.0 * bmath.pi * self.frf_1) / bmath.ns
            self.abun = bTools.Alpha_bun_phi12(self.phi12, self.phis_2)
            self.BUNe_1 = self.abun * self.st_bkt
            self.phibun2 = bTools.i_Alpha_bun(self.abun, self.phis_2)
            self.phibun1 = bTools.Phi_1_bun(self.phis_2, self.phibun2)
            self.phibun1 = bTools.proper_phi(self.phis_2, self.phibun1)
            self.phibun2 = bTools.proper_phi(self.phis_2, self.phibun2)

        bun_data = bTools.DrawBun(
            self.phis_2, self.phibun1, self.phibun2, self.A, self.B
        )
        self.BUN_x = bun_data[::2]
        self.BUN_y = bun_data[1::2]
        self.plot.addOrUpdateDataset(
            "BUN", self.BUN_x, self.BUN_y, color="blue", width=0.8
        )

        if self.BUN == "em":
            # in unit of ns
            self.phi12_1 = abs(self.phibun2 - self.phibun1)
            self.phi12t_1 = (
                abs(self.phibun2 - self.phibun1)
                / (2 * bmath.pi * self.frf_1)
                / bmath.ns
            )
            self.BUNlt = self.phi12t_1
        if self.BUN == "lt":
            # in unit of ns
            self.phi12_1 = self.phi12
            self.phi12t_1 = self.BUNlt
        if self.BUN == "lr":
            # in unit of ns
            self.phi12_1 = self.phi12
            self.phi12t_1 = self.BUNlt
        if self.BUN == "ld":
            # in unit of ns
            self.phi12_1 = self.phi12
            self.phi12t_1 = self.BUNlt

        # tie up lose ends
        self.frf_1 = self.frf_1 / bmath.mega
        self.frf = self.frf_1
        self.Bf = self.Bf_1
        self.A_1 = self.A
        self.B_1 = self.B
        self.BUNe_1 = self.BUNe_1 / bmath.A[self.species]
        self.BUNe = self.BUNe_1
        self.gammas = self.gammas_1
        self.betas_1 = self.betas
        self.st_bkt = self.st_bkt / bmath.A[self.species]
        self.m_bkt = self.m_bkt / bmath.A[self.species]
        self.fnu = math.sqrt(abs(self.A * self.B * math.cos(self.phis))) / (
            2.0 * bmath.pi
        )
        self.etas_1 = self.etas
        self.Ek_1 = self.Ek / bmath.giga / bmath.A[self.species]
        self.Ek = self.Ek / bmath.giga / bmath.A[self.species]
        self.Ekm = self.Ek * bmath.kilo
        self.Ekg = self.Ek
        self.pc_1 = self.pc / bmath.giga / bmath.A[self.species]
        self.pc = self.pc / bmath.giga / bmath.A[self.species]

        self.bline.setText("{:.3f}".format(self.Bf))
        self.frfline.setText("{:.3f}".format(self.frf))
        self.gline.setText("{:.3f}".format(self.gammas))
        self.betaline.setText("{:.3f}".format(self.betas_1))
        self.eline.setText("{:.3f}".format(self.etas_1))
        self.keline.setText("{:.3f}".format(self.Ek))
        self.bucLine.setText("{:.3f}".format(self.BKTld))
        self.nsLine.setText("{:.3f}".format(self.BKTlt))
        self.momline.setText("{:.3f}".format(self.pc))
        self.bunLine.setText("{:.3f}".format(self.BUNld))

    def bltCoor(self, x, y):
        self.phase_t1 = x / 360 / (self.frf_1 + 1.0e-20) * bmath.kilo
        self.phase_value.setText("{:.3f}".format(self.phase_t1))
        self.dE = y * 2.0 * bmath.pi * self.frf_1
        self.demev_line.setText("{:.3f}".format(self.dE))
        self.dE_Es = self.dE / self.Es * bmath.kilo * bmath.mega
        self.dees_line.setText("{:.3f}".format(self.dE_Es))
        self.dP = self.dE / self.betas_1
        self.dpmevc_val.setText("{:.3f}".format(self.dP))
        self.dP_Ps = self.dE_Es / (self.betas_1 * self.betas_1)
        self.dpps_val.setText("{:.3f}".format(self.dP_Ps))
        self.gamma_tr = self.machineValues["gamma_tr"]
        self.dR_R = self.dP_Ps / self.gamma_tr**2
        self.drr_val.setText("{:.3f}".format(self.dR_R))
        self.df_f = -self.dP_Ps * self.etas_1
        self.dff_val.setText("{:.3f}".format(self.df_f))
        self.x = x * (bmath.pi / 180)
        self.Tnu = bTools.Period_bun(self.phis, self.x, self.A_1, self.B_1)
        self.fnu_1 = 1.0 / (self.Tnu + 1.0e-20)
        self.fnu_val.setText("{:.3f}".format(self.fnu_1))
        self.dfnu = -(self.fnu - self.fnu_1)
        self.dfnu_val.setText("{:.3f}".format(self.dfnu))
        self.Vrf_i = self.A_1 * (y**2) * bmath.pi * self.h
        self.Vrf_i = (
            abs(self.Vrf_i)
            / self.e
            / bmath.kilo
            / (
                1.0e-20
                + abs(
                    math.cos(self.x)
                    - math.cos(self.phis)
                    + (self.x - self.phis) * math.sin(self.phis)
                )
            )
        )
        self.vrfkv_value.setText("{:.3f}".format(self.Vrf_i))

    def bun(self, x, y):
        self.phibun2 = x / bmath.radeg
        self.phi2 = bmath.pi - self.phis  # phi2 of bkt

        if self.etas_1 >= 0:
            if self.phibun2 <= self.phi2 and self.phibun2 > self.phis:
                self.phibun1 = bTools.Phi_1_bun(self.phis, self.phibun2)
                self.BUNlr = abs(self.phibun2 - self.phibun1)
                self.phi12 = self.BUNlr
                self.BUNlt = (
                    self.BUNlr / (2.0 * bmath.pi * self.frf * bmath.mega) / bmath.ns
                )
                self.BUNld = self.phi12 * bmath.radeg
                self.abun = bTools.Alpha_bun_phi12(self.phi12, self.phis)
                self.BUNe_1 = self.abun * self.st_bkt
                self.BUNe = self.BUNe_1
                bun_data = bTools.DrawBun(
                    self.phis, self.phibun1, self.phibun2, self.A_1, self.B_1
                )
                self.BUN_x = bun_data[::2]
                self.BUN_y = bun_data[1::2]
        else:
            if self.phibun2 >= self.phi2 and self.phibun2 < self.phis:
                self.phibun2 = bTools.Generic_phi2(self.phis, self.phibun2)
                self.phibun1 = bTools.Phi_1_bun(self.phis, self.phibun2)
                self.BUNlr = abs(self.phibun2 - self.phibun1)
                self.phi12 = self.BUNlr
                self.BUNlt = (
                    self.BUNlr / (2.0 * bmath.pi * self.frf_1 * bmath.mega) / bmath.ns
                )  # MHz
                self.BUNld = self.phi12 * bmath.radeg
                self.abun = bTools.Alpha_bun_phi12(self.phi12, self.phis)
                self.BUNe_1 = self.abun * self.st_bkt
                self.BUNe = self.BUNe_1
                self.phibun2 = bTools.i_Alpha_bun(self.abun, self.phis)
                self.phibun1 = bTools.Phi_1_bun(self.phis, self.phibun2)
                self.phibun1 = bTools.proper_phi(self.phis, self.phibun1)
                self.phibun2 = bTools.proper_phi(self.phis, self.phibun2)
                bun_data = bTools.DrawBun(
                    self.phis, self.phibun1, self.phibun2, self.A_1, self.B_1
                )
                self.BUN_x = bun_data[::2]
                self.BUN_y = bun_data[1::2]
