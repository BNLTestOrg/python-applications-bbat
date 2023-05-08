from qtpy.QtWidgets import (
    QLabel,
    QWidget,
    QGridLayout,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSlider, QTabWidget
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from cad_ui.general import CADMainWindow, PrintMenu
from cad_ui.plotting import CadPlot
import numpy as np
import blt
import os
import math
from .bmath import bmath


class BBat(CADMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(title="bbat", menubar=False, *args, **kwargs)
        self.n = 2
        self.Vrf = 0
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vn = 100
        self.Vn_k = self.Vn * bmath.kilo
        self.theta = 86
        self.phis1 = 10
        self.vzerox = [-190, 190]
        self.vzeroy = [0,0]
        self.bdot = (8.7*201*13.7)/1000
        self._g2rfsep = []
        # these are undefined
        #self._g2rflist = [self.RFV2, self.RFU2, self.vzero, self.vinf, self.phis1, self.hline]
        #v1max is in para.py
        self.e = 0
        self.h = 360
        self.A1 = 0
        self.B_1 = 0
        self.B1 = self.e/(2*bmath.pi*self.h) #B1 is not B_1, B1= Vrf*B_1;
        self.sop = 0
        self.A2Bun = 0
        self.fnu2rf = 0
        self.W2Max = 3
        self.W2Min = -3
        self.U2Max = 200
        self.U2Min = 200
        #bltvector 2RFU_x
        self._RFU2_x = []
        #bltvector 2RFU_y
        self._RFU2_y = []
        #bltvector 2RFV_x
        self._RFV2_x = []
        #bltvector 2RFV_y
        self._RFV2_y = []

        #bltvector phis_1_x
        self.phis_1_x = []
        #bltvector phis_1_y
        self.phis_1_y = []
        #bltvector vinf_x
        self.vinf_x = []
        #bltvector vinf_y
        self.vinf_y = []
        #bltvector vzero_x
        self.vzero_x = []
        #bltvector vzero_y
        self.vzero_y = []
        #bltvector hline_x
        self.hline_x = []
        #bltvector hline_y
        self.hline_y = []

        self.dE = 0
        self.frf_1 = 0
        self.dE_Es = 0
        self.Es = 0
        self.phis = 0
        self.fnu1 = 0
        self.dP = 0 
        self.dP_Ps = 0
        self.betas1 = 0
        self.lotheta = -180/self.n
        self.uptheta = 180/self.n
        self.deltav = self.Vn/(self.Vrf + 1e-20)
        self.phi2s = self.phis1 + self.theta

        # these are all commented out
        #self._phibun1 = 0
        #self._phibun2 = 0
        #self._etas1 = 0
        # pi and radeg are in the para
        #self.BUNlt = 0
        #self.BUNlr = 0
        self.species = None
        self.fnu2rf = 0
        self.A = 0
        self.ww = 1
        self.uu = 10

        self.w = None
        self.Bf = 0.30662093657
        self.frf = 28.0
        #self.gammas = 5
        self.gammas = 9.68008184872
        self.pc = 8.96644213239
        self.km = 8083.41302244
        self.kg = 8.08341302244
        self.BDot = 0.0
        self.Vrf = 100
        self.h = 360
        self.Wmax = 2
        self.Wmin = -2
        # some how e is the charge state
        # self._e = 79
        self.chargeState = 79

        self.rf_dict = {"RF Frequency (MHz)":self.frf, "B field (T)":self.Bf, "Gamma":self.gammas, "pc (GeV/c/u)":self.pc, "K (MeV/u)":self.km, "K (GeV/u)":self.kg}
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
                "Ruthenium (ru)": bmath.Q_ru,
                "Zirconium (zr)": bmath.Q_zr,
                "others": 0
            }

        # UI components
        self.setWindowTitle("bbat")
        control = self.ControlPanel()
        buttons = self.ButtonPanel()
        layout = QGridLayout()
        layout.addWidget(control, 0, 0, 17, 3)
        more = QPushButton("More Results")
        more.clicked.connect(self.show_new_window)
        layout.addWidget(more, 17, 1, 1, 1)
        layout.addWidget(buttons, 0, 3, 1, 1)

        self.plot = CadPlot()
        self.plot.plotItem.setLabel("bottom", "RF Phase (deg)")
        self.plot.plotItem.setLabel("left", "W")
        self.plot.plot_vb.setXRange(190, 370)
        self.plot.plot_vb.setYRange(self.Wmin, self.Wmax)
        self.plot.addDataset(
            "Random", np.arange(0, 10, 1), np.arange(0, 10, 1), color="b"
        )
        layout.addWidget(self.plot, 1, 3, 7, 1)

        textbox = self.textBox()
        layout.addWidget(textbox, 8, 3, 4, 1)

        wid = QWidget()
        wid.setLayout(layout)
        self.setCentralWidget(wid)

    def secondRFWindow(self):
        self.srfWid = QWidget()
        self.srfWid.setWindowTitle("Dr. BBat")
        quit = QPushButton("Close")
        quit.clicked.connect(self.srfWid.close)
        redraw = QPushButton("Redraw")
        vv = QPushButton("V")
        vv.clicked.connect(self.vplot)
        printButton = QPushButton("Print")
        pm = PrintMenu(self)
        printButton.setMenu(pm)
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
        layout.addWidget(refresh, 16, 1, 1, 3)

        srfPlot = CadPlot()
        srfPlot.plotItem.setLabel("bottom", "RF Phase (deg)")
        srfPlot.plotItem.setLabel("left", "W")

        # Create the Blt graph widget
        #second_g2rf = Blt.Graph(master, height=400, width=600)
        #second_g2rf.grid(row=0, column=0)

        # Configure the x-axis
        #second_g2rf.xaxis.configure(step=90, rotate=0, command=g2XLabels,min=-370, max=370)

        # Configure the y-axes
        #second_g2rf.yaxis.configure(rotate=0.0, min=W2min, max=W2max)
        #second_g2rf.y2axis.configure(rotate=0.0, min=U2min, max=U2max)

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
        textLayout.addWidget(phase_value, 0, 1, 1, 1)
        textLayout.addWidget(phasens, 0, 2, 1, 1)
        textLayout.addWidget(phasens_value, 0, 3, 1, 1)
        de_mev = QLabel("dE(MeV):")
        demev_line = QLabel("something")
        de_es = QLabel("dE/Es(10^-3):")
        dees_line = QLabel("something")
        textLayout.addWidget(de_mev, 1, 0, 1, 1)
        textLayout.addWidget(demev_line, 1, 1, 1, 1)
        textLayout.addWidget(de_es, 1, 2, 1, 1)
        textLayout.addWidget(dees_line, 1, 3, 1, 1)

        dp_mevc = QLabel("dp(MeV/c):")
        dpmevc_val = QLabel("something")
        dp_ps = QLabel("dP/Ps(10^-3)")
        dpps_val = QLabel("Something")
        textLayout.addWidget(dp_mevc, 2, 0, 1, 1)
        textLayout.addWidget(dpmevc_val, 2, 1, 1, 1)
        textLayout.addWidget(dp_ps, 2, 2, 1, 1)
        textLayout.addWidget(dpps_val, 2, 3, 1, 1)

        dr_r = QLabel("dR/R(10^-3):")
        drr_val = QLabel("something")
        df_f = QLabel("df/f(10^3):")
        dff_val = QLabel("something")
        textLayout.addWidget(dr_r, 3, 0, 1, 1)
        textLayout.addWidget(drr_val, 3, 1, 1, 1)
        textLayout.addWidget(df_f, 3, 2, 1, 1)
        textLayout.addWidget(dff_val, 3, 3, 1, 1)

        fnu = QLabel("fnu(Hz):")
        fnu_val = QLabel(" ")
        textLayout.addWidget(fnu, 4, 0, 1, 1)
        textLayout.addWidget(fnu_val, 4, 1, 1, 1)

        abkt = QLabel("Abkt(eVs/u):")
        abkt_val = QLabel(" ")
        abun = QLabel("Abun(eVs/u):")
        abun_val = QLabel(" ")
        textLayout.addWidget(abkt, 5, 0, 1, 1)
        textLayout.addWidget(abkt_val, 5, 1, 1, 1)
        textLayout.addWidget(abun, 6, 0, 1, 1)
        textLayout.addWidget(abun_val, 6, 1, 1, 1)

        textWid.setLayout(textLayout)
        layout.addWidget(textWid, 20, 5, 4, 4)

        self.srfWid.setLayout(layout)
        return self.srfWid
        #self.srfWid.show()

    #def init_sep(self, n):
    #    i = int(1+math.ceil(n))
        # if self._rfv2 is visible on the graph then
        # for j in self._g2rfsep:
            #try:
            #     del .second.g2rf.element[j]
            #except Exception as e:
            #    print(f"catch: {e}")
            # self._g2rfsep = []
        # while i >= 0:
            #.second.g2rf element create sep$i -xdata sep${i}_x -ydata sep${i}_y -symbol none -linewidth 0.8 -color orange
		    #lappend g2rfsep sep$i
		    #set i [expr $i-1]

    #def refresh_command(self):
        #if {[lsearch [eval .second.g2rf element show] 2RFV] != -1} {
	        #.second.g2rf element delete [concat $g2rflist $g2rfsep] }

        #.second.g2rf xaxis configure -min -180  -max 360
        #.second.g2rf element create 2RFU  -symbol none -xdata 2RFU_x -ydata 2RFU_y  -linewidth 0.8 -color blue \
	        #-mapx x2 -mapy y2
        #.second.g2rf element create 2RFV  -xdata 2RFV_x -ydata 2RFV_y -symbol none -linewidth 0.8 -color green \
	        #-mapx x2 -mapy y2

	    ##############################bdot
        #.second.g2rf element create bdot -symbol line -linewidth 0.8 -fg black  \
        #-xdata $vzerox -ydata {$bdot $bdot}

        ##############################phis_1


        #.second.g2rf element create phis_1 -symbol none -linewidth 0.8 \
        # -color brown -xdata phis_1_x -ydata phis_1_y
        # {$phis_1 $phis_1} -ydata {0 1000}

        ##############################coordinates
        #.second.g2rf element create vzero -symbol none -linewidth 0.8 \
        #-color red -xdata vzero_x -ydata vzero_y
        #.second.g2rf element create vinf -symbol none -linewidth 0.8 \
        #-color red -xdata vinf_x -ydata vinf_y
        #	-xdata vzeroy -ydata $vzerox

        ##############################hlines
        
        #.second.g2rf element create hline -symbol none -linewidth 0.8 -color orange -xdata hline_x -ydata hline_y
            #Init_Sep $n

        #bltvector phis_1_x
        #bltvector phis_1_y
        #for {set i 0} {$i < 100} {incr i} {
        #	bltvector sep${i}_x
        #	bltvector sep${i}_y
        #}

    #def draw2rfButtonCommand(self):
        """ this is the action the redraw button completes
            i believe it resets the values in the rf parameter list
            and redraws the main lines in the graph. (idk about hline, sep3,2,1,0)
        """
        #n_value = self.init_sep(self.n)
        # this creates the buttons for vrf, vn, n , and theta
	    #Draw2RF   ".second.g2rf" "2RFV" $Vrf $Vn $n $theta; 
	    #Draw2rfU  ".second.g2rf" "2RFU" $A_1 $B_1 $Vrf $Vn $n $theta $phis_1; 
	    #Drawphis  ".second.g2rf" "phis_1" $phis_1 0.5;
	    #Draw2rfSep  ".second.g2rf" "sep" $A_1 $B1 $vrf $vn $n $theta $phis_1;

        # self.redrawButton.clicked.connect(self.draw2rfButtonCommand)
        #button .second.menu.draw2rf -text "Redraw" -relief ridge -command draw2rfButtonCommand

        # self.VButton.clicked.connect(launch the other plot command here)
        #button .second.menu.v2rf -text "V" -relief ridge -command {DrawV}

        #self.help.clicked.connect( something about the help menu )
        #menubutton .second.menu.help -text Help -relief ridge -menu .second.menu.help.m 
        #menu .second.menu.help.m

    def helpText(self):
        """This help text will replace the help menu in the widget
        HLP_load .second.menu.help.m "$DIR/cSecondHelp.txt"
        """
        txt = "Dr. BBat is a Double RF Bunch and Bucket Analysis Tool\nCurves\n Green: Total RF Waveform\n"

    def setVrfValue(self, value):
        """This function is the same as Vnvalue (line 250)

        Args:
            value (_type_): _description_
        """
        self.A = value
        self.Vrf = value
        self.Vrf_k = self._Vrf * bmath.kilo
        self.deltav = self.Vn / (self.Vrf+1.e-20)
        self.lotheta = -180/self.n
        self.uptheta = 180/self.n 
        self.phi2s = self.phis1 + self.theta
        # call draw2rfButtonCommand

    def setKValue(self, value):
        self.A = value
        self.Vrf = (self.A * bmath.v1max)/100
        self.Vrf_k = self.Vrf * bmath.kilo
        self.Vn = (100.0-self.A) * bmath.v1max/100
        # .second.menu.draw2rf invoke

    def setTValue(self, value):
        self.A = value
        # .second.menu.draw2rf invoke

    #def keyPressEvent(self, a0: QtGui.QKeyEvent):
    #    if a0 == QtGui.Key.Enter:
    #        self._vn = self._Vn * Para.kilo
    #        self._vrf = self._Vrf * Para.kilo
    #        self._deltav = self._Vn/(self._Vrf+1e-20)
    #        self.init_sep(self._n)
    #        #$input2.ok invoke
	#        #.second.menu.draw2rf invoke
    #    return super().keyPressEvent(a0)

    def W2_plus(self):
        # connects to W++
        self.W2max = self.W2max + self.ww
        self.W2min = self.W2min - self.ww
        # configure the y axis to have the new max and min values

    def W2_minus(self):
        # connects to W--
        self.W2max = self.W2max - self.ww
        self.W2min = self.W2min + self.ww
        # configure the y axis to have the new max and min values

    def U2_plus(self):
        # connects to U++
        self.U2max = self.U2max + self.uu
        self.U2min = self.U2min - self.uu
        # configure the y2 axis to have the new max and min values

    def U2_minus(self):
        # connects to U--
        self.U2max = self.U2max - self.uu
        self.U2min = self.U2min + self.uu
        # configure the y2 axis to have the new max and min values

    def vplot(self):
        self.v_wid = QWidget()
        self.v_wid.setWindowTitle("Vrf form")
        self.v_layout = QGridLayout()
        self.v_plot = CadPlot()
        self.v_layout.addWidget(self.v_plot, 0, 0, 6, 10)
        self.v_ok = QPushButton("OK")
        self.v_ok.clicked.connect(self.v_wid.close)
        self.v_layout.addWidget(self.v_ok, 6, 0, 1, 10)
        self.v_wid.setLayout(self.v_layout)
        self.v_wid.show()

    #########
    # Methods for graph stuff 
    #########
    #def set2rfcoor(self, graph):
        # track mouse motion and add the cursorline
        #self.v_plot.addInfiniteLine ?
        # this might have to be custom, a constantly following mouse line

    def blt2rfCoor(self, graph, winX, winY):
        # scan [$graph invtransform $winX $winY] "%s %s" sop(x) sop(y)
        self.phase_t1 = (self.sop(x)/360) / (self.frf_1 + 1.0e-20)* bmath.kilo
        self.dE = (self.sop(y) * 2.0 * bmath.pi * self.frf_1) #frf_1 in MHz, dE in MeV
        self.dE_Es  = self.dE/self.Es*bmath.kilo* bmath.mega #Es in eV  
        self.dP = self.dE/ self.betas1
        self.dP_Ps = self.betas1*self.betas1* self.dE_Es
        # need to find the gamma tr, idk if it means the current selection
        # set to AGS for now
        self.dR_R = self.dP_Ps/ bmath.gamma_tr_ags**2
        #set dR_R [expr $dP_Ps/pow($gamma_tr,2.0)]
        self.df_f = -self.dP_Ps*self.etas_1

    #proc Set2RFBun { graph } {
    #    global bindings 
    #    set bindings(<ButtonPress-3>,$graph,coor) { 
	#    blt2rfbun %W %x %y
    #}
    #   bltResetBindings $graph <ButtonPress-3>
    #}

    def blt2RfBun(self, graph, winX, winY):
        # phi2 was not defined
        phi2 = self.sop(x)

    def MoreResults(self):
        self.morRes = QWidget()
        self.morRes.setWindowTitle("bkt")
        layout = QGridLayout()
        title = QLabel("Update these variables by clicking on either Bucket OK or Bunch OK")
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

        ok = QPushButton("OK")
        ok.clicked.connect(self.morRes.close)
        layout.addWidget(ok, 11, 0, 1, 2)

        self.morRes.setLayout(layout)
        return self.morRes
        #self.morRes.show()

    def ButtonPanel(self):
        self.buttonWid = QWidget()
        self.title_font = QFont()
        self.title_font.setPointSize(10.5)
        self.title_font.setBold(True)

        quit = QPushButton("Close")
        quit.clicked.connect(self.close)
        refresh = QPushButton("Refresh")
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

    def show_SRF(self, checked):
        self.srf = self.secondRFWindow()
        self.srf.show()

    def helpWidget(self):
        self.helpWid = QTabWidget()
        self.helpWid.setFixedSize(600,300)
        self.helpWid.setWindowTitle("Help")
        overview = QWidget()
        overLayout = QGridLayout()
        title = QLabel("BBAT Version 2 by Jennefer Maldonado (jmaldonad@bnl.gov)")
        title.setFont(self.title_font)
        body_text = QLabel("BBat pronunced b-bat is a Bunch and Bucket Analysis Tool for a single RF system.")
        overLayout.addWidget(title, 0, 0, 1, 3)
        overLayout.addWidget(body_text, 1, 0, 1, 3)
        overview.setLayout(overLayout)
        self.helpWid.addTab(overview, "Overview")

        machine = QWidget()
        machText = QLabel("Press the machine button to list the built in machines. AGS is the default machine. Choosing others, you'll need to key in some machine parameters. Click ok when you finish typing or chaning some other parameters.")
        machText.setWordWrap(True)
        machLO = QGridLayout()
        machLO.addWidget(machText, 0, 0, 1, 1)
        machine.setLayout(machLO)
        self.helpWid.addTab(machine, "Machine")

        gamWid = QWidget()
        gamText = QLabel("This button is for beam energy. The choices to enter are B field, RF frequency, Momentum (per nucleon), Kinetic Energy (per nucleon).")
        gamText.setWordWrap(True)
        gamLO = QGridLayout()
        gamLO.addWidget(gamText)
        gamWid.setLayout(gamLO)
        self.helpWid.addTab(gamWid, "Gamma")

        buckWid = QWidget()
        buckText = QLabel("Clicking the 'Bucket OK' or 'Bunch OK' button will calculate the bucket parameters and draw the bucket in phase space. Clicking on the ++W and --W buttons will change the vertical scale of the phase space.")
        buckText.setWordWrap(True)
        buckLO = QGridLayout()
        buckLO.addWidget(buckText)
        buckWid.setLayout(buckLO)
        self.helpWid.addTab(buckWid, "Bucket OK/Bunch OK")

        bunWid = QWidget()
        bunText = QLabel("Determine how to specify the bunch parameters, either the bunch length or the bunch area.")
        bunText.setWordWrap(True)
        bunLO = QGridLayout()
        bunLO.addWidget(bunText)
        bunWid.setLayout(bunLO)
        self.helpWid.addTab(bunWid, "Bunch Length")

        phaseWid = QWidget()
        phaseTitle = QLabel("How to explore phase space?")
        phaseTitle.setFont(self.title_font)
        phaseText = QLabel("The phase space coordinates are RF Phase vs. W. RF phase is measured in degrees (the vertical cursor) and W is the dE/2pi/frf (horizontal cursor) and thus measured in eVs (electron-volt second).\nWhat you initally see is a bucket and a bunch in AGS with some random parameters. The bunch length is also displayed in ns.\nWhen you move the cursor line the vertical cursor will record the RF phase. The horizontal cursor will record the dE, dP, dE/Es, dP/Ps, df/f, and dR/R.\nWhen you click the right mouse button and the vertical cursor is in the right region, the bunch will be redraw with the new phase. The bunch area is also recalculated.\nThe last item show is the voltage, which corresponds to the position of the cursor line. If you want a bunch with that phase and that bunch height, you need to supply the voltage as shown.")
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
        phase_value = QLabel("54.1689284")
        phasens = QLabel("phase(ns):")
        phasens_value = QLabel("2.45909358")
        self.mainTextLayout.addWidget(phase, 0, 0, 1, 1)
        self.mainTextLayout.addWidget(phase_value, 0, 1, 1, 1)
        self.mainTextLayout.addWidget(phasens, 0, 2, 1, 1)
        self.mainTextLayout.addWidget(phasens_value, 0, 3, 1, 1)

        bunch = QLabel("bunch length (ns): ")
        bunch_value = QLabel("2")
        self.mainTextLayout.addWidget(bunch, 1, 0, 1, 2)
        self.mainTextLayout.addWidget(bunch_value, 1, 1, 1, 2)

        de_mev = QLabel("dE(MeV):")
        demev_line = QLabel("something")
        de_es = QLabel("dE/Es(10^-3):")
        dees_line = QLabel("something")
        self.mainTextLayout.addWidget(de_mev, 2, 0, 1, 1)
        self.mainTextLayout.addWidget(demev_line, 2, 1, 1, 1)
        self.mainTextLayout.addWidget(de_es, 2, 2, 1, 1)
        self.mainTextLayout.addWidget(dees_line, 2, 3, 1, 1)

        dp_mevc = QLabel("dp(MeV/c):")
        dpmevc_val = QLabel("something")
        dp_ps = QLabel("dP/Ps(10^-3)")
        dpps_val = QLabel("Something")
        self.mainTextLayout.addWidget(dp_mevc, 3, 0, 1, 1)
        self.mainTextLayout.addWidget(dpmevc_val, 3, 1, 1, 1)
        self.mainTextLayout.addWidget(dp_ps, 3, 2, 1, 1)
        self.mainTextLayout.addWidget(dpps_val, 3, 3, 1, 1)

        dr_r = QLabel("dR/R(10^-3):")
        drr_val = QLabel("something")
        df_f = QLabel("df/f(10^3):")
        dff_val = QLabel("something")
        self.mainTextLayout.addWidget(dr_r, 4, 0, 1, 1)
        self.mainTextLayout.addWidget(drr_val, 4, 1, 1, 1)
        self.mainTextLayout.addWidget(df_f, 4, 2, 1, 1)
        self.mainTextLayout.addWidget(dff_val, 4, 3, 1, 1)

        fnu = QLabel("fnu(Hz):")
        fnu_val = QLabel(" ")
        self.mainTextLayout.addWidget(fnu, 5, 0, 1, 1)
        self.mainTextLayout.addWidget(fnu_val, 5, 1, 1, 1)
        dfnu = QLabel("dfnu(Hz):")
        dfnu_val = QLabel(" ")
        self.mainTextLayout.addWidget(dfnu, 5, 2, 1, 1)
        self.mainTextLayout.addWidget(dfnu_val, 5, 3, 1, 1)

        vrfkv = QLabel("Vrf(kV)")
        vrfkv_value = QLabel("1.20602")
        self.mainTextLayout.addWidget(vrfkv, 6, 0, 1, 2)
        self.mainTextLayout.addWidget(vrfkv_value, 6, 2, 1, 2)

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

        rfv_header = QLabel("RF voltage/turn (kV)")
        layout.addWidget(rfv_header, 4, 0, 1, 2)
        rfv_line = QLineEdit("100")
        layout.addWidget(rfv_line, 4, 2, 1, 1)

        rfh_header = QLabel("RF harmonic number")
        layout.addWidget(rfh_header, 5, 0, 1, 2)
        rfh_line = QLineEdit("360")
        layout.addWidget(rfh_line, 5, 2, 1, 1)

        w_minus = QPushButton("--W")
        w_minus.clicked.connect(self.W_minus)
        bucket_ok = QPushButton("Bucket OK")
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
                "Ruthenium (ru)",
                "Zirconium (zr)",
                "others",
            ]
        )
        
        self.species_line.setCurrentIndex(2)
        self.species_line.currentTextChanged.connect(self.updateCharge)
        layout.addWidget(self.species_line, 7, 2, 1, 1)
        charge_header = QLabel("Charge State")
        charge_header.setFont(headerFont)
        layout.addWidget(charge_header, 8, 0, 1, 2)
        self.charge_line = QLineEdit("79")
        layout.addWidget(self.charge_line, 8, 2, 1, 1)

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
        bw_minus.clicked.connect(self.W_minus)
        bunch_ok = QPushButton("Bunch OK")
        bw_plus = QPushButton("++W")
        bw_plus.clicked.connect(self.W_plus)
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

        cp_widget.setLayout(layout)
        return cp_widget

    def RFFreq_action(self):
        text = self.rff_header.currentText()
        value = self.rf_dict[text]
        # calculate something here?
        self.rff_line.setText(str(value))

        # figure out what the actual value should be
        # self.rf_dict[text] = float(self.rff_line.text())
        #for key, value in self.rf_dict.items():
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
        species = self.species_line.currentText()
        self.chargeState = self.charge_dict[species]
        self.charge_line.setText(str(self.chargeState))

    def checkSpeciesCharge(self):
        # IDK WHAT TO CONNECT THIS TO....
        if int(self.species_line.currentText()) < 0:
            print("atomic number is always > 0")
        
    def bltCoor(self, graph, winX, winY):
        pos = {}
        x, y = map(float, graph.invtransform(winX, winY).split())
        pos['x'], pos['y'] = x, y

        #self._pos = 0
        self.current_x_position = 0
        self.current_y_position = 0
        self.dE = (self.current_y_position*2.0*bmath.pi*self.frf_1)
        self.frf_1 = 0
        self.dE_Es = self.dE/self.Es*bmath.kilo*bmath.mega
        self.Es = 0
        self.A_1 = 0
        self.B_1 = 0
        self.phis = 0 
        self.fnu_1 = 1.0/(self.tnu+1.0e-20)
        self.dP = self.dE/self.betas_1
        self.dP_Ps = self.dE_Es/(self.betas_1*self.betas_1)
        self.betas_1 = 0
        self.phase_t1 = self.current_x_position/360/(self.frf_1 + 1.0e-20)*bmath.kilo
        self.Vrf_i = self.A_1*pos[y]**2* bmath.pi * self.h
        self.Vrf_i = abs(self.Vrf_i) / self.e / bmath.kilo / (1.0e-20 + abs(math.cos(x) - math.cos(self.phis) + (x - self.phis) * math.sin(self.phis)))
        self.dfnu = -(self.fnu - self.fnu_1)
        self.fnu = 0 
        self.dR_R = self.dP_Ps / self.gamma_tr**2
        self.df_f = -self.dP_Ps * self.etas_1
        # Para.pi, Para.mega, Para.kilo
        self.etas_1 = 0
        self.gamma_tr = 0

        self.x = pos['x'] * bmath.pi/180
        self.tnu = self.period_bun(self.phis, self.x, self.A_1, self.B_1)

    def bun(self, graph, winX, winY):
        pos = {}
        x, y = map(float, graph.invtransform(winX, winY).split())
        pos['x'], pos['y'] = x, y

        self.phibun1 = 0
        self.phibun2 = pos['x'] / bmath.radeg
        self.BUNe = 0 
        self.BUNld = 0
        self.BUNlr = 0
        self.phi2 = bmath.pi - self.phis # phi2 of bkt
        
        if self.etas_1 >= 0:
            if self.phibun2 <= self.phi2 and self.phibun2 > self.phis:
                self.phibun1 = Phi_1_bun(self.phis, self.phibun2)
                self.BUNlr = abs(self.phibun2-self.phibun1)
                self.phi12 = self.BUNlr
                self.BUNlt = self.BUNlr / (2.0*bmath.pi * self.frf_1 * bmath.mega) / bmath.ns
                self.BUNld = self.phi12 * bmath.radeg
                self.abun = Alpha_bun_phi12(self.phi12, self.phis)
                self.BUNe_1 = self.abun * self.sk_bkt
                self.BUNe = self.BUNe_1
		        #DrawBun $phis $phibun1 $phibun2 $A_1 $B_1
		        #.gbkt element show {BUN BKT}
        else:
            if self.phibun2 >= self.phi2 and self.phibun2 < self.phis:
                self.phibun2 = Generic_phi2(self.phis, self.phibun2)
                self.phibun1 = Phi_1_bun(self.phis, self.phibun2)
                self.BUNlr = abs(self.phibun2-self.phibun1)
                self.phi12 = self.BUNlr
                self.BUNlt = self.BUNlr/(2.0*bmath.pi*self.frf_1*bmath.mega)/bmath.ns #MHz
                self.BUNld = self.phi12*bmath.radeg
                self.abun = Alpha_bun_phi12(self.phi12, self.phis)
                self.BUNe_1 = self.abun*self.st_bkt
                self.BUNe = self.BUNe_1
                self.phibun2 = i_Alpha_bun(self.abun, self.phis)
                self.phibun1 = Phi_1_bun(self.phis, self.phibun2)
                self.phibun1 = Proper_phi(self.phis, self.phibun1)
                self.phibun2 = Proper_phi(self.phis, self.phibun2)
		        #DrawBun $phis $phibun1 $phibun2 $A_1 $B_1
		        #.gbkt element show {BUN BKT}


class cBucket():
    def __init__(self):
        self.xcir = 0
        self.xrho = 0
        self.xtr = 0

    def cMachine(self):
        if self.xcir < 0:
            print("Error: circumference is always > 0")
        if self.xrho < 0:
            print("Error: Bending radius is always > 0")
        if self.xtr < 0:
            print("Error: gamma_tr is always > 0")

    

#class Utilities:
    #def bltvector(self, n):
        #try:
            # assumes n is a global variable
        #    blt.vector(n)
        #except:
        #    pass

    #def environment(self):
        # I'm iffy about this, there's a different way to find the path
    #     env = os.environ
    #    try:
    #        path = env["PATH"]
    #    except KeyError:
    #        return 1
    #    paths = path.split(":")