class SomethingAboutConstants:
    def __init__(self) -> None:
        self.machineList = ["ags", "booster", "rhic"]

        # ags constants
        self.c_ags = 807.12  # circumference. m
        self.Ro_ags = 128.46  # average radius
        self.rho_ags = 85.37  # magnetic bending radius, m
        self.gamma_tr_ags = 8.5  # transition gamma
        self.alpha_p_ags = 0.013841  # compaction factor
        self.h_ags = 8.0  # harmonic number

        # booster constants
        self.c_booster = 201.78  # circumference, 1/4 of the AGS. m
        self.ro_booster = 32.114  # average radius
        self.rho_booster = 13.75099  # magnetic bending radius, m
        self.gamma_tr_booster = 4.88  # transition gamma
        self.alpha_p_booster = 0.0419914  # compaction factor
        self.h_booster = 2.0  # Harmonic number

        # rhic constants
        self.c_rhic = 3833.852  # circumference, 19/4 of AGS. m
        self.ro_rhic = 610.176  # average radius
        self.rho_rhic = 243.241  # magnetic bending radius, m
        self.gamma_tr_rhic = 22.8  # transition gamma in CD
        self.alpha_p_rhic = 1.0 / self.gamma_tr_rhic
        self.h_rhic = 360  # harmonic number

        self.eo_p = 938.28  # MeV/nucleon
        self.eo_d = 937.81  # MeV/nucleon
        self.eo_c = 931.25  # MeV/nucleon
        self.eo_o = 931.198  # MeV/nucleon
        self.eo_s = 930.47  # MeV/nucleon
        self.eo_si = 930.57  # MeV/nucleon
        self.eo_cu = 930.29  # MeV/nucleon
        self.eo_i = 930.68  # MeV/nucleon
        self.eo_au = 931.26  # MeV/nucleon
        self.eo_fe = 928.95  # MeV/nucleon
        self.eo_u = 931.51  # MeV/nucleon
        self.eo_ru = 930.60  # MeV/nucleon
        self.eo_zr = 930.60  # MeV/nucleon

        self.A_p = 1  # nucleon number
        self.A_d = 2
        self.A_c = 12
        self.A_o = 16
        self.A_s = 32
        self.A_si = 28
        self.A_cu = 63
        self.A_i = 127
        self.A_au = 197
        self.A_fe = 56
        self.A_u = 238
        self.A_ru = 96
        self.A_zr = 96

        self.Q_p = 1  # atomic number
        self.Q_d = 2
        self.Q_c = 6
        self.Q_o = 8
        self.Q_s = 16
        self.Q_si = 14
        self.Q_cu = 29
        self.Q_i = 53
        self.Q_au = 79
        self.Q_fe = 26
        self.Q_u = 92
        self.Q_ru = 44
        self.Q_zr = 40

        self.c = 299792458.0  # Speed of light in vacuum, m/s (exact)
        self.e_e = 1.60217733e-19  # Electronic charge, coloumb (49)
        self.h_bar = 1.05457266e-34  # Planck constant over 2, joule sec (63)
        self.eV = 1.60217733e-19  # Electron volt, Joule
        self.KeV = 1.60217733e-16  # Kiloelectron volt, Joule
        self.MeV = 1.60217733e-13  # Megaelectron volt, Joule
        self.GeV = 1.60217733e-10  # Gegaelectron volt, Joule
        self.TeV = 1.60217733e-7  # Gegaelectron volt, Joule
        self.KV = 1.0e3  # Kilovolts
        self.kV = 1.0e3  # Kilovolts
        self.ms = 1.0e-3  # seconds
        self.us = 1.0e-6  # seconds
        self.ns = 1.0e-9  # seconds
        self.u = 1.6605402e-27  # Atomic mass unit, kg, (10)
        self.tera = 1.0e12
        self.giga = 1.0e9
        self.mega = 1.0e6
        self.kilo = 1.0e3
        self.Eo_pro = 938.27231  # Proton rest mass, MeV
        self.Eo_ele = 0.51099906  # Electron rest mass, MeV
        self.gg = 9.80665  # Standard acceleration of gravity, m sec-2

        self.pi = 3.14159265358979323846
        self.pi_2 = 1.57079632679489661923
        self.pi2 = 6.28318530717958647692
        self.ROOT2 = 1.41421356237309504880
        self.epsilon = 1.0e-8

        self.radeg = 180 / self.pi
