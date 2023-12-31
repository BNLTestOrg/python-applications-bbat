import math


class bmath:
    machineList = ["ags", "booster", "rhic"]

    # ags constants
    c_ags = 807.12  # circumference. m
    ro_ags = 128.46  # average radius
    rho_ags = 85.37  # magnetic bending radius, m
    gamma_tr_ags = 8.5  # transition gamma
    alpha_p_ags = 0.013841  # compaction factor
    h_ags = 8.0  # harmonic number
    ags_values = {
        "C": c_ags,
        "Ro": ro_ags,
        "rho": rho_ags,
        "gamma_tr": gamma_tr_ags,
        "alpha_p": alpha_p_ags,
        "h": h_ags,
    }

    # booster constants
    c_booster = 201.78  # circumference, 1/4 of the AGS. m
    ro_booster = 32.114  # average radius
    rho_booster = 13.75099  # magnetic bending radius, m
    gamma_tr_booster = 4.88  # transition gamma
    alpha_p_booster = 0.0419914  # compaction factor
    h_booster = 2.0  # Harmonic number
    booster_values = {
        "C": c_booster,
        "Ro": ro_booster,
        "rho": rho_booster,
        "gamma_tr": gamma_tr_booster,
        "alpha_p": alpha_p_booster,
        "h": h_booster,
    }

    # rhic constants
    c_rhic = 3833.852  # circumference, 19/4 of AGS. m
    ro_rhic = 610.176  # average radius
    rho_rhic = 243.241  # magnetic bending radius, m
    gamma_tr_rhic = 22.8  # transition gamma in CD
    alpha_p_rhic = 1.0 / gamma_tr_rhic
    h_rhic = 360  # harmonic number
    rhic_values = {
        "C": c_rhic,
        "Ro": ro_rhic,
        "rho": rho_rhic,
        "gamma_tr": gamma_tr_rhic,
        "alpha_p": alpha_p_rhic,
        "h": h_rhic,
    }

    machineParameters = {
        "AGS": ags_values,
        "Booster": booster_values,
        "RHIC": rhic_values,
    }

    eo_p = 938.28  # MeV/nucleon
    eo_d = 937.81  # MeV/nucleon
    eo_c = 931.25  # MeV/nucleon
    # eo_o = 931.198  # MeV/nucleon
    eo_s = 930.47  # MeV/nucleon
    eo_si = 930.57  # MeV/nucleon
    eo_cu = 930.29  # MeV/nucleon
    eo_i = 930.68  # MeV/nucleon
    eo_au = 931.26  # MeV/nucleon
    eo_fe = 928.95  # MeV/nucleon
    eo_u = 931.51  # MeV/nucleon
    # eo_ru = 930.60  # MeV/nucleon
    # eo_zr = 930.60  # MeV/nucleon
    Eo = {
        "Proton (p)": eo_p,
        "Deuteron (d)": eo_d,
        "Carbon (c)": eo_c,
        "Sulfur (s)": eo_s,
        "Silicon (si)": eo_si,
        "Copper (cu)": eo_cu,
        "Iodine (i)": eo_i,
        "Gold (au)": eo_au,
        "Iron (fe)": eo_fe,
        "Uranium (u)": eo_u,
    }

    A_p = 1  # nucleon number
    A_d = 2
    A_c = 12
    # A_o = 16
    A_s = 32
    A_si = 28
    A_cu = 63
    A_i = 127
    A_au = 197
    A_fe = 56
    A_u = 238
    # A_ru = 96
    # A_zr = 96
    A = {
        "Proton (p)": A_p,
        "Deuteron (d)": A_d,
        "Carbon (c)": A_c,
        "Sulfur (s)": A_s,
        "Silicon (si)": A_si,
        "Copper (cu)": A_cu,
        "Iodine (i)": A_i,
        "Gold (au)": A_au,
        "Iron (fe)": A_fe,
        "Uranium (u)": A_u,
    }

    Q_p = 1  # atomic number
    Q_d = 2
    Q_c = 6
    # Q_o = 8
    Q_s = 16
    Q_si = 14
    Q_cu = 29
    Q_i = 53
    Q_au = 79
    Q_fe = 26
    Q_u = 92
    # Q_ru = 44
    # Q_zr = 40
    Q = {
        "Proton (p)": Q_p,
        "Deuteron (d)": Q_d,
        "Carbon (c)": Q_c,
        "Sulfur (s)": Q_s,
        "Silicon (si)": Q_si,
        "Copper (cu)": Q_cu,
        "Iodine (i)": Q_i,
        "Gold (au)": Q_au,
        "Iron (fe)": Q_fe,
        "Uranium (u)": Q_u,
    }

    c = 299792458.0  # Speed of light in vacuum, m/s (exact)
    e_e = 1.60217733e-19  # Electronic charge, coloumb (49)
    h_bar = 1.05457266e-34  # Planck constant over 2, joule sec (63)
    eV = 1.60217733e-19  # Electron volt, Joule
    KeV = 1.60217733e-16  # Kiloelectron volt, Joule
    MeV = 1.60217733e-13  # Megaelectron volt, Joule
    GeV = 1.60217733e-10  # Gegaelectron volt, Joule
    TeV = 1.60217733e-7  # Gegaelectron volt, Joule
    KV = 1.0e3  # Kilovolts
    kV = 1.0e3  # Kilovolts
    ms = 1.0e-3  # seconds
    us = 1.0e-6  # seconds
    ns = 1.0e-9  # seconds
    u = 1.6605402e-27  # Atomic mass unit, kg, (10)
    tera = 1.0e12
    giga = 1.0e9
    mega = 1.0e6
    kilo = 1.0e3
    Eo_pro = 938.27231  # Proton rest mass, MeV
    Eo_ele = 0.51099906  # Electron rest mass, MeV
    gg = 9.80665  # Standard acceleration of gravity, m sec-2

    pi = 3.14159265358979323846
    pi_2 = 1.57079632679489661923
    pi2 = 6.28318530717958647692
    ROOT2 = 1.41421356237309504880
    epsilon = 1.0e-8

    radeg = 180 / pi

    v1max = 100
    vnmax = 100
    vlomax = -200
    vupmax = 200
    lotheta = -180
    uptheta = 180

    def abkt_s(a, b):
        return 16.0 * math.sqrt(abs(b * 1.0 / a))

    def alpha_bkt(phis):
        return math.cos(phis / 2) ** 2 + (2 / math.pi) * math.sin(phis / 2) ** 2

    def formatXLabels(graph, x):
        return str(int(x)) + "\u00B0"

    def g2XLabels(graph, x):
        return str(int(x)) + "\u00B0"
