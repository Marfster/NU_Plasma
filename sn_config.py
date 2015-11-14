__author__ = 'matthias'

from iso_properties import Isotope_Masses, Isotope_Ratios, Isotopes_Mass_Range

# Masses of Isotopes

#*** Sn after IUPAC - De Laeter (2003)***#
Sn_masses = Isotope_Masses()
Sn_masses.add_masses_dict({"112" : 111.904823, "114" : 113.902783, "115" : 114.903347,
                           "116" : 115.901745, "117" : 116.902955, "118" : 117.901608,
                           "119" : 118.903311, "120" : 119.9021985, "122" : 121.9034411,
                           "124" : 123.9052745})

#*** Cd after IUPAC - De Laeter (2003)***#
Cd_masses = Isotope_Masses()
Cd_masses.add_masses_dict({"110" : 105.906458, "111" : 110.904182, "112" : 111.9027577,
                           "113" : 112.9044014, "114" : 113.9033586, "116" : 115.904756})

In_masses = Isotope_Masses()
In_masses.add_masses_dict({"113" : 112.9040624, "115" : 114.90387940})

#*** Sb after IUPAC - De Laeter (2003)***#
Sb_masses = Isotope_Masses()
Sb_masses.add_masses_dict({"121" : 120.903822226, "123" : 122.904216022})

#*** Te ***#
Te_masses = Isotope_Masses()
Te_masses.add_masses_dict({"120" : 119.904021, "122" : 121.903055, "124" : 123.902825,
                           "125" : 124.904435, "126" : 125.903310, "128" : 127.904461519,
                           "130": 129.906222921})

#*** I after IUPAC - De Laeter (2003)***#
I_masses = Isotope_Masses()
I_masses.add_masses_dict({"127" : 126.9044684})

#*** Xe after IUPAC - De Laeter (2003)***#
Xe_masses = Isotope_Masses()
Xe_masses.add_masses_dict({"124" : 123.905895421, "126" : 125.9042687, "129" : 128.90477999,
                           "130" : 129.903508911, "131" : 130.905082818, "132" : 131.904154615,
                           "134" : 133.90539459, "136": 135.9072208})

#*** Pd after IUPAC - De Laeter (2003)***#
Pd_masses = Isotope_Masses()
Pd_masses.add_masses_dict({"110" : 109.90515312})


# True Ratios

#*** Cd after Loss et al. (1990)***#
Cd_ratios = Isotope_Ratios()
Cd_ratios.add_ratios_dict("114", { "106": 0.046118, "108" : 0.032379, "110": 0.44766 , "111" : 0.45550,
                                  "112" : 0.85193, "113" : 0.42842, "116" : 0.25679})
#*** In after IUPAC (2001)***#
In_ratios = Isotope_Ratios()
In_ratios.add_ratios_dict("113", {"115" : 0.95715/0.04295})

#*** Sn after Lee et al. (1995)***#
Sn_ratios = Isotope_Ratios()
Sn_ratios.add_ratios_dict("120", {"112" : 0.029812, "114" : 0.020195, "115" : 0.010366, "116" : 0.4460,
                                  "117" : 0.235313, "118" : 0.742935, "119" : 0.26343,
                                  "122" : 0.142086, "124" : 0.177588})

#*** Sb***#
Sb_ratios = Isotope_Ratios()
Sb_ratios.add_ratios_dict("121", {"123": 0.7479})

#*** Te after Lee et al. (1995)***#
Te_ratios = Isotope_Ratios()
Te_ratios.add_ratios_dict("128", {"120" : 0.002919, "122" : 0.079603, "123" : 0.027904, "124" : 0.14853,
                                  "125" : 0.222041, "126" : 0.592264, "130" : 1.075950})

#*** Xe after Basford (1973)***#
Xe_ratios = Isotope_Ratios()
Xe_ratios.add_ratios_dict("130", {"124" : 0.02337, "126" : 0.02179, "128" : 0.47150, "129" : 6.49631,
                                  "131" : 5.21376, "132" : 6.60688, "134" : 2.56265, "136" : 2.17617})

Sn_ratios.transf_ratio("120", "119")
Sn_ratios.transf_ratio("120", "118")
Te_ratios.transf_ratio("128", "126")
Te_ratios.transf_ratio("128", "125")
Xe_ratios.transf_ratio("130", "129")
Cd_ratios.transf_ratio("114", "111")

database = {"Sn" : {"Masses" : Sn_masses, "Ratios" : Sn_ratios}, "Cd" : {"Masses" : Cd_masses, "Ratios" : Cd_ratios},
            "In" : {"Masses" : In_masses, "Ratios" : In_ratios}, "Te" : {"Masses" : Te_masses, "Ratios" : Te_ratios},
            "Xe" : {"Masses" : Xe_masses, "Ratios" : Xe_ratios}, "Sb" : {"Masses" : Sb_masses, "Ratios" : Sb_ratios},
            "I" : {"Masses" : I_masses}, "Pd" : {"Masses" : Pd_masses}}

# Isotopes on Mass range 110 to 131
cycles1_mass_range = Isotopes_Mass_Range()
cycles1_mass_range.add_mass_range_dict({"110" : ["Pd", "Cd"], "111" : ["Cd"], "112" : ["Sn", "Cd"],
                                        "113" : ["Cd", "In"], "114" : ["Sn", "Cd"], "115" : ["Sn", "In"],
                                        "116" : ["Sn", "Cd"], "117" : ["Sn"], "118" : ["Sn"], "119" : ["Sn"],
                                        "120" : ["Sn", "Te"], "122" : ["Sn", "Te"], "124" : ["Sn", "Te", "Xe"],
                                        "126" : ["Te", "Xe"]})

cycles2_mass_range = Isotopes_Mass_Range()
cycles2_mass_range.add_mass_range_dict({"110" : ["Pd", "Cd"], "111" : ["Cd"], "112" : ["Sn", "Cd"],
                                        "113" : ["Cd", "In"], "114" : ["Sn", "Cd"], "115" : ["Sn", "In"],
                                        "116" : ["Sn", "Cd"], "117" : ["Sn"], "118" : ["Sn"], "119" : ["Sn"],
                                        "120" : ["Sn", "Te"], "121" : ["Sb"], "122" : ["Sn", "Te"],
                                        "123" : ["Sb", "Te"], "124" : ["Sn", "Te", "Xe"], "125" : ["Te"],
                                        "126" : ["Te", "Xe"], "127" : ["I"], "129" : ["Xe"], "131" : ["Xe"]})

cycle_Sb_mass_range = Isotopes_Mass_Range()
cycle_Sb_mass_range.add_mass_range_dict({"116" : ["Sn", "Cd"], "117" : ["Sn"], "118" : ["Sn"], "119" : ["Sn"],
                                         "120" : ["Sn", "Te"], "121" : ["Sb"], "122" : ["Sn", "Te"],
                                         "123" : ["Sb", "Te"], "124" : ["Sn", "Te", "Xe"], "125" : ["Te"]})

cycle_Sb_2_mass_range = Isotopes_Mass_Range()
cycle_Sb_2_mass_range.add_mass_range_dict({"116" : ["Sn", "Cd"], "117" : ["Sn"], "118" : ["Sn"], "119" : ["Sn"],
                                         "120" : ["Sn", "Te"], "121" : ["Sb"], "122" : ["Sn", "Te"],
                                         "123" : ["Sb", "Te"], "124" : ["Sn", "Te", "Xe"], "125" : ["Te"],
                                          "127" : ["I"], "129" : ["Xe"]})


# Isotopes used for Interference correction
#cycles1_corr = {}
cycles1_corr = {"Cd" : "111", "Te" : "126"}
#cycles2_corr = {}
cycles2_corr = {"Cd" : "111","Te" : "125", "Xe" : "129"}
cycle_Sb_corr = {"Te" : "125"}
# Isotopes used for Interference correction on denominator
denom_corr_ratio = {"isotope_nom" : "117", "isotope_denom" : "119"}

# cup configuration

cycles1 = {"cycle1" : {"H9 (1)" : "126", "H8 (1)" : "124", "H7 (1)" : "122", "H6 (1)" : "120", "H5 (1)" : "119",
                       "H4 (1)" : "118", "H3 (1)" : "117", "H2 (1)" : "116", "H1 (1)" : "115", "Ax (1)" : "114",
                       "L1 (1)" : "113", "L2 (1)" : "112", "L3 (1)" : "111", "L4 (1)" : "110"},
           "zero1" : {"H9 (Z1)" : "126", "H8 (Z1)" : "124", "H7 (Z1)" : "122", "H6 (Z1)" : "120", "H5 (Z1)" : "119",
                       "H4 (Z1)" : "118", "H3 (Z1)" : "117", "H2 (Z1)" : "116", "H1 (Z1)" : "115", "Ax (Z1)" : "114",
                       "L1 (Z1)" : "113", "L2 (Z1)" : "112", "L3 (Z1)" : "111", "L4 (Z1)" : "110"}}

cycles2 = {"cycle1" : {"H8 (1)" : "124", "H7 (1)" : "122", "H6 (1)" : "120", "H5 (1)" : "119",
                       "H4 (1)" : "118", "H3 (1)" : "117", "H2 (1)" : "116", "H1 (1)" : "115", "Ax (1)" : "114",
                       "L1 (1)" : "113", "L2 (1)" : "112", "L3 (1)" : "111", "L4 (1)" : "110"},
           "cycle2" : {"H8 (2)" : "131", "H7 (2)" : "129", "H6 (2)" : "127", "H5 (2)" : "126",
                       "H4 (2)" : "125", "H3 (2)" : "124", "H2 (2)" : "123", "H1 (2)" : "122", "Ax (2)" : "121",
                       "L1 (2)" : "120", "L2 (2)" : "119", "L3 (2)" : "118", "L4 (2)" : "117"},
           "zero1" : {"H8 (Z1)" : "124", "H7 (Z1)" : "122", "H6 (Z1)" : "120", "H5 (Z1)" : "119",
                       "H4 (Z1)" : "118", "H3 (Z1)" : "117", "H2 (Z1)" : "116", "H1 (Z1)" : "115", "Ax (Z1)" : "114",
                       "L1 (Z1)" : "113", "L2 (Z1)" : "112", "L3 (Z1)" : "111", "L4 (Z1)" : "110"}}

cycle_Sb = {"cycle1" : {"H6 (1)" : "125", "H5 (1)" : "124",
                       "H4 (1)" : "123", "H3 (1)" : "122", "H2 (1)" : "121", "H1 (1)" : "120", "Ax (1)" : "119",
                       "L1 (1)" : "118", "L2 (1)" : "117", "L3 (1)" : "116"},
             "zero1" : {"H6 (Z1)" : "125", "H5 (Z1)" : "124",
                       "H4 (Z1)" : "123", "H3 (Z1)" : "122", "H2 (Z1)" : "121", "H1 (Z1)" : "120", "Ax (Z1)" : "119",
                       "L1 (Z1)" : "118", "L2 (Z1)" : "117", "L3 (Z1)" : "116"}}

cycle_Sb_2 = {"cycle1" : {"H8 (1)" : "129", "H7 (1)" : "127", "H6 (1)" : "125", "H5 (1)" : "124",
                       "H4 (1)" : "123", "H3 (1)" : "122", "H2 (1)" : "121", "H1 (1)" : "120", "Ax (1)" : "119",
                       "L1 (1)" : "118", "L2 (1)" : "117", "L3 (1)" : "116"},
             "zero1" : {"H8 (Z1)" : "129", "H7 (Z1)" : "127", "H6 (Z1)" : "125", "H5 (Z1)" : "124",
                       "H4 (Z1)" : "123", "H3 (Z1)" : "122", "H2 (Z1)" : "121", "H1 (Z1)" : "120", "Ax (Z1)" : "119",
                       "L1 (Z1)" : "118", "L2 (Z1)" : "117", "L3 (Z1)" : "116"}}