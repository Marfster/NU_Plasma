__author__ = 'matthias'

''' Sn properties - Masses, Ratios, Cup Configurations, ... '''
from sn_config import *

'''Classes for Reading in the Data and applying Internal Normalisation'''
from nu_data_reduction import NU_data_read, int_norm, evaluation

''' Input parameter '''
path = "/Volumes/friebelm/PhD/NU Plasma/Measurements/2015-05-03/"
files_1 = range(1518, 1729, 1)

''' Configuration '''
# cup configuration
cup_config = cycle_Sb
# Isotopes used for Interference correction
corr_isotopes_1 = {"Cd" : "111", "Te" : "125"}
corr_isotopes_2 = {"Cd" : "111","Te" : "125", "Xe" : "129"}
corr_isotopes_Sb = {"Te": "125"}
# Mass Range of cup configuration
mass_range = cycle_Sb_mass_range
# to evaluating isotopes
isotopes = [["116", "117", "118", "119", "122", "124"]]

#number of iterations for beta
iter_beta = 10

# Interference_corr on the denominator isotope
isotope_denom_corr = False
denom_corr_ratio = {"isotope_nom" : "117", "isotope_denom" : "119"}
# Isotope used for correcting line2
isotope_line2_corr = "119"
# Isotope used for referencing
isotope_denom = "120"
'''OPTIONS'''
bgd_corr = False
line2_corr = False
eval_method = 1


def data_process(path, files, cup_config, isotopes, mass_range, corr_isotopes, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, option, iter_beta, isotope_denom):
    # Empty dataframes
    data_sample = {}
    sample_names = {}
    avg_ratio_sample_all = {}
    sd2_ratio_sample_all = {}

    if bgd_corr:
        for i in files:
            files[i] = [i, [i-1, i+1]]


    for sample in files:
        df = NU_data_read(path, sample, cup_config)
        cycles = range(1, len(df.data_read(sample).index)+1)

        if bgd_corr == True:
            df_zero = df.data_zero_corr(sample)
            blk_1 = NU_data_read(path, files[sample][1][0], cup_config)
            blk_2 = NU_data_read(path, files[sample][1][1], cup_config)
            df_bgd_1 = blk_1.data_zero_corr(files[sample][1][0])
            df_bgd_2 = blk_2.data_zero_corr(files[sample][1][1])
            new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)
            new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)

        df_zero = df.data_zero_corr(sample)
        new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)

        if line2_corr == True:
            new_corr.line2_corr(df_zero, isotope_line2_corr)

        methods = {1 : new_corr.raw_signals, 2: new_corr.raw_ratios, 3: new_corr.raw_ratios_corr,
                   4: new_corr.mass_fractionation, 5: new_corr.internal_norm_1, 6: new_corr.internal_norm_2,
                   7: new_corr.external_norm_Sb}

        if option < 2:
            data_sample[sample] = methods[option]()

        elif option >= 2 and option < 4:
            data_sample[sample] = methods[option](isotope_denom)

        else:
            data_sample[sample] = methods[option](isotope_denom, iter_beta)

        avg_ratio_sample_all[sample] = new_corr.avg_to_df(data_sample, sample)
        sd2_ratio_sample_all[sample] = new_corr.SD_to_df(data_sample, sample)
        sample_names[sample] = df.extract_metadata(sample, "Sample Name")

        data_all = new_corr.to_df_all(sample_names, avg_ratio_sample_all, sd2_ratio_sample_all, ratios = False)

    return data_all

#data_raw_signals = data_process(path, files_1, cup_config, isotopes, mass_range, corr_isotopes_Sb, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, 1, iter_beta, isotope_denom)

data_raw_ratios = data_process(path, files_1, cup_config, isotopes, mass_range, corr_isotopes_Sb, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, 3, iter_beta, isotope_denom)

#data_external_sb = data_process(path, files_1, cup_config, isotopes, mass_range, corr_isotopes_Sb, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, 7, iter_beta, isotope_denom)

print data_raw_ratios
#print data_external_sb
