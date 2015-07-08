
''' Sn properties - Masses, Ratios, Cup Configurations, ... '''
from sn_config import *
from collections import Counter
import pandas as pd

'''Classes for Reading in the Data and applying Internal Normalisation'''
from nu_data_reduction import NU_data_read, int_norm, evaluation

path = "/Volumes/friebelm/PhD/NU Plasma/Measurements/2015-05-03/"

# First and Last Datafile Number of the measurement
#files = range(50, 57)
#files = [31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45]
#files = {101 : [101, [100, 102]], }
#dict with files and blank files, eg {31 : [31, [30, 32]], }
#blank files
#files = range(1122, 1249, 1)
#for i in range(261, 301, 2):
#    files[i] = [i, [i-1, i+1]]
#print files

# cup configuration
cup_config = cycle_Sb
# Isotopes used for Interference correction
corr_isotopes_1 = {"Cd" : "111", "Te" : "125"}
corr_isotopes_2 = {"Cd" : "111","Te" : "125", "Xe" : "129"}
corr_isotopes_Sb = {"Te": "125"}
# Mass Range of cup configuration
mass_range = cycle_Sb_mass_range
#isotopes = [["111", "112", "114", "115", "116", "117", "118", "119", "122", "124", "125", "126", "129"], ["117", "118", "119", "122", "124"]]
#isotopes = [["112", "114", "115", "116", "117", "118", "119", "122", "124"]]
isotopes = [["116", "117", "118", "119", "122", "124"]]

#number of iterations for beta
iter_beta = 10

# Interference_corr on the denominator isotope
isotope_denom_corr = False


#df = NU_data_read(path, files[33][0], cup_config)
#blk_1 = NU_data_read(path, files[33][1][0], cup_config)
#blk_2 =  NU_data_read(path, files[33][1][1], cup_config)

data_sample = {}
sample_names = {}
avg_ratio_sample_all = {}
sd2_ratio_sample_all = {}
avg_ratio_sample_all_2 = {}
sd2_ratio_sample_all_2 = {}

#df_zero = df.data_zero_corr(33)
#df_bgd_1 = blk_1.data_zero_corr(31)
#df_bgd_2 = blk_2.data_zero_corr(32)

#data_sample[sample] = {}
#cycles = range(1, 21)

#new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)
#print new_corr.get_df()
#bgd_corr = new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)
#raw = new_corr.raw_ratios()
#print raw

files_1 = range(1518, 1729, 1)
#files_1.remove(1146)


for sample in files_1:
    df = NU_data_read(path, sample, cup_config)
    cycles = range(1, len(df.data_read(sample).index)+1)
    #blk_1 = NU_data_read(path, files[sample][1][0], cup_config)
    #blk_2 = NU_data_read(path, files[sample][1][1], cup_config)

    df_zero = df.data_zero_corr(sample)
    #df_bgd_1 = blk_1.data_zero_corr(files[sample][1][0])
    #df_bgd_2 = blk_2.data_zero_corr(files[sample][1][1])

    new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes_Sb, denom_corr_ratio)
    new_corr.line2_corr(df_zero, "119")
    #new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)

    data_sample[sample] = new_corr.external_norm_Sb()
    #data_sample[sample] = new_corr.mass_fractionation()

    avg_ratio_sample_all[sample] = new_corr.avg_to_df(data_sample, sample)
    sd2_ratio_sample_all[sample] = new_corr.SD_to_df(data_sample, sample)
    sample_names[sample] = df.extract_metadata(sample, "Sample Name")

data_all_1 = new_corr.to_df_all(sample_names, avg_ratio_sample_all, sd2_ratio_sample_all, ratios = False)

#data_all_1.columns = ["Sample Name", "111Cd (V)", "112SnCd (V)", "113CdIn (V)", "114SnCd (V)", "115SnIn (V)", "116SnCd (V)", "117Sn (V)", "118Sn (V)", "119Sn (V)","120Sn (V)", "122Sn (V)", "124Sn (V)", "125Te (V)", "126TeXe (V)", "129Xe (V)"]



# for sample in files_2:
#     df = NU_data_read(path, sample, cup_config)
#     #blk_1 = NU_data_read(path, files[sample][1][0], cup_config)
#     #blk_2 = NU_data_read(path, files[sample][1][1], cup_config)
#
#     df_zero = df.data_zero_corr(sample)
#     #df_bgd_1 = blk_1.data_zero_corr(files[sample][1][0])
#     #df_bgd_2 = blk_2.data_zero_corr(files[sample][1][1])
#
#     new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes_2, denom_corr_ratio)
#
#     #new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)
#
#     data_sample[sample] = new_corr.internal_norm_1()
#
#
#     avg_ratio_sample_all_2[sample] = new_corr.avg_to_df(data_sample, sample)
#     sd2_ratio_sample_all_2[sample] = new_corr.SD_to_df(data_sample, sample)
#
# data_all_2 = new_corr.to_df_all(avg_ratio_sample_all_2, sd2_ratio_sample_all_2)
#
# data_all = data_all_1.append(data_all_2)
# print data_all
data_all_1.to_csv('/Volumes/friebelm/PhD/NU Plasma/Measurements/2015-05-03/' + "2015_05_03_Sn_H6_L3_Sb_corr" + ".csv")


# Sample-Standard-Bracketing
#print ((data_all.loc[265]/((data_all.loc[263]+data_all.loc[267])/2))-1)*1000

    # signals = ["111", "112", "113", "114", "115", "116", "117", "118", "119", "120", "122", "124", "126"]
    # data_sample[sample] = {}
    # for cycle in cycles:
    #     data_sample[sample][cycle] = {}
    #     corr = int_norm(df_zero, cycle, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)
    #     for signal in signals:
    #         data_sample[sample][cycle][signal] = corr.lookup_signal(signal)
