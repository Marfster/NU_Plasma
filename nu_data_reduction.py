__author__ = 'matthias'

import pandas as pd
import re
from collections import Counter
from math import log
import numpy as np
import math



class NU_data_read(object):

    def __init__(self, path, datafile_list, cup_configuration):

        '''******************************************************************************
           Values for Reading in the data and storing it accordingly to the cup settings
           ****************************************************************************** '''
        self.path = path #Path of Datafiles
        self.files = datafile_list  # List of datafile numbers
        #self.cycles = no_cycles # isotope measurement lines
        self.cups = cup_configuration # "dictionary with cup configuration
        self.zero_cycles = [k for k in cup_configuration.keys() if 'zero' in k] # cycles for one measurements
    '''***********************************
       Functions used for Reading the data
       ***********************************'''

    # Extracts sample name, date, ... of csv file
    def extract_metadata(self, filenumber, metadata_string):
        datafile = "Data_" + str(filenumber) + ".csv"

        df = pd.read_table(self.path + datafile, dtype=str, header=None)
        for row in range(len(df.index)):
            line = str(df[0][row])
            if metadata_string in line:
                metadata = re.search(r""+metadata_string+",(.*)", line).group(1)
        return metadata

    # Reads in the data
    def data_read(self, filenumber):
        datafile = "Data_" + str(filenumber) + ".csv"


        data = pd.read_csv(self.path + datafile, skiprows=58, index_col='Cycle') # before line 58 only headers
        return data

    # Return the measured signals in the cycle of each sample
    def data_signals(self, filenumber):
        data = self.data_read(filenumber)
        columns = data.columns.values.tolist()
        dict_signals = {}
        for cycle in self.cups: # reads through cup configuration and writes signals per cup into a dictionary
            dict_signals[cycle] = {}
            for cup in self.cups[cycle]:
                if (cup in columns):
                    dict_signals[cycle][cup] = data[cup].to_dict()
        return dict_signals

    # Applies the Zero correction of each sample
    def data_zero_corr(self, filenumber):
        data = self.data_signals(filenumber)
        # Outlier correction 2SD
        #cyclex_signal = cyclex_signal[np.abs(cyclex_signal-cyclex_signal.mean())<=(2*cyclex_signal.std())]
        # Zero correction
        zero_signals = {}
        counter = 0
        for cycle in self.cups:
            if (cycle not in self.zero_cycles):
                zero_signals[cycle] = {}
                if (len(self.zero_cycles) == 1):
                    for cup in self.cups[cycle]:
                        match = [k for k in self.cups[self.zero_cycles[counter]].keys() if cup[:-4] in k] #screen throug zero1 dict and match cycle cups
                        if (len(match) > 0):
                            cup_cycle = Counter(data[cycle][cup])
                            cup_zero = Counter(data[self.zero_cycles[counter]][match[0]])
                            cup_cycle.subtract(cup_zero)
                            zero_signals[cycle][cup] = cup_cycle
                else:
                    counter += 1
                    for cup in self.cups[cycle]:
                        match = [k for k in self.cups[self.zero_cycles[counter]].keys() if cup[:-4] in k] #screen throug zero1 dict and match cycle cups
                        if (len(match) > 0):
                            cup_cycle = Counter(data[cycle][cup])
                            cup_zero = Counter(data[self.zero_cycles[counter]][match[0]])
                            cup_cycle.subtract(cup_zero)
                            zero_signals[cycle][cup] = cup_cycle
        return zero_signals


    # Background correction for one sample with one background
    def data_bgd_corr(self, filenumber, filenumber_bgd):
        zero_data = self.data_zero_corr(filenumber)
        zero_data_bgd = self.data_zero_corr(filenumber_bgd)
        bgd_signals = {}
        for cycle in zero_data:
            bgd_signals[cycle] = {}
            for cup in zero_data[cycle]:
                if (cup in zero_data_bgd[cycle]):
                    cup_cycle = Counter(zero_data[cycle][cup])
                    cup_cycle_bgd = Counter(zero_data_bgd[cycle][cup])
                    cup_cycle.subtract(cup_cycle_bgd)
                    bgd_signals[cycle][cup] = cup_cycle
        return bgd_signals

class Neptune_data_read(object):

    def __init__(self,path, datafile_list, cup_configuration):

        '''******************************************************************************
        Values for Reading in the data and storing it accordingly to the cup settings
        ****************************************************************************** '''
        self.path = path  # Path of Datafiles
        self.files = datafile_list  # List of datafile numbers
        #self.cycles = no_cycles # isotope measurement lines
        self.cups = cup_configuration # "dictionary with cup configuration
        '''***********************************
        Functions used for Reading the data
        ***********************************'''

    # Extracts sample name
    def extract_metadata(self, filex, metadata_string):

        df = pd.read_table(filex, dtype=str, header=None)
        for row in range(len(df.index)):
            line = str(df[0][row])
            if metadata_string in line:
                metadata = re.search(r"" + metadata_string + ":(.*)", line).group(1)
        return metadata

    # Reads in the data
    def data_read(self, filex):
        data = pd.read_table(filex, skiprows=16, index_col = "Cycle")
        data = data.ix[:,:"124Sn"]
        data = data.drop(["***"],axis=0)
        data.ix[:,"117Sn":] = data.ix[:,"117Sn":].astype(np.float64)
        data = data.dropna(axis=0)
        data.index = data.index.astype(np.int64)
        return data

    # Return the measured signals in the cycle of each sample
    def data_signals(self, filex):
        data = self.data_read(filex)
        columns = data.columns.values.tolist()
        dict_signals = {}
        for cycle in self.cups:
            dict_signals[cycle] = {}
            for cup in self.cups[cycle]:
                if (cup in columns):
                    dict_signals[cycle][cup] = data[cup].to_dict()
        return dict_signals


    # Background correction for one sample with one background
    def data_bgd_corr(self, filenumber, filenumber_bgd):
        zero_data = self.data_zero_corr(filenumber)
        zero_data_bgd = self.data_zero_corr(filenumber_bgd)
        bgd_signals = {}
        for cycle in zero_data:
            bgd_signals[cycle] = {}
            for cup in zero_data[cycle]:
                if (cup in zero_data_bgd[cycle]):
                    cup_cycle = Counter(zero_data[cycle][cup])
                    cup_cycle_bgd = Counter(zero_data_bgd[cycle][cup])
                    cup_cycle.subtract(cup_cycle_bgd)
                    bgd_signals[cycle][cup] = cup_cycle
        return bgd_signals

class int_norm(object):

    def __init__(self, data_dict, cycle_no, cup_configuration, database, mass_range, isotopes_for_corr, denom_corr_ratio=None):

        self.data_dict = data_dict                                                      #dictionary with data (usually zero and/or bgd corrected)
        self.cycle_no = cycle_no                                                        # number of cycles per measurement
        self.cups = cup_configuration                                                   # cup configuration used for measurement (see "sn_config.py")
        self.zero_cycles = [k for k in self.cups.keys() if 'zero' in k]                 # number of cycles per zero measurement
        self.database = database                                                        # refers to database in "sn_config.py", which contains the used Isotope masses and ratios
        self.mass_range = mass_range                                                    # refers to defined mass range of isotopes for calculation (see "sn_config.py")
        self.isotopes_for_corr = isotopes_for_corr                                      # defines Isotopes used for Interferences correction (see "sn_config.py or ipython notebook")
        self.graph_of_corr = self.mass_range.get_graph_of_corr(self.isotopes_for_corr)  # creates a graph of dependencies for topologic sorting
        self.order_of_corr = self.mass_range.get_order_of_corr(self.isotopes_for_corr)  # creates a order for topologic sorting
        self.denom_corr_ratio = denom_corr_ratio                                        # alternative (interference-free) Isotope used for Interference correction on denominator isotope
    '''***********************************
       Functions used for internal normalisation
       ***********************************'''

    # Lookup Isotope Signal
    def lookup_signal(self, isotope, isotope_from_line1 = True):
        invers_dict = {}
        cup_list = []
        for line in self.data_dict:
            invers_dict[line] = dict(zip(self.cups[line].values(), self.cups[line]))
            if (isotope in invers_dict[line]):
                cup_list.append((line, invers_dict[line][isotope]))
        if isotope_from_line1 == False and len(cup_list) > 1:
            return self.data_dict[cup_list[1][0]][cup_list[1][1]][self.cycle_no]
        else:
            return self.data_dict[cup_list[0][0]][cup_list[0][1]][self.cycle_no]

    # Raw Isotope Ratios

    def raw_ratio(self, isotope_nom, isotope_denom, isotope_from_line1 = True):
        return self.lookup_signal(isotope_nom, isotope_from_line1)/self.lookup_signal(isotope_denom, isotope_from_line1)

    # Mass fractionation correction
    def int_norm(self, element_nom,  element_denom, isotope_nom, isotope_denom, beta, isotope_from_line1 = True, corr_isotope_denom = None):
        signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
        if corr_isotope_denom == None:
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1 = True)
        else:
            signal_isotope_denom = corr_isotope_denom
        mass_isotope_nom = self.database[element_nom]["Masses"].get_Isotope_mass(isotope_nom)
        mass_isotope_denom = self.database[element_denom]["Masses"].get_Isotope_mass(isotope_denom)

        return (signal_isotope_nom/signal_isotope_denom) * (mass_isotope_nom/mass_isotope_denom)**beta

    # Interference correction for one interference on the nominator isotope
    def interference_corr_all(self, element_denom, isotope_denom, beta, corr_isotope_denom = None, isotope_from_line1 = True):
        corr_dict = {}
        counter = 0
        for level in self.order_of_corr:
            if counter == 0:
                for isotope in level:
                    element_nom = self.mass_range.get_isotopes(isotope)[0]
                    corr_dict[isotope] = self.int_norm(element_nom, element_denom, isotope, isotope_denom, beta,
                                                       isotope_from_line1=isotope_from_line1,
                                                       corr_isotope_denom=corr_isotope_denom)
                    #print isotope
                    #print element_nom
                    #print element_denom
                    #print isotope_denom
                    #print corr_dict[isotope]
                counter += 1
            elif counter > 0:
                for isotope in level:
                    corr_isotopes = list(self.graph_of_corr[isotope])

                    #print isotope

                    for corr_isotope in corr_isotopes:
                        #print self.mass_range.get_isotopes(isotope)
                        element_nom = self.mass_range.get_isotopes(isotope)
                        # corr_isotope_set = set(corr_isotopes)
                        # element_nom = element_nom.difference(self.mass_range.get_isotopes(corr_isotope))
                        element_nom = list(element_nom)[0]
                        element_corr = set(self.mass_range.get_isotopes(corr_isotope))
                        element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope))
                        element_corr = list(element_corr)

                        if isotope in corr_dict:
                            None
                        else:
                            corr_dict[isotope] = self.int_norm(element_nom, element_denom, isotope, isotope_denom, beta,
                                                               isotope_from_line1=isotope_from_line1,
                                                               corr_isotope_denom=corr_isotope_denom)

                        true_ratio_corr = self.database[element_corr[0]]["Ratios"].get_ratio(isotope, corr_isotope)
                        corr_dict[isotope] -= corr_dict[corr_isotope] * true_ratio_corr

                        #print isotope
                        #print corr_isotope
                        #print element_corr
                        #print corr_dict[corr_isotope]
                        #print corr_dict[corr_isotope] * true_ratio_corr
                        #print corr_dict[isotope]

                counter += 1

        return corr_dict

    # Interference correction
    def interference_corr_ratio(self, element, isotope_nom, isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True):

        if (isotope_nom in self.order_of_corr[0] == True) and (isotope_denom in self.order_of_corr[0] == True):
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False and (isotope_denom_corr == True):
            new_isotope_denom = self.denom_corr_ratio['isotope_denom']
            isotope_ratio = self.interference_corr_all(element, new_isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_denom] # eg(Sn120/Sn119)
            mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)
            mass_new_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_denom)
            isotope_ratio_raw = isotope_ratio * (mass_isotope_denom/mass_new_isotope_denom)**(-beta)
            isotope_signal = isotope_ratio_raw * self.lookup_signal(new_isotope_denom, isotope_from_line1 = isotope_from_line1)
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_signal, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False:
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        else:
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        return isotope_ratio

    def beta(self, iterations, element, isotope_nom, isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True):

        if not list(self.graph_of_corr[isotope_nom]):
            signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1)
            true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom)
            mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom)
            mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)
            beta_temp = log(true_ratio / (signal_isotope_nom/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom)

        else:

            if isotope_denom_corr == False:
                signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
                signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1)
                true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom)
                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom)
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)

                corr_isotopes = list(self.graph_of_corr[isotope_nom])
                corr_isotope = corr_isotopes[0]
                element_corr = set(self.mass_range.get_isotopes(corr_isotope))
                element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope_nom))
                element_corr = list(element_corr)[0]
                true_ratio_corr = self.database[element_corr]["Ratios"].get_ratio(isotope_nom, corr_isotope)
                signal_corr_isotope = self.lookup_signal(corr_isotope, isotope_from_line1)

                isotope_nom_guess = signal_isotope_nom - (signal_corr_isotope * true_ratio_corr)

                beta_temp = log(true_ratio / (isotope_nom_guess/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom)

                for i in range(iterations):
                    isotope_ratio = self.interference_corr_all(element, isotope_denom, beta_temp)[isotope_nom]
                    if isotope_ratio > 0:
                        isotope_ratio_raw = isotope_ratio * (mass_isotope_nom/mass_isotope_denom)**(-beta_temp)
                        beta_update = log(true_ratio / isotope_ratio_raw) / log(mass_isotope_nom/mass_isotope_denom)

                        convergence = (beta_temp - beta_update) / beta_temp
                        beta_temp = beta_update
                    else:
                        print "Iterative Beta Correction Failed! --> beta_temp before iteration used"
                        beta_temp = log(true_ratio / (isotope_nom_guess/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom)
                        break

            elif isotope_denom_corr == True:
                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) #e.g. M(116Sn)
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)# e.g. M(120Sn)
                true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom)

                new_isotope_nom = self.denom_corr_ratio['isotope_nom'] #e.g. 117Sn
                signal_new_isotope_nom = self.lookup_signal(new_isotope_nom, isotope_from_line1)
                mass_new_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_nom)

                new_isotope_denom = self.denom_corr_ratio['isotope_denom'] #e.g. 119Sn
                signal_new_isotope_denom = self.lookup_signal(new_isotope_denom, isotope_from_line1)
                mass_new_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_denom)

                new_true_ratio = self.database[element]["Ratios"].get_ratio(new_isotope_nom, new_isotope_denom)
                beta_temp = log(new_true_ratio / (signal_new_isotope_nom / signal_new_isotope_denom)) / log(mass_new_isotope_nom / mass_new_isotope_denom)

                for i in range(iterations):
                    new_isotope_ratio_nom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_nom] #e.g. 116Sn/119Sn
                    new_isotope_ratio_denom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_denom] # e.g. 120Sn/119Sn

                    new_isotope_ratio_nom_raw = new_isotope_ratio_nom * (mass_isotope_nom/mass_new_isotope_denom)**(-beta_temp)
                    new_isotope_ratio_denom_raw = new_isotope_ratio_denom * (mass_isotope_denom/mass_new_isotope_denom)**(-beta_temp)

                    isotope_ratio_raw = new_isotope_ratio_nom_raw / new_isotope_ratio_denom_raw

                    beta_update = log(true_ratio / isotope_ratio_raw) / log(mass_isotope_nom/mass_isotope_denom)

                    convergence = (beta_temp - beta_update) / beta_temp
                    beta_temp = beta_update

        return beta_temp


    def interference_corr_all_raw(self, element_denom, isotope_denom, corr_isotope_denom = None, isotope_from_line1 = True):
        corr_dict = {}
        counter = 0
        for level in self.order_of_corr:
            if counter == 0:
                for isotope in level:
                    corr_dict[isotope] = self.raw_ratio(isotope, isotope_denom, isotope_from_line1 = isotope_from_line1)
                counter += 1
            elif counter > 0:
                for isotope in level:
                    corr_isotopes = list(self.graph_of_corr[isotope])
                    if (len(corr_isotopes) > 1):
                        for mass_no in corr_isotopes:
                            if mass_no in self.order_of_corr[counter-1]:
                                corr_isotope = mass_no
                    else:
                        corr_isotope = corr_isotopes[0]

                    element_corr = set(self.mass_range.get_isotopes(corr_isotope))
                    element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope))
                    element_corr = list(element_corr)[0]

                    true_ratio_corr = self.database[element_corr]["Ratios"].get_ratio(isotope, corr_isotope)
                    corr_dict[isotope] = self.raw_ratio(isotope, isotope_denom, isotope_from_line1 = isotope_from_line1)
                    corr_dict[isotope] -= corr_dict[corr_isotope] * true_ratio_corr
                counter += 1

        return corr_dict

    # Interference correction
    def interference_corr_ratio_raw(self, element, isotope_nom, isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True):

        if (isotope_nom in self.order_of_corr[0] == True) and (isotope_denom in self.order_of_corr[0] == True):
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False and (isotope_denom_corr == True):
            new_isotope_denom = self.denom_corr_ratio['isotope_denom']
            isotope_ratio = self.interference_corr_all_raw(element, new_isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_denom] # eg(Sn120/Sn119)
            isotope_signal = isotope_ratio * self.lookup_signal(new_isotope_denom, isotope_from_line1 = isotope_from_line1)
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_signal, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False:
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        else:
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        return isotope_ratio

class evaluation(object):

    def __init__(self, data_dict, cycles, isotopes ,cup_configuration, database, mass_range, isotopes_for_corr, denom_corr_ratio=None):
        self.data_dict = data_dict
        self.cycles = cycles
        self.isotopes = isotopes
        self.cups = cup_configuration
        self.database = database
        self.mass_range = mass_range
        self.isotopes_for_corr = isotopes_for_corr
        self.denom_corr_ratio = denom_corr_ratio

    def get_df(self):
        return self.data_dict

    def line2_corr(self, df, isotope):
        invers_dict = {}
        if "cycle2" in df:
            for cycle in df:
                invers_dict[cycle] = dict(zip(self.cups[cycle].values(), self.cups[cycle]))
                if cycle == "cycle2":
                    keys = df[cycle].keys()
                    keys.append(keys.pop(keys.index(invers_dict["cycle2"][isotope])))
                    for cup in keys:
                        for meas_point in df[cycle][cup]:
                                if (df["cycle1"][invers_dict["cycle1"][isotope]][meas_point]) and  df["cycle2"][invers_dict["cycle2"][isotope]][meas_point]:
                                    line2_corr = df["cycle1"][invers_dict["cycle1"][isotope]][meas_point]/df["cycle2"][invers_dict["cycle2"][isotope]][meas_point]
                                    df[cycle][cup][meas_point] = df[cycle][cup][meas_point] * line2_corr
                                    self.data_dict = df
        return self.data_dict

    #background correction
    def data_bgd_corr(self, df_bgd_1, df_bgd_2):
        if df_bgd_2:
            bgd_signals = {}
            for cycle in df_bgd_1:
                bgd_signals[cycle] = {}
                for cup in df_bgd_1[cycle]:
                    if (cup in df_bgd_2[cycle]):
                        bgd_signals[cycle][cup] = {}
                        for meas_point in df_bgd_1[cycle][cup]:
                            if (meas_point in df_bgd_2[cycle][cup]):
                                avg_cup_cycle_bgd = np.divide((np.add(df_bgd_1[cycle][cup][meas_point], df_bgd_2[cycle][cup][meas_point])),2)
                                bgd_signals[cycle][cup][meas_point] = avg_cup_cycle_bgd
        else:
            bgd_signals = df_bgd_1

        df_bgd_corr = {}
        for cycle in self.data_dict:
            df_bgd_corr[cycle] = {}
            for cup in self.data_dict[cycle]:
                if (cup in bgd_signals[cycle]):
                        names = ['id','data']
                        formats = ['float','float']
                        dtype = dict(names = names, formats=formats)
                        cup_cycle = np.array(self.data_dict[cycle][cup].items(), dtype=dtype)
                        cup_cycle_bgd = np.array(bgd_signals[cycle][cup].items(), dtype=dtype)
                        mean = np.nanmean(cup_cycle_bgd["data"])
                        x2 = np.full(len(self.data_dict[cycle][cup]),(mean))
                        cup_cycle = np.array(cup_cycle)
                        cup_cycle["data"] = cup_cycle["data"] - x2
                        df_bgd_corr[cycle][cup] = Counter(dict(enumerate(cup_cycle["data"], 1)))
        self.data_dict = df_bgd_corr
        return self.data_dict

    def data_bgd_corr_2(self, ):
        # method for blank correction
            # arbitrary blank positions
            #blk1 = [item for item in blk_ls if item < sample]
            #blk2 = [item for item in blk_ls if item > sample]
            #blk1 = blk1[-1]
            #blk2 = blk2[0]
            #blk_corr_sample == True

            # blank measurement surrounds sample
            #if ((sample-1) in blk_ls) and ((sample+1) in blk_ls):
            #    blk1 = sample-1
            #    blk2 = sample+1
            #    blk_corr_sample = True
            #else:
            #    blk_corr_sample = False
            if (blk_corr == True) and (blk_corr_sample == True):
                blk_1 = NU_data_read(path, blk1, cup_config)
                blk_2 = NU_data_read(path, blk2, cup_config)
                df_zero = df.data_zero_corr(sample)
                df_bgd_1 = blk_1.data_zero_corr(blk1)
                df_bgd_2 = blk_2.data_zero_corr(blk2)
                new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes_Sb , denom_corr_ratio)
                new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)
            else:
                df_zero = df.data_zero_corr(sample)
                new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes_Sb , denom_corr_ratio)

    #raw signals
    def raw_signals(self, isotope_denom):
        data_sample = {}
        isotopes_0 = self.isotopes[0]
        isotopes_0.append(isotope_denom)

        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            data_sample[n] = {}

            for isotope in isotopes_0:
                data_sample[n][isotope] = corr.lookup_signal(isotope, isotope_from_line1=True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.lookup_signal(isotope, isotope_from_line1=False)

        return data_sample
    #raw ratios
    def raw_ratios(self, isotope_denom):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            data_sample[n] = {}

            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.raw_ratio(isotope, isotope_denom)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.raw_ratio(isotope, isotope_denom)
        return data_sample

    def raw_ratios_corr(self, isotope_denom):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio_raw("Sn", isotope , isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio_raw("Sn", isotope , isotope_denom, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # internal normalisation

    # 1 only mass fractionation
    def mass_fractionation(self, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            beta = corr.beta(iter, "Sn", "116", "120", isotope_denom_corr = False)
            data_sample[n] = {}

            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.int_norm("Sn", "Sn", isotope, isotope_denom, beta)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.int_norm("Sn", "Sn", isotope, isotope_denom, beta, isotope_from_line1 = False)

        return data_sample
    # 2 no corr on isotope_denom,
    def internal_norm_1(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            beta = corr.beta(iter, "Sn", norm_ratio[0], norm_ratio[1], isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = False, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = False, isotope_from_line1 = False)
        return data_sample


    # 3 corr on isotope_denom
    def internal_norm_2(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            beta = corr.beta(iter, "Sn", norm_ratio[0], norm_ratio[1], isotope_denom_corr = True)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # External normalisation with Sb
    def external_norm_Sb(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            beta = corr.beta(iter, "Sb", norm_ratio[0], norm_ratio[1], isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    def beta(self, element, norm_ratio, iter):
        beta = {}
        for n in self.cycles:
            corr = int_norm(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            beta[n] = corr.beta(iter, element, norm_ratio[0], norm_ratio[1], isotope_denom_corr=False)

        return beta

    #outlier detection from https://github.com/joferkington/oost_paper_code/blob/master/utilities.py
    def mad_based_outlier(self, points, thresh=3.5):
        """
        Returns a boolean array with True if points are outliers and False
        otherwise.

        Parameters:
        -----------
            points : An numobservations by numdimensions array of observations
            thresh : The modified z-score to use as a threshold. Observations with
                a modified z-score (based on the median absolute deviation) greater
                than this value will be classified as outliers.

        Returns:
        --------
            mask : A numobservations-length boolean array.

        References:
        ----------
            Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
            Handle Outliers", The ASQC Basic References in Quality Control:
            Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
        """
        if len(points.shape) == 1:
            points = points[:,None]
        median = np.median(points, axis=0)
        diff = np.sum((points - median)**2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)
     
        modified_z_score = 0.6745 * diff / med_abs_deviation

        return modified_z_score > thresh

    #outlier rejection
    def mad_outlier_rejection(self, df):
        data_sample_outlier = pd.DataFrame()
        columns = df.columns.tolist()
        for column in columns:
            data_sample_outlier[column] = df[column].where(~self.mad_based_outlier(df[column]), other=np.NaN)
        return data_sample_outlier

    def mad_outlier_rejection_dict(self, dictionary):
        data_sample_outlier = {}
        for cycle in dictionary:
            temp = pd.DataFrame.from_dict(dictionary[cycle], orient = 'columns')
            temp = self.mad_outlier_rejection(temp)
            data_sample_outlier[cycle] = temp.to_dict()
        for cup in data_sample_outlier[cycle]:
            data_sample_outlier[cycle][cup] = Counter(data_sample_outlier[cycle][cup])
        return data_sample_outlier

    def z_score_outlier_rejection(self, df):
        df = df[np.abs(df - df.mean())<=(2*df.std())]

        return df


    #data processing wrap up
    def data_process(self, path, files, cup_config, isotopes, mass_range, corr_isotopes, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, option, iter_beta, isotope_denom):
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

            if option < 4:
                data_sample[sample] = methods[option](isotope_denom)

            else:
                data_sample[sample] = methods[option](isotope_denom, iter_beta)

            avg_ratio_sample_all[sample] = new_corr.avg_to_df(data_sample, sample)
            sd2_ratio_sample_all[sample] = new_corr.SD_to_df(data_sample, sample)
            sample_names[sample] = df.extract_metadata(sample, "Sample Name")

            data_all = new_corr.to_df_all(sample_names, avg_ratio_sample_all, sd2_ratio_sample_all, ratios = False)

        return data_all


    # standard-sample bracketing

    # to_dataframe
    def all_cycles_to_df(self, df, sample):
        df_zero = pd.DataFrame(df_zero["cycle1"])
        sample_name = df.extract_metadata(sample, "Sample Name")
        date = df.extract_metadata(sample, "Date")
        starttime = df.extract_metadata(sample, "Start Time")
        df_zero["date"] = datetime.strptime(date+starttime, '%d/%m/%Y%H:%M')
        df_zero["sample"] = sample_name
        df_zero["cycle"] = df_zero.index.values
        df_zero.index = df_zero["date"]
        df_zero = df_zero[columns]
        df_all = df_all.append(df_zero, ignore_index = True)

    def avg_to_df(self, data_sample, sample):

        df_ratio_sample = pd.DataFrame.from_dict(data_sample[sample], orient = 'index')
        df_ratio_sample = df_ratio_sample.sort_index(axis=1)
        df_ratio_sample = df_ratio_sample[np.abs(df_ratio_sample-df_ratio_sample.mean())<=(2*df_ratio_sample.std())]

        avg_ratio_sample= df_ratio_sample.mean().to_dict()
        return avg_ratio_sample

    def SD_to_df(self, data_sample, sample):

        df_ratio_sample = pd.DataFrame.from_dict(data_sample[sample], orient = 'index')
        df_ratio_sample = df_ratio_sample.sort_index(axis=1)
        df_ratio_sample = df_ratio_sample[np.abs(df_ratio_sample-df_ratio_sample.mean())<=(2*df_ratio_sample.std())]

        sd2_ratio_sample =  (df_ratio_sample.std()/math.sqrt(len(df_ratio_sample.index))).to_dict()
        return sd2_ratio_sample

    def to_df_all(self, sample_names, avg_ratio_sample_all, sd2_ratio_sample_all, ratios = True):
        avg_df_ratio_sample = pd.DataFrame.from_dict(avg_ratio_sample_all, orient = 'index')
        avg_df_ratio_sample= avg_df_ratio_sample.sort_index(axis=1)
        sd2_df_ratio_sample = pd.DataFrame.from_dict(sd2_ratio_sample_all, orient = 'index')
        sd2_df_ratio_sample= sd2_df_ratio_sample.sort_index(axis=1)
        sample_name_sample = pd.DataFrame.from_dict(sample_names, orient = 'index')
        sample_name_sample = sample_name_sample.sort_index(axis=1)
        final_df = avg_df_ratio_sample.merge(sd2_df_ratio_sample, on=avg_df_ratio_sample.index, how="outer")
        final_df = final_df.join(sample_name_sample, on="key_0", how="left")
        #final_df =  pd.merge(avg_df_ratio_sample, sd2_df_ratio_sample, on=avg_df_ratio_sample.index, how="outer")
        final_df = final_df.sort_index(axis=1)
        final_df.index = final_df['key_0']
        final_df = final_df.drop('key_0', axis=1)
        if ratios == True:
            if len(self.isotopes) > 1:
                final_df.columns = ["Sample Name", "112/120", "1SE", "114/120", "1SE", "115/120", "1SE", "116/120", "1SE", "117_2/120_2", "1SE", "117/120", "1SE", "118_2/120_2", "1SE", "118/120", "1SE", "119_2/120_2", "1SE", "119/120", "1SE", "122_2/120_2", "1SE", "122/120", "1SE", "124_2/120_2", "1SE", "124/120", "1SE"]
                #final_df.to_csv("2cycle_zero_120Sn_corr_std_error" + ".csv")
            else:
                final_df.columns = ["Sample Name", "112/120", "1SE", "114/120", "1SE", "115/120", "1SE", "116/120", "1SE", "117/120", "1SE", "118/120", "1SE", "119/120", "1SE", "122/120", "1SE", "124/120", "1SE",]
                #final_df.to_csv("2cycle_zero" + ".csv")
        return final_df