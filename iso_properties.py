__author__ = 'matthias'

from toposort import toposort, toposort_flatten

class Isotope_Masses():

    def __init__(self):
        self.isotope_masses = {}

    def add_masses_dict(self, masses_dict):
        self.isotope_masses = masses_dict

    def add_single_isotope(self, isotope_mass_no, isotope_mass):
        self.isotope_masses[isotope_mass_no] = isotope_mass

    def get_Isotope_mass(self, isotope_mass_no):
        return self.isotope_masses[isotope_mass_no]

    def get_all_Isotope_masses(self):
        return self.isotope_masses

class Isotope_Ratios():

    def __init__(self):
        self.isotope_ratios = {}
        self.isotope_abundances = {}

    def add_isotope_ratio(self, isotope_nom, isotope_denom, value):
        if isotope_denom in self.isotope_ratios:
            self.isotope_ratios[isotope_denom][isotope_nom] = value
        else:
            self.isotope_ratios[isotope_denom] = {isotope_nom : value}

        self.isotope_ratios[isotope_denom][isotope_nom] = value

    def add_ratios_dict(self, isotope_denom, ratios_dict):
        self.isotope_ratios[isotope_denom] = ratios_dict

    def remove_ratio(self, isotope_nom, isotope_denom):
        self.isotope_ratios[isotope_denom].pop(isotope_nom)

    def transf_ratio(self, isotope_denom_old, isotope_denom_new):
        for isotope in self.isotope_ratios[isotope_denom_old]:
            trans_ratio = self.isotope_ratios[isotope_denom_old][isotope] / self.isotope_ratios[isotope_denom_old][isotope_denom_new]
            self.add_isotope_ratio(isotope, isotope_denom_new, trans_ratio)

        invers_ratio = 1 / self.isotope_ratios[isotope_denom_old][isotope_denom_new]
        self.add_isotope_ratio(isotope_denom_old, isotope_denom_new, invers_ratio)

        self.remove_ratio(isotope_denom_new, isotope_denom_new)

    def calc_abundances(self, isotope_denom):
        sum_ratios = 0
        for isotope in self.get_all_ratios(isotope_denom):
            sum_ratios += self.get_ratio(isotope, isotope_denom)
        for isotope in self.get_all_ratios(isotope_denom):
            self.isotope_abundances[isotope] = 100/(sum_ratios + 1) * self.get_ratio(isotope, isotope_denom)

        abundances_isotope_denom = 100
        for isotope in self.isotope_abundances:
            abundances_isotope_denom -= self.isotope_abundances[isotope]
        self.isotope_abundances[isotope_denom] = abundances_isotope_denom

    def get_ratio(self, isotope_nom, isotope_denom):
        return self.isotope_ratios[isotope_denom][isotope_nom]

    def get_all_ratios(self, isotope_denom):
        return self.isotope_ratios[isotope_denom]

    def get_abundance(self, isotope, isotope_denom):
        self.calc_abundances(isotope_denom)
        return self.isotope_abundances[isotope]

    def get_all_abundances(self, isotope_denom):
        self.calc_abundances(isotope_denom)
        return self.isotope_abundances

class Isotope_Abundances(object):

    def __init__(self):
        self.isotope_abundances = {}

    def calc_abundances(self, cls, isotope_denom):
        sum_ratios = 0
        for isotope in cls.get_all_ratios(isotope_denom):
            sum_ratios += cls.get_ratio(isotope, isotope_denom)
        for isotope in cls.get_all_ratios(isotope_denom):
            self.isotope_abundances[isotope] = 100/(sum_ratios + 1) * cls.get_ratio(isotope, isotope_denom)

        abundances_isotope_denom = 100
        for isotope in self.isotope_abundances:
            abundances_isotope_denom -= self.isotope_abundances[isotope]
        self.isotope_abundances[isotope_denom] = abundances_isotope_denom

    def get_abundance(self, isotope):
        return self.isotope_abundances[isotope]

    def get_all_abundances(self):
        return self.isotope_abundances

class Isotopes_Mass_Range():

    def __init__(self):
        self.isotopes_mass_range = {}

    def add_mass_range_dict(self, mass_range_dict):
        self.isotopes_mass_range = mass_range_dict

    def add_single_mass_isotope(self, isotope_mass_no, isotopes):
        self.isotopes_mass_range[isotope_mass_no] = isotopes

    def get_isotopes(self, isotope_mass_no):
        return self.isotopes_mass_range[isotope_mass_no]

    def get_mass_range(self):
        return self.isotopes_mass_range

    def get_interferences(self, isotope_mass_no, Element):
        isotopes = self.isotopes_mass_range[isotope_mass_no]
        isotopes.remove(Element)
        return isotopes

    # Graph of dependencies
    def get_graph_of_corr(self, corr_isotopes):
        graph = {}
        for mass_no in self.isotopes_mass_range:
            isotopes = self.isotopes_mass_range[mass_no]
            graph[mass_no] = set()
            for isotope in isotopes:
                if isotope in corr_isotopes and set(isotopes).issubset(corr_isotopes) == False:
                    graph[mass_no].add(corr_isotopes[isotope])
                elif isotope not in corr_isotopes:
                    None
                elif len(self.isotopes_mass_range[corr_isotopes[isotope]]) == 1 and len(self.isotopes_mass_range[mass_no]) > 1:
                    graph[mass_no].add(corr_isotopes[isotope])
        return graph
    # Order the dependencies - directed topology
    def get_order_of_corr(self, corr_isotopes):
        return list(toposort(self.get_graph_of_corr(corr_isotopes)))
