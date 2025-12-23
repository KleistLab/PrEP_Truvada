import configparser
import os
from enum import Enum

import pandas as pd

# TODO: unit of IC50: TDF and FTC in umol, others in nmol, please unify the unit!!


class DrugClass(Enum):
    """
    Enum class to represent the class of drug
    """
    PI = 'PI'
    InI = 'InI'
    NRTI = 'NRTI'
    RTI = 'RTI'
    CRA = 'CRA'


class Drug(Enum):
    """
    Enum class to represent drug, contains drug class, molecular weight, IC50 and hill coefficient of
    each drug.
    Attributes:
        drug_class: DrugClass
        molecular_weight: float
            molecular weight (g/mol)
        IC_50: float
            drug concentration at the target site (blood or intracellular, ect.)
            where the target process is inhibited by 50%  (nM)
        m: float
            hill coefficient
    """
    MVC = (DrugClass.CRA, None, 5.06, 0.61)
    EFV = (DrugClass.RTI, 315.67, 10.8, 1.69)  # nM for IC50
    NVP = (DrugClass.RTI, None, 116, 1.55)
    DLV = (DrugClass.RTI, None, 336, 1.56)
    ETR = (DrugClass.RTI, None, 8.59, 1.81)
    RPV = (DrugClass.RTI, None, 7.73, 1.92)
    RAL = (DrugClass.InI, None, 25.5, 1.1)
    EVG = (DrugClass.InI, None, 55.6, 0.95)
    DTG = (DrugClass.InI, 419.38, 89, 1.3)      # nM for IC50
    ATV = (DrugClass.PI, None, 23.9, 2.69)
    APV = (DrugClass.PI, None, 262, 2.09)
    DRV = (DrugClass.PI, None, 45, 3.61)
    IDV = (DrugClass.PI, None, 130, 4.53)
    LPV = (DrugClass.PI, None, 70.9, 2.05)
    NFV = (DrugClass.PI, None, 327, 1.81)
    SQV = (DrugClass.PI, None, 88, 3.68)
    TPV = (DrugClass.PI, None, 483, 2.51)
    TDF = (DrugClass.NRTI, 635.5, 0.1, 1)       # uM for IC50
    TDF_vag = (DrugClass.NRTI, 635.5, 0.156909539547242, 1)    # used to change the ic50 of TDF
    TDF_col = (DrugClass.NRTI, 635.5, 0.0941663630292149, 1)
    FTC = (DrugClass.NRTI, 247.2, 0.85, 1)      # uM for IC50
    FTC_col = (DrugClass.NRTI, 247.2, 0.38, 1)      # used to change the IC50 of FTC-TP
    ISL = (DrugClass.NRTI, 293.258, 0.440029, 1)        # uM for IC50
    # 3TC = (DrugClass.NRTI, None, None, 1)

    def __init__(self, drug_class, molecular_weight, ic50, m):
        self.drug_class = drug_class
        self.molecular_weight = molecular_weight
        self.IC_50 = ic50
        self.m = m


class ViralDynamicParameter:
    """
    Class which contains the viral dynamic parameters
    """
    CL = None
    rho = None
    beta_T = None
    lambda_T = None
    delta_T = None
    delta_T1 = None
    delta_T2 = None
    delta_PIC_T = None
    k_T = None
    N_T = None
    CL_T = None
    T_u = None
    beta_M = None
    lambda_M = None
    delta_M = None
    delta_M1 = None
    delta_M2 = None
    delta_PIC_M = None
    k_M = None
    N_M = None
    CL_M = None
    M_u = None
    P_Ma4 = None
    P_La5 = None
    alpha = None
    zeta = None
    delta_l = None

    @classmethod
    def set_vd_parameters(cls, file='config.ini'):
        """
        Read the viral dynamic parameters from a file.
        :param
        file: str, name of file
        """
        config = configparser.ConfigParser()
        filepath = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '../config', file))
        config.read(filepath)
        cls.CL = float(config['viraldynamics']['CL'])
        cls.rho = float(config['viraldynamics']['rho'])
        cls.beta_T = float(config['viraldynamics']['beta_T'])
        cls.lambda_T = float(config['viraldynamics']['lambda_T'])
        cls.delta_T = float(config['viraldynamics']['delta_T'])
        cls.delta_T1 = float(config['viraldynamics']['delta_T1'])
        cls.delta_T2 = float(config['viraldynamics']['delta_T2'])
        cls.delta_PIC_T = float(config['viraldynamics']['delta_PIC_T'])
        cls.k_T = float(config['viraldynamics']['k_T'])
        cls.N_T = float(config['viraldynamics']['N_T'])
        cls.CL_T = (1 / cls.rho - 1) * cls.beta_T     # clearance rate of unsuccessful infected virus
        cls.T_u = cls.lambda_T / cls.delta_T        # steady state level of uninfected T-cells
        cls.beta_M = float(config['viraldynamics']['beta_M'])
        cls.lambda_M = float(config['viraldynamics']['lambda_M'])
        cls.delta_M = float(config['viraldynamics']['delta_M'])
        cls.delta_M1 = float(config['viraldynamics']['delta_M1'])
        cls.delta_M2 = float(config['viraldynamics']['delta_M2'])
        cls.delta_PIC_M = float(config['viraldynamics']['delta_PIC_M'])
        cls.k_M = float(config['viraldynamics']['k_M'])
        cls.N_M = float(config['viraldynamics']['N_M'])
        cls.CL_M = (1 / cls.rho - 1) * cls.beta_M  # clearance rate of unsuccessful infected virus
        cls.M_u = cls.lambda_M / cls.delta_M  # steady state level of uninfected T-cells
        cls.P_Ma4 = float(config['viraldynamics']['P_Ma4'])
        cls.P_La5 = float(config['viraldynamics']['P_La5'])
        cls.alpha = float(config['viraldynamics']['alpha'])
        cls.zeta = float(config['viraldynamics']['zeta'])
        cls.delta_l = float(config['viraldynamics']['delta_l'])


class ExtinctionCalculator:
    """
    Class to calculate the extinction probability with absence of drug.
    """
    @classmethod
    def get_pe_basic(cls):
        """
        Calculate the extinction probability with absence of drug for the basic viral dynamics.
        :return:
        PE: float
            extinction probability with absence of drug.
        """
        return 0.90177743

def calculate_propensity_constant(if_macrophage=False, if_reservoir=False):
    """
    Calculate the reaction propensities a1-a6 with absence of drug considering replicative capacity [/hr]
    if macrophage is considered, the reactions will be extended to 1-12
    """
    a = dict()
    a[1] = (ViralDynamicParameter.CL + ViralDynamicParameter.CL_T * ViralDynamicParameter.T_u) / 24    # v -> *
    a[2] = (ViralDynamicParameter.delta_PIC_T + ViralDynamicParameter.delta_T1) / 24  # T1 -> *
    a[3] = ViralDynamicParameter.delta_T2 / 24  # T2 -> *
    a[4] = ViralDynamicParameter.beta_T * ViralDynamicParameter.T_u / 24  # V -> T1
    a[5] = ViralDynamicParameter.k_T / 24  # T1 -> T2
    a[6] = ViralDynamicParameter.N_T / 24  # T2 -> T2 + V
    if if_reservoir and not if_macrophage:
        a[7] = ViralDynamicParameter.P_La5 * ViralDynamicParameter.k_T / 24         # T1 -> Tl
        a[5] = ViralDynamicParameter.k_T * (1 - ViralDynamicParameter.P_La5) / 24
    elif if_macrophage:
        a[1] = (ViralDynamicParameter.CL + ViralDynamicParameter.CL_T * ViralDynamicParameter.T_u +
                ViralDynamicParameter.CL_M * ViralDynamicParameter.M_u) / 24
        a[5] = ViralDynamicParameter.k_T * (1 - ViralDynamicParameter.P_La5) / 24
        a[7] = ViralDynamicParameter.P_La5 * ViralDynamicParameter.k_T / 24       # T1 -> Tl
        a[8] = ViralDynamicParameter.beta_M * ViralDynamicParameter.M_u / 24      # V -> M1
        a[9] = (ViralDynamicParameter.delta_PIC_M + ViralDynamicParameter.delta_M1) / 24      # M1 -> *
        a[10] = ViralDynamicParameter.delta_M2 / 24        # M2 -> *
        a[11] = ViralDynamicParameter.k_M / 24             # M1 -> M2
        a[12] = ViralDynamicParameter.N_M / 24             # M2 -> M2 + V
        if if_reservoir:
            a[13] = ViralDynamicParameter.delta_l / 24      # Tl -> *
            a[14] = ViralDynamicParameter.alpha / 24        # Tl -> T2
            a[15] = ViralDynamicParameter.zeta / 24         # Tl -> Tl + Tl
    return a

def calculate_propensities_for_drug_class(propensity_dict, eta_complement, drugclass, if_macrophage=False, if_reservoir=False):
    """
    Calculate the reaction propensities for the given drug effect eta for a certain drug class. The results will be replaced 
    in place in propensity_dict.
    eta_complement: the array of 1 - eta
    """
    if drugclass is DrugClass.InI:
        a5 = eta_complement * ViralDynamicParameter.k_T / 24
        if if_reservoir and not if_macrophage:
            propensity_dict[5] = a5 * (1 - ViralDynamicParameter.P_La5)
            propensity_dict[7] = a5 * ViralDynamicParameter.P_La5a
        elif if_macrophage:
            propensity_dict[7] = a5 * ViralDynamicParameter.P_La5
            propensity_dict[5] = a5 * (1 - ViralDynamicParameter.P_La5)
            propensity_dict[11] = eta_complement * ViralDynamicParameter.k_M / 24
        else:
            propensity_dict[5] = a5
    elif drugclass is DrugClass.CRA:    # TODO: wrong form of a1!!
        a4 = eta_complement * ViralDynamicParameter.beta_T * ViralDynamicParameter.T_u / 24
        a1 = (ViralDynamicParameter.CL + (1 / ViralDynamicParameter.rho - 1) * eta_complement * 
                  ViralDynamicParameter.beta_T * ViralDynamicParameter.T_u) / 24
        if if_macrophage:
            propensity_dict[1] = ((1 / ViralDynamicParameter.rho - 1) * eta_complement *  ViralDynamicParameter.beta_M * ViralDynamicParameter.M_u) / 24 + a1
            propensity_dict[4] = a4
            propensity_dict[8] = eta_complement * ViralDynamicParameter.beta_M * ViralDynamicParameter.M_u / 24
        else:
            propensity_dict[1], propensity_dict[4] = a1, a4
    elif drugclass is DrugClass.RTI or DrugClass.NRTI:
        a4 = eta_complement * ViralDynamicParameter.beta_T * ViralDynamicParameter.T_u / 24
        a1 = (ViralDynamicParameter.CL + (1 / ViralDynamicParameter.rho - eta_complement) *
                  ViralDynamicParameter.beta_T * ViralDynamicParameter.T_u) / 24
        if if_macrophage:
            propensity_dict[1] = ((1 / ViralDynamicParameter.rho - eta_complement) * ViralDynamicParameter.beta_M * ViralDynamicParameter.M_u) / 24 + a1
            propensity_dict[4] = a4
            propensity_dict[8] = eta_complement * ViralDynamicParameter.beta_M * ViralDynamicParameter.M_u / 24
        else:
            propensity_dict[1], propensity_dict[4] = a1, a4
    elif drugclass is DrugClass.PI:
        propensity_dict[6] = eta_complement * ViralDynamicParameter.N_T / 24
        if if_macrophage:
            propensity_dict[12] = eta_complement * ViralDynamicParameter.N_M / 24