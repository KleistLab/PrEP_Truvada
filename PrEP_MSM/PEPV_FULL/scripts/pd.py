import os

import pandas as pd
import torch
from scipy.interpolate import griddata

from .utils import DrugClass, calculate_propensities_for_drug_class, calculate_propensity_constant


class PharmacoDynamicsInterface(object):
    """
    Class to compute the reaction propensities
    Parameters:
    pk_objects: array-like
        an array contains all PK objects (corresponding to each drug)
    reservoir: boolean
        if reservoir is considered. If the case, add a7 and reduce a5 respectively
    macrophage: boolean
        if macrophage and latent is considered. If the case, add a7-a12 and modify a1 a5.
    file: str
        file name of possible file for PD computation. Currently only used for Truvada: file of MMOA matrix
    """
    # for Truvada  and Truvada + other RTI e.g. EFV
    drugCombination = {'FTC': {'TDF': 'Truvada'}, 'EFV':{'FTC': {'TDF':['Truvada', 'EFV']}}}

    def __init__(self, pk_objects, file=None, reservoir=False, macrophage=False):
        self._pk_objects = pk_objects
        self._file = file
        self._reservoir = reservoir
        self._macrophage = macrophage
        self._pd_objects = self._map_pk_to_pd()
        self._combine_pd()
        

    @staticmethod
    def check_drug_combination(drug_names):
        """
        Check if a certain combination exists for two (or more) drugs (e.g. Truvada)

        :parameter:
        drug_names: array-like
            list of str of drug names
        :return:
        name: str
            name of drug combination (None if no such combination exists)
        """
        # TODO: return the name of drug combination or None if no combination available (using self.drugCombination)
        drug_names.sort()
        n_drugs = len(drug_names)
        i = 0
        combi = PharmacoDynamicsInterface.drugCombination
        while i < n_drugs:
            if drug_names[i] in combi:
                combi = combi[drug_names[i]]
                i += 1

            else:
                return None
        return combi

    def _get_pd_class(self, pk_objects):
        """
        get the corresponding PD object for a drug/ drug combination represented by _pk_objects

        :parameter:
        _pk_objects: array-like
            array of pk objects
        :return:
        Pharmacodynamics_ : Pharmacodynamics object
            corresponding PD object
        """
        def _find_corresponding_pk_objects(name, pk_objects):
            # find the corresponding pk objects according to the drug name
            res = list()
            if name == 'Truvada':
                for object in pk_objects:
                    if object.regimen.get_drug_name() == 'TDF':
                        res.append(object)
                    if object.regimen.get_drug_name() == 'FTC':
                        res.append(object)
            else:
                for object in pk_objects:
                    if object.regimen.get_drug_name() == name:
                        res.append(object)
            return res

        def _return_pd_object_drug_class(pk_object, _reservoir, _macrophage):
            # return the corresponding PD object of each drug class for one PK object as a list
            if pk_object[0].regimen.get_drug_class() is DrugClass.CRA:
                return PharmacodynamicsCRA(pk_object, _reservoir, _macrophage)
            elif pk_object[0].regimen.get_drug_class() is DrugClass.InI:
                return PharmacodynamicsInI(pk_object, _reservoir, _macrophage)
            elif pk_object[0].regimen.get_drug_class() is DrugClass.PI:
                return PharmacodynamicsPI(pk_object, _reservoir, _macrophage)
            else:
                return PharmacodynamicsRTI(pk_object, _reservoir, _macrophage)


        if len(pk_objects) > 1:
            name = PharmacoDynamicsInterface.check_drug_combination([pk.regimen.get_drug_name() for
                                                                     pk in pk_objects])
            if not name:
                raise SystemExit('Data for this drug combination does not exist\n')
            else:
                # TODO call the corresponding PD for combination
                if isinstance(name, list):   # for Truvada + other RTI e.g. EFV
                    tmp = list()
                    for n in name:
                        pk_list = _find_corresponding_pk_objects(n, pk_objects)
                        if n == 'Truvada':
                            tmp.append(PharmacodynamicsTruvada(pk_list, self._file, self._reservoir, self._macrophage))
                        else: 
                            tmp.append(_return_pd_object_drug_class(pk_list, self._reservoir, self._macrophage))
                    return tmp
                elif name == 'Truvada':
                    return [PharmacodynamicsTruvada(pk_objects, self._file, self._reservoir, self._macrophage)]
        else:
            return [_return_pd_object_drug_class(pk_objects, self._reservoir, self._macrophage)]

    def _map_pk_to_pd(self):
        """
        generate a list which contains all PD objects corresponds to self._pk_objects.
        (if drug combination exists, generate one corresponding PD object)

        :return:
        pd_list: array-like
            list of pd objects
        """
        pd_list = list()
        pk_dict = dict()
        for pk_object in self._pk_objects:
            if pk_object.regimen.get_drug_class() is DrugClass.InI:
                if 5 not in pk_dict:
                    pk_dict[5] = list()         # drug that impact a5
                pk_dict[5].append(pk_object)
            elif pk_object.regimen.get_drug_class() is DrugClass.PI:
                if 6 not in pk_dict:
                    pk_dict[6] = list()         # drug that impact a5
                pk_dict[6].append(pk_object)
            else:
                if 1 not in pk_dict:
                    pk_dict[1] = list()  # drug that impact a1
                pk_dict[1].append(pk_object)
        for pk_list in pk_dict.values():
            pd_list += self._get_pd_class(pk_list)
        return pd_list

    def _combine_pd(self):
        """
        Update self._propensity_dict based on self._pd_objects, i.e. the pd objects for each drug/drug combination
        :return:
        a0: dict
            dictionary of propensities
        """
        # take the propensities of first PD object as base
        self._propensity_dict = self._pd_objects[0].get_propensities().copy()
        drugclass = self._pd_objects[0].get_drug_class()
        eta_complement = None
        pt_pd_sameclass_start = 0           # pointer to mark the place of same-type drugs 
        pt_pdsameclass_end = 0
        for pd_object in self._pd_objects[1:]:    # iterating from the second PD object (currently max. 2 objects)
            propensities = pd_object.get_propensities()
            for key, value in propensities.items():
                if isinstance(value, torch.Tensor):   # time-varying values
                    if not isinstance(self._propensity_dict[key], torch.Tensor):    # corresponding propensity in base is constant
                        self._propensity_dict[key] = value
                        pt_pd_sameclass_start += 1
                    else:   # in this case drugs belong to same class and work on same reaction propensities
                    # !! the input regimens have to place the drugs of same class together, and after all other drugs,
                    #  otherwise the eta can not be combined correctly. 
                    # raise SystemExit('Cannot compute drug effect for this combination\n')     
                        pt_pdsameclass_end += 1
        if pt_pd_sameclass_start < pt_pdsameclass_end:  # there are drugs of same classes
            eta_list = [pd_object.get_drug_effect() for pd_object in self._pd_objects[pt_pd_sameclass_start:pt_pdsameclass_end+1]]
            eta_complement = 1 - eta_list[0]
            for i in range(1, len(eta_list)):
                eta_complement *= 1 - eta_list[i]
            calculate_propensities_for_drug_class(self._propensity_dict, eta_complement, self._pd_objects[pt_pd_sameclass_start].get_drug_class(), self._macrophage, self._reservoir)
                        

    def get_propensity(self):
        return self._propensity_dict

    def get_pd_objects(self):
        return self._pd_objects


class AbstractPharmacodynamics(object):
    """
    Abstract class of pharmacodynamics. Aim is to compute the propensities a1-a6 according to each drug type.
     A strain with certain genotype can also be considered by giving the change of IC50, m and replicative capacity.

    Parameters
    -----------
    pk_objects: array-like
        array containing pk objects, each object contains the pk profile for one drug
        if drug combination is taken (actually only one drug will be accepted and handled now,
        this parameter is for further extension for drug combinations consist of more than
        one drug of same drug class, e.g. Truvada (two RTIs)
    reservoir: boolean
        indicate if latent reservoir is included. If True, there will be 7 reactions (viral dynamics can be extended)
        and new approach will be used for cumulative probability (PGS class: PgSystemReservoirNewApproach)
    macrophage: boolean
        indicate if the macrophage is included. If True, there will be 12 reactions
   Attributes
   ------------
   _pk_objects: array like
   _propensities: dict
        contain the propensity array a1-a6
    """
    def __init__(self, pk_objects, reservoir=False, macrophage=False):
        self._pk_objects = pk_objects
        self._macrophage = macrophage
        self._reservoir = reservoir
        self._propensities = calculate_propensity_constant(self._macrophage, self._reservoir)
        if len(self._pk_objects) == 1:
            pk_object = self._pk_objects[0]
            m = pk_object.regimen.get_hill_coefficient()
            ic_50 = pk_object.regimen.get_ic50()
            pk_profile = pk_object.get_concentration()
            self._eta = pk_profile[..., 0] ** m / (ic_50 ** m + pk_profile[..., 0] ** m)

    def _compute_distinct_propensities(self):
        """
        Compute the propensities that are impacted by the drug
        """
        raise NotImplementedError

    def get_drug_effect(self):
        return self._eta

    def get_propensities(self):
        return self._propensities

    def get_drug_class(self):
        return self._pk_objects[0].get_regimen().get_drug_class()


class PharmacodynamicsInI(AbstractPharmacodynamics):
    """
    Class of pharmacodynamics model of InI
    """
    def __init__(self, pk_objects, reservoir=False, macrophage=False):
        super().__init__(pk_objects, reservoir, macrophage)
        self._compute_distinct_propensities()

    def _compute_distinct_propensities(self):
        """
        Compute the propensity for InI
        """
        calculate_propensities_for_drug_class(self._propensities, 1 - self._eta, DrugClass.InI, self._macrophage, self._reservoir)
        

class PharmacodynamicsCRA(AbstractPharmacodynamics):
    """
    Class of pharmacodynamics model of CRA, for one mutant strain against CRA
    """
    def __init__(self, pk_objects, reservoir=False, macrophage=False):
        super().__init__(pk_objects, reservoir, macrophage)
        self._compute_distinct_propensities()

    def _compute_distinct_propensities(self):
        """
        Compute the propensity for CRA
        """
        calculate_propensities_for_drug_class(self._propensities, 1 - self._eta, DrugClass.CRA, self._macrophage, self._reservoir)
        


class PharmacodynamicsRTI(AbstractPharmacodynamics):
    """
    Class of pharmacodynamics model of RTI, for one mutant strain
    """
    def __init__(self, pk_objects, reservoir=False, macrophage=False):
        super().__init__(pk_objects, reservoir, macrophage)
        self._compute_distinct_propensities()

    def _compute_distinct_propensities(self):
        """
        Compute the propensity for RTI
        """
        calculate_propensities_for_drug_class(self._propensities, 1 - self._eta, DrugClass.RTI, self._macrophage, self._reservoir)
        


class PharmacodynamicsPI(AbstractPharmacodynamics):
    """
    Class of pharmacodynamics model of PI, for one mutant strain against PI
    """
    def __init__(self, pk_objects, reservoir=False, macrophage=False):
        super().__init__(pk_objects, reservoir, macrophage)
        self._compute_distinct_propensities()

    def _compute_distinct_propensities(self):
        """
        Compute the propensity for PI
        """
        calculate_propensities_for_drug_class(self._propensities, 1 - self._eta, DrugClass.PI, self._macrophage, self._reservoir)
        

class PharmacodynamicsTruvada(AbstractPharmacodynamics):
    """
    Class of pharmacodynamics model of Truvada
    """

    def __init__(self, pk_objects, mmoa_file=None, reservoir=False, macrophage=False):
        super().__init__(pk_objects, reservoir, macrophage)
        self._compute_distinct_propensities(file=mmoa_file)

    def _compute_distinct_propensities(self, file=None, interpolation=True):
        if self._pk_objects[0].regimen.get_drug_name() == 'FTC':
            c_ftc_new = self._pk_objects[0].get_concentration()[..., 0]
            c_tfv_new = self._pk_objects[1].get_concentration()[..., 0]
        else:
            c_ftc_new = self._pk_objects[1].get_concentration()[..., 0]
            c_tfv_new = self._pk_objects[0].get_concentration()[..., 0]
        if interpolation:
            mmoa_file = file or os.path.join(os.path.split(__file__)[0], '../Data/modMMOA_FTC_TDF_zero_extended.csv')
            try:
                conc_effect_matrix = pd.read_csv(mmoa_file)
            except FileNotFoundError:
                raise SystemExit('MMOA file for Truvada not found')
            c_ftc = torch.tensor(conc_effect_matrix['FTC'], dtype=torch.float32)
            c_tfv = torch.tensor(conc_effect_matrix['TFV'], dtype=torch.float32)
            c = torch.stack((c_ftc, c_tfv), dim=1)
            eta = torch.tensor(conc_effect_matrix['eps'], dtype=torch.float32)
            eta_new = griddata(c, eta, (c_ftc_new, c_tfv_new), method='cubic')
        else:
            c = torch.log10(torch.stack((c_ftc_new, c_tfv_new), dim=-1))
            model_file = file or os.path.join(os.path.split(__file__)[0], '../Data/model.pt')
            model = torch.nn.Sequential(
                torch.nn.Linear(2, 500),
                torch.nn.Sigmoid(),
                torch.nn.Linear(500, 100),
                torch.nn.Sigmoid(),
                torch.nn.Linear(100, 50),
                torch.nn.Sigmoid(),
                torch.nn.Linear(50, 1)
            )
            try:
                model.load_state_dict(torch.load(model_file))
            except FileNotFoundError:
                raise SystemExit('NN model not found')
            eta_new = model(c)
            eta_new = torch.clamp(eta_new, min=0, max=1)
        self._eta = 1 - torch.tensor(eta_new, dtype=torch.double)
        calculate_propensities_for_drug_class(self._propensities, 1 - self._eta, DrugClass.RTI, self._macrophage, self._reservoir)
