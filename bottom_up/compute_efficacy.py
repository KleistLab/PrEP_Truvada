#!/usr/bin/env python3
# Lanxin Zhang
import sys

from utils import *
from pe_truvada import *


def process_pe_to_phi(indicator, pe_array): 
   """
    Converts extinction probabilities into prophylactic efficacy estimates based on a given scenario.

    Parameters
    ----------
    indicator : list of bool
        A boolean list of 5 elements indicating the hypothesis/scenario:
        - The first four elements correspond to different assumptions (see `compute_pe_truvada`).
        - The fifth element indicates whether the scenario assumes 100% RAI.

    pe_array : list of float
        A list containing extinction probabilities for RAI and RVI returned by `compute_pe_truvada`.
    """
    r_rvi = 9.1e-5  # success rate of rvi for inoculum size (see Supplementary text 2)
    r_rai = 3.7e-3  # success rate of rai
    dict_v2p_rvi = {i: p_inoculum(r_rvi, i)[0] for i in range(20)}     # probability (value) of each inoculum size (key)
    dict_v2p_rai = {i: p_inoculum(r_rai, i)[0] for i in range(20)}
    pe0 = 0.9017        # extinction probability without PrEP
    pe_rai_0 = expected_pe(pe0, dict_v2p_rai) 
    pe_rvi_0 = expected_pe(pe0, dict_v2p_rvi)
    phi = None
    if indicator[3]:
        def mean_value_exposure_averaged(p_rai, p_rvi, ratio):  # calculate the averaged PE value based on RAI ratio
            return p_rai * ratio + (1 - ratio) * p_rvi
        ratio = 0.04    # ratio of rai/rvi among all sexual activities
        p0 = mean_value_exposure_averaged(pe_rai_0, pe_rvi_0, ratio)
        pe_rai = expected_pe(pe_array[0], dict_v2p_rai)
        pe_rvi = expected_pe(pe_array[1], dict_v2p_rvi)
        pe = mean_value_exposure_averaged(pe_rai, pe_rvi, ratio)
        phi = 1 - (1 - pe) / (1 - p0)
    elif indicator[4]: 
        pe = expected_pe(pe_array[0], dict_v2p_rai)
        phi = 1 - (1 - pe) / (1 - pe_rai_0)
    else:
        pe = expected_pe(pe_array[0], dict_v2p_rvi)
        phi = 1 - (1 - pe) / (1 - pe_rvi_0)
    return phi


def main():
    """
    Compute the PrEP efficacy of one certain dose for all 8 hypotheses
    """
    sys.stdout.write('If other dosing regimen is desired, please replace dosing.csv in ../data/ \n')
    dose = int(input('Please enter the number of doses per week: '))
    tend = int(input('Please enter the duration of regimen in days: '))
    adh = pd.read_csv('../data/dosing.csv', index_col=0).iloc[[7 - dose]].values
    indicators = [[1, i, j, 0] for i in range(2) for j in range(2)]
    for indicator in indicators[:4]:
        pe_rai, pe_rvi = compute_pe_truvada(indicator, adh_pattern=adh, tend=tend)
        key = ''.join(map(str, indicator))
        phi = process_pe_to_phi(indicator, [pe_rvi.numpy()])
        np.save("phi_{}_d{}".format(key, dose), phi)
        indicator2 = indicator[:3] + [1]
        key2 = ''.join(map(str, indicator2))
        phi2 = process_pe_to_phi(indicator2, [pe_rai.numpy(), pe_rvi.numpy()])
        np.save("phi_{}_d{}".format(key2, dose), phi2)


if __name__ == '__main__':
    main()
    print('Done')
