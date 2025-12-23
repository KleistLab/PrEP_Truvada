import numpy as np

trial_parameters_drug_det_men ={

    'iPrEx':{
        'individuals': 624, 
        'avg_obs_time': 1.33,
        'n_inf':3,
        'incidence': (3/830),
        'py':830
        
        },

    'IPERGAY':{
        'individuals': 171, 
        'avg_obs_time': 1.11,
        'n_inf':0,
        'incidence': (0/190),
        'py':190},
    'HPTN083':{
        'individuals': 1932, 
        'avg_obs_time': 1.39,
        'n_inf':6,
        'incidence': (6/2685),
        'py':2685
        }, #excluding TGW?
        
    'PURPOSEII':{
    'individuals': 727, 
    'incidence': 2/654, 
    'avg_obs_time': 0.90,
    'py':654,
    'n_inf': 2
         },
    'DISCOVER':{
    'individuals': 2555, 
    'incidence': 1/4198, 
    'avg_obs_time': 4373/2661,
    'py':4198,
    'n_inf': 1
            }
        
        }


trial_parameters_drug_undet_men = {

    'iPrEx':{
        'individuals': 600, 
        'avg_obs_time': 1.33,
        'n_inf':31,
            'incidence': (31/798),
        'py':798
        
        },
 #excluding TGW?
    'IPERGAY':{
        'individuals': 28, 
        'avg_obs_time': 1.11,
        'n_inf':2,
        'incidence': (2/31),
        'py':31},

    'HPTN083':{
        'individuals': 315, 
        'avg_obs_time': 1.39,
        'n_inf':33,
        'incidence': (33/438),
        'py':438
        },
    'PURPOSEII':{
    'individuals': 193, 
    'incidence': 6/174, 
    'avg_obs_time': 0.90,
            'n_inf':6,
            'py':174
        },
    'DISCOVER':{
    'individuals': 106, 
    'incidence': 10/175, 
    'avg_obs_time': 4373/2661,
    'py': 175,
    'n_inf': 10  }
        }


df_035 = [30.97981547619050,56.921959077381000,74.13786160714280,
          85.41745684523810,92.36096726190500,96.44264360119090,98.67986086309570]

df_001 = [58.10755877976190,83.71472767857150,93.01222767857160,
          96.99397247023840,98.6081287202387,99.43083928571490,99.896758928572]

hptn083_pAdh = np.array([0.09967141, 0.05421687, 0.05421687, 0.14129244, 0.14129244,
        0.14129244, 0.36801752])
purpose2_pAdh = np.array([0.11982083, 0.04143337, 0.04143337, 0.19932811, 0.19932811,
        0.19932811, 0.19932811])
discover_pAdh = np.array([0.04380475594493116, 0.07759699624530664/2, 0.07759699624530664/2, 0.8785982478097623/4,0.8785982478097623/4,0.8785982478097623/4,0.8785982478097623/4])

ipergay_pDet_001, iprex_pDet_035 = np.array(df_001)/sum(df_001), np.array(df_035)/sum(df_035)
