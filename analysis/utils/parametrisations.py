'''
    Dictionary of used parametrisations
'''

class ITS:
    '''
        Parametrisation for ITS response
    '''
    def __init__(self):
        pass

    # Parameter for proton expected cluster size 
    # cl size = kp1 / (betagamma)^kp2 + kp3
    pr_exp_params = {
            'kp1': 1.18941, 
            'kp2': 1.53792,
             'kp3': 1.69961,
        }
    
    # Parameter for proton resolution
    # res = res1 * erf((betagamma - res2) / res3)
    pr_res_params = {
            'kp1': 1.94669e-01, 
            'kp2': -2.08616e-01,
             'kp3': 1.30753,
        }
    
    # Parameter for helium expected cluster size 
    # cl size = kp1 / (betagamma)^kp2 + kp3
    he_exp_params = {
            'kp1': 2.35117, 
            'kp2': 1.80347,
             'kp3': 5.14355,
        }
    
    # Parameter for helium resolution
    # res = res1 * erf((betagamma - res2) / res3)
    he_res_params = {
            'res1': 8.74371e-02, 
            'res2': -1.82804,
             'res3': 5.06449e-01,
        }

class TOF:
    '''
        Parametrisation for TOF response
    '''
    def __init__(self):
        pass

    # Parameters for proton TOF resolution
    # res = res1 * exp(res2 * abs(pt))
    pr_res_params = {
            'res1': 1.22204e-02,
            'res2': 7.48467e-01,
        }

class PIDforTracking:
    '''
        Parametrisation ised to correct the pt of tracks associated with the wrong particle and thus receiving a wrong TPC correction
    '''

    def __init__(self):
        pass
    
    # Parameters to correct the pt of He3 wrongly tracked as H3
    # [kp1, kp2]
    # pt_corr = pt - pt*(kp1 + kp2*pt)
    # correction applied to He3 wrongly tracked as H3 with pt < ptmax
    he_params = {
        
        # param 2023
        'kp1': 0.3262,
        'kp2': -0.1263,
        'ptmax': 2.5
        
        # param 2024
        #'kp1': 0.1593,
        #'kp2': -0.0445,
        #'ptmax': 2.5
        }
    