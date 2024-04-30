

particleMass = { 'He3': 2.80923,
                 'Pr': 0.938272,
                }

PID = { 0: {'part': 'E',
            'label': 'e'},            # electron
        1: {'part': 'Mu',
            'label': '#mu'},           # muon
        2: {'part': 'Pi',
            'label': '#pi'},           # pion
        3: {'part': 'K',
            'label': 'K'},            # kaon
        4: {'part': 'Pr',
            'label': 'p'},           # proton
        5: {'part': 'Deu',
            'label': 'd'},          # deuteron
        6: {'part': 'Tri',
            'label': '^{3}H'},          # triton
        7: {'part': 'He3',
            'label': '^{3}He'},          # helium3
        8: {'part': 'He4',
            'label': '^{4}He'},          # alpha
        9: {'part': 'Pi0',
            'label': '#pi^{0}'},          # pi0
        10: {'part': 'Photon',
            'label': '#gamma'},      # photon 
        11: {'part': 'K0',
            'label': 'K_{0}'},          # K0
        12: {'part': 'Lambda',
            'label': '#Lambda'},      # Lambda
        13: {'part': 'HyperTri',
            'label': '^{3}_{#Lambda}H'},    # HyperTriton
        14: {'part': 'HyperH4',
            'label': '^{4}_{#Lambda}H'},     # Hyperhydrog4
        15: {'part': 'Xi-',
            'label': '#Xi^{-}'},         # XiMinus
        16: {'part': 'Omega-',
            'label': '#Omega^{-}'},      # OmegaMinus
   }

PIDlabels = {key: value['label'] for key, value in PID.items()}