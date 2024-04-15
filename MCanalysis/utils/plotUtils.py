from ROOT import TH1F, TH2F

from .particle import PID


def setPIDlabels(hist, axis:str):
    '''
    Set the particle labels on the x or y axis of a histogram
    (visualize particle labels instead of PID indices).
    '''

    if axis == 'x':
        for PIDidx, PIDinfo in PID.items():
            hist.GetXaxis().SetBinLabel(PIDidx+1, PIDinfo['label'])
    elif axis == 'y':
        for PIDidx, PIDinfo in PID.items():
            hist.GetYaxis().SetBinLabel(PIDidx+1, PIDinfo['label'])
    else:   
        print('Invalid axis. Choose "x" or "y".')
    return hist