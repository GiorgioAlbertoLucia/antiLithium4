import yaml

import sys
sys.path.append('..')
from framework.src.plotter import Plotter

if __name__ == '__main__':

    #input_file = '/home/galucia/antiLithium4/analysis/config/cfgPlot.yml'
    #input_file = '/home/galucia/antiLithium4/analysis/config/cfgPlotPresentation.yml'
    #input_file = '/home/galucia/antiLithium4/analysis/figures/PAG_07112024/cfgPAG_07112024.yml'
    #input_file = '/home/galucia/antiLithium4/analysis/figures/03-01-2025-Non-Normalised/cfgNonNormalised.yml'
    #input_file = '/home/galucia/antiLithium4/analysis/figures/03-01-2025-Normalised-Grid/cfgNormalised.yml'
    #input_file = '/home/galucia/antiLithium4/analysis/figures/03-01-2025-Corrected-InvMass/cfgCorrected.yml'
    input_file = '/home/galucia/antiLithium4/analysis/figures/03-01-2025-Non-Corrected-InvMass/cfgNonCorrected.yml'

    with open(input_file, 'r') as f:
        config = yaml.safe_load(f)

    plotter = Plotter(config['outPath'])

    for plot in config['plots']:

        plotter.createCanvas(plot['axisSpecs'], **plot['canvas'])
        plotter.createMultiGraph(plot['axisSpecs'])
        if plot['legend']['bool']:  
            position = [plot['legend']['xmin'], plot['legend']['ymin'], plot['legend']['xmax'], plot['legend']['ymax']]
            plotter.createLegend(position, **plot['legend']['kwargs'])

        if 'graphs' in plot:
            for graph in plot['graphs']:
                plotter.addGraph(graph['inPath'], graph['graphName'], graph['graphLabel'], **graph['kwargs'])

        if 'hists' in plot.keys():
            for hist in plot['hists']:
                plotter.addHist(hist['inPath'], hist['histName'], hist['histLabel'], **hist['kwargs'])

        if 'ROIs' in plot.keys():
            for roi in plot['ROIs']:
                plotter.addROI(roi['lineSpecs'], roi['boxSpecs'], **roi['kwargs'])

        if 'lines' in plot.keys():
            for line in plot['lines']:
                plotter.addLine(line['lineSpecs'], **line['kwargs'])

        if 'funcs' in plot:
            for func in plot['funcs']:
                plotter.addFunc(func['inPath'], func['funcName'], func['funcLabel'], **func['kwargs'])
        
        if 'multigraph' in plot:
            plotter.drawMultiGraph(**plot['multigraph']['kwargs'])
        if plot['legend']['bool']:  
            plotter.drawLegend(**plot['legend']['kwargs'])

        plotter.save(plot['outPDF'])
    
    plotter.outFile.Close()