import yaml

import sys
sys.path.append('..')
from framework.src.plotter import Plotter

if __name__ == '__main__':

    with open('/home/galucia/antiLithium4/analysis/config/cfgPlot.yml', 'r') as f:
        config = yaml.safe_load(f)

    plotter = Plotter(config['outPath'])

    for plot in config['plots']:

        plotter.createCanvas(plot['axisSpecs'])

        for hist in plot['hists']:
            plotter.addHist(hist['inPath'], hist['histName'], hist['histLabel'], **hist['kwargs'])

        for roi in plot['ROIs']:
            plotter.addROI(roi['lineSpecs'], roi['boxSpecs'], **roi['kwargs'])
        
        if plot['legend']['bool']:  
            position = [plot['legend']['xmin'], plot['legend']['ymin'], plot['legend']['xmax'], plot['legend']['ymax']]
            plotter.addLegend(position, **plot['legend']['kwargs'])

        plotter.save(plot['outPDF'])
    
    plotter.outFile.Close()