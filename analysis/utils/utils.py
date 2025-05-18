'''
    Add to torchic
'''

from ROOT import TPaveText, TGraphErrors, TCanvas

def write_params_to_text(params: tuple, coordinates:tuple=[0.7, 0.9, 0.7, 0.9]) -> TPaveText:
    '''
        Write the parameters of the fit to the canvas

        Args:
            params: dictionary of the parameters to write (string: RooRealVar)
    '''

    text = TPaveText(coordinates[0], coordinates[1], coordinates[2], coordinates[3], 'NDC')
    text.SetFillColor(0)
    text.SetBorderSize(0)
    for param in params:
        if param.isConstant():
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} {param.getUnit()} (fixed)')
        else:
            text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
    return text


def create_graph(df, x: str, y: str, ex, ey, name:str='', title:str='') -> TGraphErrors:
        '''
            Create a TGraphErrors from the input DataFrame

            Parameters
            ----------
            x (str): x-axis variable
            y (str): y-axis variable
            ex (str): x-axis error
            ey (str): y-axis error
        '''

        # eliminate None values on x, y
        #df = df.filter(df[x].is_not_null())
        #df = df.filter(df[y].is_not_null())

        if len(df) == 0:
            return TGraphErrors()
        graph = TGraphErrors(len(df[x]))
        for irow, row in df.iterrows():
            graph.SetPoint(irow, row[x], row[y])
            xerr = row[ex] if ex != 0 else 0.
            yerr = row[ey] if ey != 0 else 0.
            graph.SetPointError(irow, xerr, yerr)
        
        graph.SetName(name)
        graph.SetTitle(title)

        return graph

class PdfCanvas:
     
    def __init__(self, pdf_path: str):
        self.pdf_path = pdf_path
        self.canvas = TCanvas('canvas', 'canvas', 800, 600)
        self.canvas.SetLogy(False)
        self.canvas.Print(f'{self.pdf_path}(')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.canvas.Clear()
        self.canvas.Print(f'{self.pdf_path})')
        self.canvas.Close()

    def draw_object(self, obj, **kwargs):
        
        if 'logy' in kwargs:
            self.canvas.SetLogy(kwargs['logy'])
        self.canvas.cd()
        self.canvas.Clear()
        obj.Draw(kwargs.get('draw_option', ''))
        self.canvas.Update()
        self.canvas.Print(self.pdf_path)
