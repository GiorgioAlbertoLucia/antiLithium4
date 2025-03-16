from ROOT import TFile, TCanvas

if __name__ == '__main__':

    infile = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/purity.root')
    tofdir = infile.Get('Pr_TOF')

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    incanvases = []

    list = tofdir.GetListOfKeys()
    entries = list.GetEntries()
    for ientry in range(entries):

        key = list.At(ientry)
        if 'cNSigmaTOF' not in key.GetName():
            continue

        obj = tofdir.Get(key.GetName())
        list_prim = obj.GetListOfPrimitives()
        c_list = {}

        for ikey in range(list_prim.GetEntries()):
            print(f'  {list_prim.At(ikey).GetName()=}')
            prim = list_prim.At(ikey)
            if 'frame' in prim.GetName():
                pt_name = obj.GetName().replace('cNSigmaTOF', '')
                pt_name = pt_name[3:7]
                pt = float(pt_name)
                title = f'{(pt-0.05):.2f} < #it{{p}}_T < {(pt+0.05):.2f} GeV/c; n#sigma_{{TOF}}'
                prim.SetTitle(title)
            if 'paramBox' in prim.GetName():
                prim.SetX1NDC(0.45)
                prim.SetX2NDC(0.7)
                prim.SetY1NDC(0.15)
                prim.SetY2NDC(0.6)
                prim.SetTextSize(0.03)
                prim.SetTextFont(40)
            c_list[prim.GetName()] = prim

        incanvases.append(c_list)

    #canvas.Divide(2, entries // 2 + 1)
    for i, c_list in enumerate(incanvases):
        
        print(f'\nCanvas {i + 1}')
        #canvas.cd(i + 1)
        #ipad = TPad(f'pad{i}', f'pad{i}', 0, 0, 1, 1)
        #ipad.Draw()
        #ipad.cd()
        canvas.Clear()
        canvas.cd()
        for _name, prim in c_list.items():
            print(f'  {prim.GetName()}')
            if 'frame' in prim.GetName():
                prim.Draw()
            else:
                prim.Draw('same')
        #ipad.SetLogy()
        canvas.SetLogy()
        #ipad.Update()
        if i == 0:
            canvas.Print('tof_purity.pdf(')
        elif i == len(incanvases) - 1:
            canvas.Print('tof_purity.pdf)')
        else:
            canvas.Print('tof_purity.pdf')

    canvas.Update()
    #canvas.SaveAs('tof_purity.pdf')