
from torchic.utils.terminal_colors import TerminalColors as tc
import pandas as pd
import uproot


if __name__ == '__main__':

    infile = '/data/galucia/lithium_local/MC/LHC25a4.root'
    table_names = ['O2he3hadtable', 'O2he3hadtablemc', 'O2he3hadmult']
    base = 'DF'

    outfile = uproot.recreate('input/mc.root')

    f = uproot.open(infile)
    keys = list(f.keys())
    _file_folders = [folder for folder in keys if (folder.startswith(base) and '/' not in folder)]
    file_folders_duplicated = [folder.split(';')[0] for folder in _file_folders] # list with potentially duplicated folders
    seen = {}
    for idx, val in enumerate(file_folders_duplicated):
        if val not in seen:
            seen[val] = idx
    file_folders = [_file_folders[idx] for idx in seen.values()]

    for folder in file_folders:
        dfs = []
        for table_name in table_names:
            table_path = f'{infile}:{folder}/{table_name}'
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+table_path+tc.RESET)
            dfs.append(uproot.open(table_path).arrays(library='pd'))

        df = pd.concat(dfs, axis=1)
        folder_clean = folder.split(';')[0]  # Clean folder name
        outfile[f'{folder_clean}/{table_names[0]}'] = df