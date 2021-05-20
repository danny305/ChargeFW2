import gemmi
import argparse
import json
from pathlib import Path
from collections import Counter
from pprint import pprint



def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('-c','--cif-file', required=True, type=str)

    parser.add_argument('-d', '--cif-dir', required=False, type=Path, 
        default=Path('./input')
    )
    parser.add_argument('-j','--params-file', required=False, type=str, 
        default='EEM_10_Cheminf_b3lyp_aim.json'
    )
    parser.add_argument('-p','--params-dir', required=False, type=Path, 
        default=Path('../data/parameters')
    )

    args = parser.parse_args()

    args.cif_file = args.cif_dir / args.cif_file

    args.params_file = args.params_dir / args.params_file

    if args.cif_file.suffix != '.cif':
        raise ValueError('cif_file parameter must end with a .cif extension.')

    if args.params_file.suffix != '.json':
        raise ValueError('params_file parameter must end with a .json extension.')

    return args



def elements_from_params(params):

    if params.exists():
        with params.open('rt') as f:
            params_data = json.load(f)
    else:
        raise FileNotFoundError(f'Parameter file passed does not exist: {params}')

    atom_data = params_data['atom']['data']

    params_elements = {param['key'][0].upper() for param in atom_data}

    return params_elements



def append_rm_atom_site_category(block, rm_atoms, tags):

    prefix = '_rm_atom_site.'

    table = block.find(prefix, tags)
    if table:
        loop = table.get_loop()
    else: 
        loop = block.init_loop(prefix, tags)

    for atom in rm_atoms:
        loop.add_row(atom)



def filter_atoms_to_remove(atom_site, allowed_elements=['C','H','N','O','S','P','F','CL','BR','I']):

    rm_idx, rm_rows = zip(*[
        (idx, atom) 
        for idx, atom in enumerate(atom_site) 
            if not atom[2] in allowed_elements
    ])

    return rm_idx, rm_rows



def rm_atoms_from_atom_site(atom_site, rm_idx):

    for idx in rm_idx[::-1]:
        atom_site.remove_row(idx)



if __name__=="__main__":

    args = cli()

    param_elements = elements_from_params(args.params_file)

    print(f'Elements found in params file: {param_elements}')

    doc = gemmi.cif.read_file(str(args.cif_file))

    block = doc.sole_block()

    atom_site = block.find_mmcif_category("_atom_site.")

    rm_idx, rm_atoms = filter_atoms_to_remove(atom_site, param_elements)

    rm_elem_counts = Counter([row[2] for row in rm_atoms])

    append_rm_atom_site_category(block, rm_atoms, atom_site.tags)

    rm_atoms_from_atom_site(atom_site, rm_idx)

    if rm_elem_counts:
        print(f'Atoms rows removed by element: ', rm_elem_counts)

        rm_elem_str = '_rm_' + '_'.join([elem.capitalize() for elem in rm_elem_counts.keys()]) + '.'

        pdb_name, pdb_ext = str(args.cif_file).split('.', 1)

        out_cif = pdb_name + rm_elem_str + pdb_ext

        doc.write_file(out_cif, gemmi.cif.Style.Pdbx)

        print(f'writing out: {out_cif}')

    else: 
        print('No Elements were removed.')
