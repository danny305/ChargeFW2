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

    # TODO finish these and integrate them into main
    parser.add_argument('--keep-cations', required=False,)

    parser.add_argument('--keep-water', required=False,)

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
        loop = table.loop
    else: 
        loop = block.init_loop(prefix, tags)

    for atom in rm_atoms:
        loop.add_row(atom)



def filter_element_atoms_to_remove(atom_site, allowed_elements=['C','H','N','O','S','P','F','CL','BR','I']):
    '''
        atom_site row list indice
        index 5  = _atom_site.label_comp_id
        index 17 = _atom_site.auth_comp_id
    
    '''
    filtered_rows = [
        (idx, atom) 
        for idx, atom in enumerate(atom_site) 
            if not atom[2] in allowed_elements
    ]

    rm_idx, rm_rows = [], []
    if filtered_rows:
        rm_idx, rm_rows = zip(*filtered_rows)

    return rm_idx, rm_rows


def filter_residue_atoms_to_remove(atom_site, rm_residues, label_comp_id=False):
    '''
        atom_site row list indice
        index 5  = _atom_site.label_comp_id
        index 17 = _atom_site.auth_comp_id
    
    '''
    if isinstance(rm_residues, (str, tuple, list, set)):
        if isinstance(rm_residues, str):
            rm_residues = [rm_residues]
        if not isinstance(rm_residues, set):
            rm_residues = set(rm_residues)
    else:
        raise TypeError(
            f'rm_residues parameter must be either a str of 1 comp_id" \
            "or a list/tuple/set of comp_id. rm_residues type: {type(rm_residues)}'
        )

    label_idx = 5
    auth_idx  = 17

    if label_comp_id:
        comp_idx = label_idx
    else: 
        comp_idx = auth_idx

    filtered_rows = [
        (idx, atom) 
        for idx, atom in enumerate(atom_site) 
            if atom[comp_idx] in rm_residues
    ]

    rm_idx, rm_rows = [], []
    if filtered_rows:
        rm_idx, rm_rows = zip(*filtered_rows)

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



    # Remove cation atoms
    rm_elem_idx, rm_elem_atoms = filter_element_atoms_to_remove(atom_site, param_elements)

    rm_elem_counts = Counter([row[2] for row in rm_elem_atoms])

    append_rm_atom_site_category(block, rm_elem_atoms, atom_site.tags)

    rm_atoms_from_atom_site(atom_site, rm_elem_idx)



    # Remove water atoms
    rm_water_idx, rm_water_atoms = filter_residue_atoms_to_remove(atom_site, ['HOH'])

    label_idx = 5
    auth_idx  = 17

    comp_idx = auth_idx
    rm_water_counts = Counter([row[comp_idx] for row in rm_water_atoms])

    append_rm_atom_site_category(block, rm_water_atoms, atom_site.tags)

    rm_atoms_from_atom_site(atom_site, rm_water_idx)



    # Write out a file if atoms were removed
    if rm_elem_counts or rm_water_atoms:
        rm_list = []

        atoms_removed = Counter()
        if rm_elem_counts:
            atoms_removed.update(rm_elem_counts)

            rm_list.extend([elem.capitalize() for elem in rm_elem_counts.keys()])

        else: 
            print('No cations were removed.')

        if rm_water_counts:
            atoms_removed.update(rm_water_counts)

            rm_list.extend([comp.upper() for comp in rm_water_counts.keys()])

        else: 
            print('No waters were removed.')

        print(f'Total Number of rows removed from atom_site: {atoms_removed}')

        rm_elem_str = '_rm_' + '_'.join(rm_list) + '.'

        pdb_name, pdb_ext = str(args.cif_file).split('.', 1)

        out_cif = pdb_name + rm_elem_str + pdb_ext

        doc.write_file(out_cif, gemmi.cif.Style.Pdbx)

        print(f'writing out: {out_cif}')
