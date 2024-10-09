import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from Bio import PDB
import seaborn as sns
import argparse

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments including PDB and PLY file paths.
    """
    parser = argparse.ArgumentParser(description='Process PDB and PLY files.')
    parser.add_argument('--pdb', type=str, required=True,
                        help='Path to the PDB file to be updated.')
    parser.add_argument('--ply', type=str, required=True,
                        help='Path to the PLY file to parse.')
    return parser.parse_args()

# Dictionary to convert three-letter codes to one-letter codes
three_to_one_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HIE': 'H', 'HID': 'H', 'ASH': 'D', 'GLH': 'E', 'CYX': 'C',
    'NME': 'X', 'ACE': 'X'
}

def parse_ply_file(filename):
    """
    Parse a PLY file to extract atom data.

    Args:
        filename (str): Path to the PLY file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed atom data.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    data = []
    header_ended = False

    for line in lines:
        if header_ended:
            components = line.strip().split()
            if len(components) == 9:
                x, y, z = map(float, components[:3])
                red, green, blue = map(float, components[3:6])
                values = float(components[6])
                atom = int(components[7])
                blurred_lvl = float(components[8])
                data.append([x, y, z, red, green, blue, values, atom, blurred_lvl])
        elif line.startswith('end_header'):
            header_ended = True

    df = pd.DataFrame(data, columns=['x', 'y', 'z', 'red', 'green', 'blue', 'values', 'atom', 'blurred_lvl'])
    return df

def update_pdb_bfactor(pdb_filename, clusters):
    """
    Update the B-factor column in a PDB file based on cluster information.

    Args:
        pdb_filename (str): Path to the PDB file to be updated.
        clusters (pd.DataFrame): DataFrame containing cluster information.
    """
    with open(pdb_filename, 'r') as file:
        pdb_lines = file.readlines()

    atom_to_cluster = {}
    for _, row in clusters.iterrows():
        cluster_id = int(row['cluster_id'].split('_')[1])
        for atom in row['atom']:
            atom_to_cluster[atom] = cluster_id

    updated_lines = []
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_number = int(line[6:11].strip())
            if atom_number in atom_to_cluster:
                cluster_number = atom_to_cluster[atom_number]
                bfactor = f"{cluster_number:.2f}".rjust(6)
                updated_line = line[:60] + bfactor + line[66:]
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)

    with open('updated_' + pdb_filename, 'w') as file:
        file.writelines(updated_lines)

def parse_pdb(pdb_file):
    """
    Parse a PDB file and extract residue information with maximum B-factor.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        list: List of dictionaries containing residue information.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                if PDB.is_aa(residue, standard=True) or res_name in three_to_one_dict:
                    b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
                    if b_factors:
                        max_b_factor = round(max(b_factors))
                    else:
                        max_b_factor = 0
                    residues.append({
                        'residue': three_to_one_dict.get(res_name, '?'),
                        'residue_number': residue.get_id()[1],
                        'chain_id': chain.get_id(),
                        'b_factor': max_b_factor
                    })

    return residues

# Main execution
if __name__ == "__main__":
    args = parse_args()

    df = parse_ply_file(args.ply)
    cluster_columns = ['red', 'green', 'blue']
    clusters = df.groupby(cluster_columns)['atom'].apply(lambda atoms: list(set(atoms))).reset_index()
    clusters['size'] = clusters['atom'].apply(len)
    clusters = clusters.sort_values(by='size', ascending=False).reset_index(drop=True)
    clusters['cluster_id'] = ['cluster_' + str(i) for i in range(len(clusters))]

    update_pdb_bfactor(args.pdb, clusters)

    pdb_file = 'updated_' + args.pdb
    residues = parse_pdb(pdb_file)

    chunk_size = 50
    chunks = [residues[i:i + chunk_size] for i in range(0, len(residues), chunk_size)]

    bfactor_values = sorted(set(residue['b_factor'] for residue in residues))
    color_palette = sns.color_palette("YlOrBr_r", len(bfactor_values) - 1)
    bfactor_color_map = {0: 'white'}
    bfactor_color_map.update(dict(zip(bfactor_values[1:], color_palette)))

    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=bfactor_color_map[b], markersize=10, label=f'Patch {b}') for b in bfactor_values]

    plt.figure(figsize=(2, 2))
    plt.legend(handles=legend_handles, title='B-factor', loc='center', frameon=False)
    plt.axis('off')
    plt.savefig('legend.png', bbox_inches='tight', dpi=300)
    plt.close()

    for idx, chunk in enumerate(chunks):
        sequence = [residue['residue'] for residue in chunk]
        sequence_numbers = [residue['residue_number'] for residue in chunk]
        b_factors = [residue['b_factor'] for residue in chunk]

        plt.figure(figsize=(12, 2))
        scatter = plt.scatter(sequence_numbers, [1]*len(sequence_numbers), c=[bfactor_color_map[b] for b in b_factors], marker='o')

        for residue, number in zip(sequence, sequence_numbers):
            plt.text(number, 1.02, residue, ha='center', va='bottom', fontsize=8)

        plt.yticks([])
        plt.xlabel('Residue Number')
        plt.title(f'Protein Sequence Part {idx + 1}')
        plt.grid(axis='x')

        tick_interval = 5
        plt.xticks(np.arange(min(sequence_numbers), max(sequence_numbers) + tick_interval, tick_interval))

        plt.tight_layout()
        plt.savefig(f'chunk_{idx + 1}.png', dpi=300)

    df = pd.DataFrame(residues)
    df_in_patch = df[df['b_factor'] != 0]
    df_not_in_patch = df[df['b_factor'] == 0]

    fig, ax = plt.subplots(figsize=(10, 6))

    def plot_stacked_text(ax, df, title):
        clusters = df.groupby('b_factor')

        for b_factor, group in clusters:
            y_position = 0
            for _, row in group.iterrows():
                ax.text(b_factor, y_position, f"{row['residue']} {row['residue_number']}",
                        ha='center', va='bottom', fontsize=8, color='black')
                y_position += 1

        ax.set_title(title)
        ax.set_xlabel('Patch ID)')
        ax.set_xlim(0, df['b_factor'].max() + 1)
        ax.set_xticks(range(1, df['b_factor'].max() + 1))
        #ax.set_ylim(0.0, 30) # adjust for individual systems!

    plot_stacked_text(ax, df_in_patch, 'Found in Patch')

    plt.tight_layout()
    plt.savefig('patch_res.png', dpi=300)

    df_not_in_patch.to_csv('nopatch.csv')
    df_in_patch.to_csv('patch.csv')

