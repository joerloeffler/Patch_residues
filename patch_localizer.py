import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import PDB
import seaborn as sns
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process PDB and PLY files with patch size filtering (atom-based B-factor writing).')
    parser.add_argument('--pdb', type=str, required=True, help='Path to the PDB file to be updated.')
    parser.add_argument('--ply', type=str, required=True, help='Path to the PLY file to parse.')
    parser.add_argument('--chunk-size', type=int, default=50, help='Residues per chunk for sequence plotting.')
    parser.add_argument('--patch-size', type=int, default=10, help='Minimum number of residues required for a patch.')
    return parser.parse_args()

three_to_one_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HIE': 'H', 'HID': 'H', 'ASH': 'D', 'GLH': 'E', 'CYX': 'C',
    'NME': 'X', 'ACE': 'X'
}

def parse_ply_file(filename):
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
    return pd.DataFrame(data, columns=['x', 'y', 'z', 'red', 'green', 'blue', 'values', 'atom', 'blurred_lvl'])

def update_pdb_bfactor_atom_level_filtered(pdb_filename, atom_to_cluster, valid_patches):
    with open(pdb_filename, 'r') as file:
        pdb_lines = file.readlines()
    updated_lines = []
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_number = int(line[6:11].strip())
            cluster_id = atom_to_cluster.get(atom_number, 0)
            bfactor_value = cluster_id if cluster_id in valid_patches else 0
            bfactor = f"{bfactor_value:.2f}".rjust(6)
            updated_line = line[:60] + bfactor + line[66:]
            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)
    with open('updated_' + pdb_filename, 'w') as file:
        file.writelines(updated_lines)

def parse_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                if PDB.is_aa(residue, standard=True) or res_name in three_to_one_dict:
                    b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
                    max_b_factor = round(max(b_factors)) if b_factors else 0
                    chain_id = chain.get_id() or ' '
                    residues.append({
                        'residue': three_to_one_dict.get(res_name, '?'),
                        'residue_number': residue.get_id()[1],
                        'chain_id': chain_id,
                        'b_factor': max_b_factor
                    })
    return pd.DataFrame(residues)

def count_residues_per_patch(atom_to_cluster, pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    residue_patch_counts = {}

    for model in structure:
        for chain in model:
            chain_id = chain.get_id() or ' '
            for residue in chain:
                res_seq = residue.get_id()[1]
                res_id = (chain_id, res_seq)
                atom_numbers = [atom.get_serial_number() for atom in residue.get_atoms()]
                clusters = [atom_to_cluster.get(atom_number, 0) for atom_number in atom_numbers if atom_number in atom_to_cluster]
                patch_ids = set([cluster for cluster in clusters if cluster != 0])
                for patch_id in patch_ids:
                    residue_patch_counts.setdefault(patch_id, 0)
                    residue_patch_counts[patch_id] += 1

    return residue_patch_counts

if __name__ == "__main__":
    args = parse_args()
    df = parse_ply_file(args.ply)
    clusters = df.groupby(['red', 'green', 'blue'])['atom'].apply(lambda atoms: list(set(atoms))).reset_index()
    clusters['size'] = clusters['atom'].apply(len)
    clusters = clusters.sort_values(by='size', ascending=False).reset_index(drop=True)
    clusters['cluster_id'] = ['cluster_' + str(i) for i in range(len(clusters))]

    atom_to_cluster = {atom: int(row['cluster_id'].split('_')[1]) for _, row in clusters.iterrows() for atom in row['atom']}

    residue_patch_counts = count_residues_per_patch(atom_to_cluster, args.pdb)
    valid_patches = [patch_id for patch_id, count in residue_patch_counts.items() if count >= args.patch_size]

    update_pdb_bfactor_atom_level_filtered(args.pdb, atom_to_cluster, valid_patches)

    pdb_file = 'updated_' + args.pdb
    df_residues = parse_pdb(pdb_file)

    filtered_df_in_patch = df_residues[df_residues['b_factor'].isin(valid_patches)].sort_values(by='residue_number').reset_index(drop=True)
    df_not_in_patch = df_residues[~df_residues['b_factor'].isin(valid_patches)]

    filtered_df_in_patch.to_csv('patch.csv', index=False)
    df_not_in_patch.to_csv('nopatch.csv', index=False)
    small_patches = pd.DataFrame({'b_factor': list(residue_patch_counts.keys()), 'residue_count': list(residue_patch_counts.values())})
    small_patches = small_patches[~small_patches['b_factor'].isin(valid_patches)]
    small_patches.to_csv('small_patches.csv', index=False)

    bfactor_values = sorted(set(filtered_df_in_patch['b_factor']))
if len(bfactor_values) == 0:
    print("No valid patches survived the filtering. Exiting without plotting.")
    exit(0)

chunk_size = args.chunk_size
chunks = [filtered_df_in_patch.iloc[i:i + chunk_size] for i in range(0, len(filtered_df_in_patch), chunk_size)]
color_palette = sns.color_palette("YlOrBr_r", len(bfactor_values))
bfactor_color_map = dict(zip(bfactor_values, color_palette))

legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=bfactor_color_map[b], markersize=10, label=f'Patch {b}') for b in bfactor_values]

plt.figure(figsize=(2, 2))
plt.legend(handles=legend_handles, title='Patch ID', loc='center', frameon=False)
plt.axis('off')
plt.savefig('legend.png', bbox_inches='tight', dpi=300)
plt.close()

for idx, chunk in enumerate(chunks):
    sequence = chunk['residue']
    sequence_numbers = chunk['residue_number']
    b_factors = chunk['b_factor']

    plt.figure(figsize=(12, 2))
    plt.scatter(sequence_numbers, [1]*len(sequence_numbers), c=[bfactor_color_map[b] for b in b_factors], marker='o')
    for residue, number in zip(sequence, sequence_numbers):
        plt.text(number, 1.02, residue, ha='center', va='bottom', fontsize=8)
    plt.yticks([])
    plt.xlabel('Residue Number')
    plt.title(f'Protein Sequence Part {idx + 1} (Filtered Patches)')
    plt.grid(axis='x')
    tick_interval = 5
    plt.xticks(np.arange(min(sequence_numbers), max(sequence_numbers) + tick_interval, tick_interval))
    plt.tight_layout()
    plt.savefig(f'chunk_{idx + 1}.png', dpi=300)

fig, ax = plt.subplots(figsize=(10, 6))
clusters = filtered_df_in_patch.groupby('b_factor')
y_offset = 0
for b_factor, group in clusters:
    sorted_group = group.sort_values(by='residue_number')
    for _, row in sorted_group.iterrows():
        ax.text(b_factor, y_offset, f"{row['residue']} {row['residue_number']}", ha='center', va='bottom', fontsize=8, color='black')
        y_offset += 1
    y_offset += 2

ax.set_title('Residues Found in Filtered Patches')
ax.set_xlabel('Patch ID')
ax.set_xlim(min(bfactor_values) - 1, max(bfactor_values) + 1)
ax.set_ylim(-1, y_offset)
plt.tight_layout()
plt.savefig('patch_res.png', dpi=300)
