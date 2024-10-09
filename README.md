# Patch_residues
Map hydrophobic patches to residues using the ply files generated from PEP-patch

This script processes Protein Data Bank (PDB) files and PLY files to update the B-factor of atoms based on clustering information extracted from the PLY file. The script also generates visualizations of the protein's residue sequences, highlighting residue positions with non-zero B-factors (i.e., those that belong to identified clusters).
Dependencies

The script requires the following Python packages:

    pandas
    numpy
    matplotlib
    seaborn
    argparse
    biopython

    pip install pandas numpy matplotlib seaborn biopython

The script accepts two command-line arguments:

    --pdb: Path to the PDB file that will be updated with cluster information.
    --ply: Path to the PLY file containing cluster data.

Main Execution Flow

    1) PLY File Parsing: The script reads the PLY file to extract atom information such as coordinates and cluster-related data.
    2) Cluster Grouping: Atoms are grouped into clusters based on their color (RGB) values. Each cluster is assigned an ID and its size is calculated.
    3) PDB File Update: The script updates the B-factors in the original PDB file based on the cluster assignments and saves the new PDB file with the updated B-factors.
    4) Residue Visualization: The updated PDB file is parsed, and the residues are visualized in chunks, colored according to their B-factors. Each chunk is saved as a separate image (chunk_X.png).
    5) B-Factor Legend: A legend for the B-factor color mapping is generated and saved as legend.png.
    6) Residue Statistics: Residues with non-zero B-factors are grouped by patch ID and plotted in a stacked text format. The output is saved as patch_res.png.
    7) CSV Output: Two CSV files are generated:
        patch.csv: Contains information about residues found in a patch (i.e., with a non-zero B-factor).
        nopatch.csv: Contains information about residues not found in a patch (i.e., with a zero B-factor).
