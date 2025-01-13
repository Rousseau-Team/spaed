<p align="center">
  <img src="img/title.png" border="0"/>
  <h2 align="center">Segmentation of PhAge Endolysin Domains</h2>
</p>

SPAED is a tool to identify domains in phage endolysins. It takes as input the PAE file(s) obtained from AlphaFold and outputs a csv file with delineations.

Additional scripts are provided to visualize predicted domains with PyMOL and to obtain their amino acid sequences. 

## Installation & usage

First create a virtual environment, then: 

**From pypi**:
```
pip install numpy, pandas, scipy, spaed

python
from spaed import spaed
```

ex. `spaed(pae_path, output_file="spaed_predictions.csv")`


**From source**:
```
git clone https://github.com/Rousseau-Team/spaed.git

pip install numpy pandas scipy
```

ex. `python spaed/src/spaed/spaed.py pae_path`

Optional dependency for structure visualisation: pymol (`conda install -c conda-forge -c schrodinger pymol-bundle`)
ex. `python spaed/src/spaed/pymol_vis.py pred_path pdb_path --output_folder pymol_output --output_type {pse|png|both}`

## Advanced usage
**Positional arguments**:
- **pae_path** - Folder of or singular PAE file in json format as outputted by Alphafold2/Colabfold.


**Optional arguments**:
- **output_file** - File to save table of segmented domains in csv format. (default spaed_predictions.csv)
- **fasta_path** - Path to fasta file or folder containing fasta files. If specified, spaed will save the sequences corresponding to predicted domains into a new fasta file named  "spaed_predicted_domains.faa" in the same output folder as output_file. Ensure fasta names or headers correspond to entries in pae files.
- **RATIO_NUM_CLUSTERS** - Maximum number of clusters initially generated by hierarchical clustering corresponds to len(protein) // RATIO_NUM_CLUSTERS. (Default 10). For a protein 400 residues long, 40 clusters will be generated.
- **MIN_DOMAIN_SIZE** - Minimum size a domain can have. (default 30).
- **PAE_SCORE_CUTOFF** - Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions. Residues with PAE score < PAE_SCORE_CUTOFF are considered close together. (default = 4).
- **MIN_DISORDERED_SIZE** - Minimum size a disordered region can be to be considered a separate entity from the domain it is next to (default 20).
- **FREQ_DISORDERED** - For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values <MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted). (default 6).
- **PROP_DISORDERED** - Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict the region as disordered. (default 80%).
- **FREQ_LINKER** - For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered as part of the linker. Values < MIN_DOMAIN_SIZE are logical as they are less than the expected size of the nearest domain. Increasing leads to a more leniant assignment of residues as part of the linker. (default 20).

MAX_CLUSTERS refers to the # of clusters the hierarchical clustering (in scipy) should try to create. A higher number than the number of domains you expect to observe should be used. This is to allow flexible regions (like the ends of the sequences or the linker regions) to be assigned to their own clusters. By default we use one tenth of the length of the protein, but you can manually set this value to any integer.

If you are interested in looking at the disordered regions in N- or C-terminal, consider increasing FREQ_DISORDERED ([4-30]), decreasing MIN_DISORDERED_SIZE ([10-30]) or decreasing PROP_DISORDERED ([50-95]). This will result in more disordered regions being detected, but also many false positives. I would not change them all at the same time as this will probably increase the sensitivity too much.

If you are interested in linkers or have a protein that is less well folded, consider modifying the FREQ_LINKER parameter ([4-30]). This value is used to adjust the boundaries of the linkers and as such, a higher value will result in longer linkers. However, linkers that were missed will still not be detected.

## Citation

For now, please cite this repository when you use SPAED. A preprint will be available shortly.
