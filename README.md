# CLUSTAL : ALIGNEMENT MULTIPLE HEURISTIQUE PAR LA METHODE CLUSTAL
[desctiption du projet]

Objectif : Réalisez un programme reprenant la méthode décrite dans l'article. L’algorithme de construction de l’arbre pourra être remplacé par la construction d’un arbre par embranchement séquentiel. L’algorithme heurisitique d’alignement séquentiel pourra être remplacé par un alignement basé sur la programmation dynamique.

Référence : Desmond G. Higgins , Paul M. Sharp, Fast and sensitive multiple sequence alignments on a microcomputer, Bioinformatics, Volume 5, Issue 2, April 1989, Pages 151–
153, https://doi.org/10.1093/bioinformatics/5.2.151


## Installation
Clone the project:
```
git clone https://github.com/KarinDuong/clustal_alignement.git
cd clustal_alignment
```

Create and activate a conda environment:
```
conda env create -f environment.yml
conda activate clustal_env
```

## Usage as command line tool
Run clustal_alignment on a test file
```
python src/clustal_alignment.py --input [FASTA filename]
```

Run clustal_alignment on a test file, with a different gap score value. By default set at -8.0.
```
python src/clustal_alignment.py --input [FASTA filename] --gap_score [float value]
```

Run clustal_alignment on a test file, with the information of the sequence type in the input file. By default set True for protein sequence.
```
python src/clustal_alignment.py --input [FASTA filename] --gap_score [float value] --protein_sequence [boolean value]
```


