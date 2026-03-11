# Paving the way for advanced genomic resources in the European red algae Porphyra dioica and Porphyra linearis

This repository includes all the scripts used to generate chromosome-level genome assemblies and high-quality genomic resources for the species *Porphyra dioica* and *Porphyra linearis*, facilitating comparative genomics and evolutionary studies within the Biodiversity Genomics Europe (BGE) framework.

# Authors 
Jordi Morcillo-Baeza, Phycology Research Group, Ghent University, Krijgslaan 281 S8, 9000 Ghent, Belgium

Jessica Knoop, Phycology Research Group, Ghent University, Krijgslaan 281 S8, 9000 Ghent, Belgium

Olivier De Clerck, Phycology Research Group, Ghent University, Krijgslaan 281 S8, 9000 Ghent, Belgium

# Software Versions Used
| Tool | Version | Purpose |
| :--- | :--- | :--- |
| Flye | 2.9.6 | Long-read de novo assembly |
| Racon | 1.5.0 | Long-read polishing (3 rounds) |
| Medaka | 1.12.1 | Consensus polishing (Model: r1041_e82_400bps_sup_v5.0.0) |
| fastp | 0.23.4 | Hi-C adapter trimming |
| HapHiC | 1.0.1 | Hi-C Scaffolding |
| Juicebox | 1.8.8 | Manual curation |
| Tiara | 1.0.2 | Taxonomic classification |
| SeqKit | 2.8.2 | Coverage and assembly statistics |
| Quast | 5.0.2 | Assembly quality assessment |
| BUSCO | 5.8.2 | Assembly quality assessment |

# Acknowledgements
This work was supported by the SIB Swiss Institute of Bioinformatics under the Biodiversity Genomics Europe project.

