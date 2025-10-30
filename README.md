Analysis of macromolecular complexes data in the PDB
==

## Disclaimer
This is a repository of a [Jupyter notebook](https://github.com/PDBe-KB/pdbe-assemblies-analysis/blob/main/assemblies_analysis.ipynb), which is supplementary material for the publication:

**Publication**<br>
Appasamy, S.D., Berrisford, J., Gaborova, R., Nair, S., Anyango, S., Grudinin, S., Deshpande, M., Armstrong, D., Pidruchna, I., Ellaway, J.I.J., Leines, G.D., Gupta, D., Harrus, D., Varadi, M. and Velankar, S. (2023) Annotating macromolecular complexes in the Protein Data Bank: improving the FAIRness of structure data. Scientific Data, 10, 853. https://doi.org/10.1038/s41597-023-02778-9
.

The code in this notebook reproduces the analysis presented in the publication.

> **Update (October 2025):**  
> The data CSV files in the `data/` directory have been updated to align with the analyses prepared for the **PDBe-KB Complexes** draft manuscript.  

## Google Colab
An interactive version of this Jupyter notebook is [available here](https://colab.research.google.com/github/PDBe-KB/pdbe-assemblies-analysis/blob/main/assemblies_analysis.ipynb).

## Background
Macromolecular complexes are crucial functional units in virtually all cellular processes. Their atomic-level understanding is vital to understanding molecular mechanisms and affects applications, such as developing new therapeutics. The Protein Data Bank (PDB) is the central repository for experimentally determined macromolecular structures. However, it can be challenging to find all instances that represent the same assembly in the PDB due to the current PDB annotation practices, which do not include the annotation of assemblies. This study highlights the importance of annotations for macromolecular complexes and the need for more robust methods to identify and classify these complexes across the PDB. We propose a new approach that uses external resources such as the Complex Portal and Gene Ontology to describe assemblies accurately and put them into their biological contexts. We anticipate that the new approach to identifying and classifying complexes will improve the usability and utility of the PDB for researchers in the field of structural biology and drug discovery.
