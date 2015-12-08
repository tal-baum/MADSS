# README

**MADSS v1.1**
Modular Assembly of Drug Safety Subnetworks
Updated December 7, 2015

This folder contains all of the files necessary to run MADSS. For a description of the algorithm, please see the following publication:

Lorberbaum, T., Nasir, M., Keiser, M., Vilar, S., Hripcsak, G. and Tatonetti, N. (2015), Systems Pharmacology Augments Drug Safety Surveillance. *Clinical Pharmacology & Therapeutics*, 97: 151–158. doi: 10.1002/cpt.2
http://onlinelibrary.wiley.com/doi/10.1002/cpt.2/abstract
First author: tal.lorberbaum_columbia.edu
Corresponding author: nick.tatonetti_columbia.edu

To run MADSS, run the MADSS.py Python script along with the name of the adverse event/ phenotype of interest. For example, to investigate acute myocardial infarction (MI):
`python MADSS.py MI`

Seed sets are available for the following phenotypes:
- Acute myocardial infarction (MI)
- Upper gastrointestinal bleeding (Gastro)
- Acute liver failure (Liver)
- Acute kidney failure (Kidney)
- Long QT Syndrome (LQTS)

For other phenotypes, please specify a seed set using Ensembl protein IDs in the MADSS.py file.

The following Python modules are needed to run MADSS:
networkx (https://networkx.github.io/)
sklearn (http://scikit-learn.org/)
MySQLdb (to query from MySQL databases)


> **Note**: Included in `string700_data.p` and `/stored_vals` are pickled files allowing the user to run MADSS using a pruned PPI network from STRING v9.1 (http://string91.embl.de/, see `/madss_libraries/madss_interactome.py`). In `/stored_vals` we additionally include drug targets from DrugBank v3 for drugs in the acute MI gold standard so the user can generate output from MADSS without needing to connect to an external database (DrugBank, http://www.drugbank.ca/). To investigate other drugs and phenotypes, the user will have to manually compile a list of drug targets or create a local version of the DrugBank database to query.

MADSS is released under a Creative Commons BY-NC-SA 4.0 license. For complete details see LICENSE.txt or visit http://creativecommons.org/licenses/by-nc-sa/4.0/

![CC BY-NC-SA 4.0](https://upload.wikimedia.org/wikipedia/commons/thumb/1/12/Cc-by-nc-sa_icon.svg/100px-Cc-by-nc-sa_icon.svg.png)