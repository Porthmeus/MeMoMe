# MeMoMe Databases

# Purpose
The databases in this repository are more or less direct copies of the original databases, only standardized for easier handling in MeMoMe. In general, all files contain the following columns:

| Column | Description | Example (HMDB) |
|:---|:-----|---:|
| id | Contains the accession id in the original data base. This is a single value for each metabolite | "HMDB0000158" |
|names| contains the names, synonyms and aliases (including IUPAC names, if given) for the compound. The data is a single string seperating each individual name by \_\|\_.| "L-Tyrosine\_\|\_Tyr\_\|\_(2S)-2-amino-3-(4-hydroxyphenyl)propanoic acid |
| annotation | Contains the accession ids to other metabolite databases which are cross referenced formatted as python dictionary. Note that this column contains only entries which are registered with identifiers.org and have a relevant id. However, some special cases exists like metanetx. | "{'drugbank':['DB00135'], 'biocyc':['TYR'], 'chebi':['17895']} |
| inchi | Contains the (standardized) inchi string | "InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4-2-6/h1-4,8,11H,5,10H2,(H,12,13)/t8-/m0/s1" |

The rest of the database is either saved as separate columns following the described ones or it is removed to keep the file reasonably small (e.g. for HMDB, most information is removed). Sometimes, additional information is derived from external databases. For example, BiGG provides InChI-Keys, but no InChI-strings for the metabolites in the database, here we used the pubchem database to infer the strings from the keys.

# Licensing
The data in general is licenced under the GLP3.0. However, since most of the data is not new, the original licence applies.

|Database|License|
|:-------|:------|
| BiGG | BiGG-License |
| VMH | CC-BY2.0 |
| HMDB | Custom-License, but free for non-commercial use |
| gapseq | GLP3.0 |
| ModelSeed | CC-BY2.0 |

