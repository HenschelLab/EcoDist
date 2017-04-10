# EcoDist
Ecosystem distribution and entropy of Operational Taxonomic Units as derived from a comprehensive 16S rRNA profile database with Environmental Ontology annotation. OTUs are closed reference from GreenGenes 13.5 with 97% sequence identity. 

## Dependencies: Python 2.7
*matplotlib
*networkx
*MySQLdb

GlobalMicroBiome database as published in Henschel et. al, PLoS Computational Biology, 2015


Some parts also need nltk (possibly that can be avoided)

The central script here is ecoDist.py. Other versions of the script deal with environments where I do not have access to the SQL database.
