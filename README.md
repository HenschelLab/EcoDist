# EcoDist
Ecosystem distribution and entropy of Operational Taxonomic Units as derived from a comprehensive 16S rRNA profile database with Environmental Ontology annotation. OTUs are closed reference from GreenGenes 13.5 with 97% sequence identity. 

## Dependencies: Python 2.7

* matplotlib
* networkx
* MySQLdb
* scipy
* numpy
* biom

SQL MicroBiome database as published in [Henschel et. al, PLoS Computational Biology, 2015](http://dx.doi.org/10.1371/journal.pcbi.1004468)

Some parts also need nltk (possibly that can be avoided)

The central script here is ecoDist.py.
It contains a
It uses precalculated ecosystem distribution information per OTU from data/usedOTUs.pcl

ecoDistSQL*.py deal with SQL info, incl. dumping
