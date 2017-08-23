
"""
Creates an ecodistribution plot for a sample similar to qiime barplots plus a second dimension:
for each OTU, the distribution over ecosystems is color coded/visualized
Moreover ecosystem distribution entropy for each OTU in the sample is calculated

for each otu: create a stats of Ecosystem occurrences
SQL version: OTU distribution is calculated with respect to global database (as opposed to Biom table distribution)
This way: faster, more comprehensive, but also more adaptable to otu tables with few or single samples (e.g. Sabkha, Mangroves biom table etc)

Input OTU table: Biom format
Output: ecodistribution plots for each sample in the otu table
TODO: entropy stats: Stronlgy varying entropies = mixture? all high entropies?
"""

import numpy as np
import dump
import pdb
import MySQLdb
from MySQLdb.cursors import DictCursor

def clean(p):
    return p.strip().strip('"').split("__")[-1]

class OTU:
    def __init__(self, otuID, count):
        self.id = otuID
        self.count = count
        query = "SELECT * FROM OTUS_unified_NR WHERE otu_id=%s" % otuID
        curs.execute(query)
        result = curs.fetchone()
        lineage = result["lineage"] if result else "?"
        self.lineage = ";".join([clean(p) for p in lineage.split(";")])
        self.phylum = ";".join(self.lineage.split(";")[:3])
        query = "SELECT ecosystem, COUNT(*) AS freq FROM `OTUS_samples_unified` NATURAL JOIN samples_EnvO_annotation_unified NATURAL JOIN envoColors WHERE otu_id='%s' GROUP BY ecosystem" % self.id
        curs.execute(query)
        self.ecoDistribution = np.zeros(len(ecosystems))
        for rec in curs.fetchall():
            position = ecosystemsIndex[rec["ecosystem"]]
            self.ecoDistribution[position] = rec["freq"]
        self.hasData = self.ecoDistribution.sum() > 0

        from scipy.stats.distributions import entropy
        self.entropy = entropy(self.ecoDistribution/self.ecoDistribution.sum())

    def dumpData(self,dir):
        dump.dump((self.lineage, self.count, self.phylum, self.ecoDistribution, self.entropy), "%s/otu_%s.pcl"%(dir, self.id))        
    def __cmp__(self, o):
        return cmp(self.lineage, o.lineage)

if __name__ == "__main__":
    ## Settings
    datadir = "/data/EarthMicrobiomeProject/OTUdata"

    ## MySQL connection
    conn = MySQLdb.connect(db="GlobalMicroBiome", host="cis1-db", user="ahenschel", passwd="angi4rf")
    curs = conn.cursor(DictCursor)
    ecosystems = ['Air', 'Plant', 'Hypersaline', 'Soil', 'Animal/Human', 'Freshwater', 'Marine', 'Geothermal', 'Anthropogenic', 'Biofilm']
    ecosystemsIndex = dict([(eco,idx) for idx,eco in enumerate(ecosystems)])
    try:
        curs.execute("SELECT DISTINCT(otu_id) FROM `OTUS_unified_NR`")
        for rec in curs.fetchall(): #allOTUs():            
            otu = OTU(rec["otu_id"], 0)
            otu.dumpData(datadir)
    finally:
        conn.close()
