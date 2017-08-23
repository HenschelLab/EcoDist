
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
from itertools import chain, groupby
import networkx as nx
import sys,os
import numpy as np
from scipy.spatial.distance import squareform
from biom.parse import parse_biom_table
from collections import defaultdict
import matplotlib.pyplot as plt
from htmlcols import htmlcols
import dump
import pdb
import MySQLdb
from MySQLdb.cursors import DictCursor
import time

def ecodistributionPlot(otudata, width, colors, names, phyla, phylaNames, entropy, filename):
    from matplotlib.patches import ConnectionPatch
    """colors is an array = len(otudata[0])"""
    otudataN = np.array([row/float(row.sum()) for row in otudata]) ## normalizing row-wise

    ind = np.hstack(([0], width.cumsum()[:-1]))  #np.arange(len(otudataN))
    left, wid = 0.1, 0.8
    fig = plt.figure(facecolor='white', figsize=(20,20))
    figlegend = plt.figure()
    ax1 = fig.add_axes([left, 0.78, wid, 0.15]) # entropy
    ax2 = fig.add_axes([left, 0.28, wid, 0.5], sharex=ax1)  # ecosystem distribution 
    ax3 = fig.add_axes([left, 0.23, wid, 0.04], sharex=ax1) # phylo distribution
    ax4 = fig.add_axes([left, 0.15, wid, 0.04]) # phylo distribution legend
    #ax5 = fig.add_axes([left, 0.0, wid, 0.04]) # phylo distribution legend
    
    bottom = np.zeros(len(otudataN))
    legendbars = []
    ax1.bar(ind, entropy, width, linewidth=0)
    ax1.set_xticks([])
    for idx, habitat in enumerate(otudataN.T):
        color = colors[idx]
        b = ax2.bar(ind, habitat, width, linewidth=0.1, color=color, bottom=bottom)        
        legendbars.append(b[0])
        bottom += habitat
    ax2.set_xticks([])
    #ax2.legend(legendbars, names)
    ax2.set_ylim(0, 1)
    ind = np.hstack(([0], phyla.cumsum()[:-1]))
    legendbarsPhyla = []
    for idx, (start, width) in enumerate(zip(ind, phyla)):
        col = phylaColorDict[phylaNames[idx]]
        rect = ax3.bar(start, 1, width, color=col, linewidth=0)
        ax4.bar(idx+0.1, 1, 0.8, color=col, linewidth=0)
        legendbarsPhyla.append(rect[0])
        con = ConnectionPatch(xyA=(start+width/2.,0), xyB=(idx+0.5, 1), axesA=ax3, axesB=ax4, arrowstyle="->", coordsA="data", coordsB="data", shrinkB=1)
        ax3.add_artist(con)    #ax3.legend(legendbarsPhyla, phylaNames)
    ax3.set_xticks([])
    ax3.set_yticks([])
    for spine in ['right', 'top', 'left', 'bottom']:
        ax3.spines[spine].set_color('none')
        ax4.spines[spine].set_color('none')

    ax4.set_xticks(np.arange(len(phyla)) + 0.5)
    ax4.set_xticklabels([p.split(";")[-1] for p in phylaNames], rotation=-90)
    ax4.set_yticks([])
    figlegend.legend(legendbarsPhyla, [p.split(";")[-1] for p in phylaNames], 'center')
    fig.savefig("%s.svg"%filename)
    fig.savefig("%s.png"%filename)
    figlegend.savefig("%s_leg.svg"%filename)
    figlegend.savefig("%s_leg.png"%filename)
    fig.clf()
    figlegend.clf()
    plt.close(fig)
    plt.close(figlegend)
                               
def clean(p):
    return p.strip().strip('"').split("__")[-1]

def extractPhylum(lineage, level=3):
    lineage = ";".join([clean(p) for p in lineage.split(";")])
    return ";".join(lineage.split(";")[:level])


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
        t1 = time.time()
        self.addEcoDistributionSQL()
        t2 = time.time()
        self.calculateEntropy()
    def addEcoDistributionSQL(self):
        query = "SELECT ecosystem, COUNT(*) AS freq FROM `OTUS_samples_unified` NATURAL JOIN samples_EnvO_annotation_unified NATURAL JOIN envoColors_unified WHERE otu_id='%s' GROUP BY ecosystem" % self.id
        curs.execute(query)
        self.ecoDistribution = np.zeros(len(ecosystems))
        for rec in curs.fetchall():
            position = ecosystemsIndex[rec["ecosystem"]]
            self.ecoDistribution[position] = rec["freq"]
        self.hasData = self.ecoDistribution.sum() > 0

    def calculateEntropy(self):
        from scipy.stats.distributions import entropy
        self.entropy = entropy(self.ecoDistribution/self.ecoDistribution.sum())
        #curs.execute("UPDATE OTUS_unified_NR SET entropy =  '%s' WHERE  otu_id =  '%s';" % (self.entropy, self.id))
        #normed = self.ecoDistribution/totalEcoDistribution
        #self.entropyNormalized = entropy(normed/normed.sum())
    def dumpData(self,dir):
        dump.dump((self.lineage, self.count, self.phylum, self.ecoDistribution, self.entropy), "%s/otu_%s.pcl"%(dir, self.id))        
    def __cmp__(self, o):
        return cmp(self.lineage, o.lineage)
    
def getOTU(idx, sampleData):
    """lookup in cache, if OTU not there, create it and put in otucache"""
    otuid = otuTable.ObservationIds[idx]
    if otucache.has_key(otuid):
        return otucache[otuid]
    otu = OTU(otuid, sampleData[idx])
    otucache[otuid] = otu
    return otu

if __name__ == "__main__":
    ## Settings
    ecosystemColors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal/Human': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Air': 'AliceBlue', 'Hypersaline':'r'}
    ecosystems = ecosystemColors.keys()
    ecosystemsIndex = dict([(eco,idx) for idx,eco in enumerate(ecosystems)])
    datadir = "/home/zain/Projects/KnowYourEnv/Data/"
    resultdir = "/home/zain/Projects/KnowYourEnv/Results/EcoPhylPlotsAll"

    ## MySQL connection
    conn = MySQLdb.connect(db="GlobalMicroBiome", host="cis1-db", user="ahenschel", passwd="angi4rf")
    curs = conn.cursor(DictCursor)

    try:
        reuse = True
        if reuse:
            phylaColorDict = dump.load("/home/zain/Projects/KnowYourEnv/Data/phylaColorDict.pcl")
            otucache = dump.load("/home/zain/Projects/KnowYourEnv/Data/otucache.pcl")
        else:
            curs.execute("SELECT DISTINCT(lineage) FROM OTUS_unified_NR")
            usedPhyla = set([extractPhylum(rec["lineage"]) for rec in curs.fetchall()])
            phylaColorDict = dict([(phylum, htmlcols[np.random.randint(len(htmlcols))]) for phylum in usedPhyla])
            otucache = {}
        phylaColorDict['?'] = '#f0f0f0'
        cols4plot = [ecosystemColors[ecosystem] for ecosystem in ecosystems] + ['b']
        ## Biom table
        #biom = "/home/handreas/Data/SequencesMITDec2013/FLASH_merged/SplittedLibraries_Corrected/Areej/otu_table.biom
        for biom in [#"/home/handreas/Data/SequencesMITDec2013/FLASH_merged/SplittedLibraries_Corrected/Areej/otu_table_Day0andE2Soil.biom",
                     "/home/handreas/Data/SequencesMITDec2013/FLASH_merged/SplittedLibraries_Corrected/Areej/otutable_crabgut.biom",
                     "/home/handreas/SHAIMA/Final/otu_table.biom",
                     "/data/EarthMicrobiomeProject/BetaDiversity/All_against_all_update_r2000.biom"]: ## Note, this is a smaller table!
            otuTable = parse_biom_table(open(biom, 'U'))
            print "Parsed otu table %s"%biom

            mincount = 0
            for sampleData, sampleName, m in otuTable.iterSamples():
                if True:
                    t1 = time.time()
                    print "Creating %s OTU objects for sample %s " %(len( np.where(sampleData > mincount)[0]), sampleName)
                    selectedOTUs = [getOTU(idx, sampleData) for idx in np.where(sampleData > mincount)[0]]

                    selectedOTUs = [otu for otu in selectedOTUs if otu.hasData]
                    t2 = time.time()
                    print sampleName, t2-t1
                    selectedOTUs.sort()
                    otudata = [otu.ecoDistribution for otu in selectedOTUs]
                    widths = np.array([otu.count for otu in selectedOTUs]) #np.ones(len(otudata)) * 1 ## possibly stretch by occ.
                    phylaNames, phylaCount = zip(*[(k, sum([e.count for e in g])) for k, g in groupby(selectedOTUs, lambda o:o.phylum)])
                    entropies =  [otu.entropy for otu in selectedOTUs]
                    ecodistributionPlot(otudata, widths, cols4plot, ecosystems + ['Misc'], np.array(phylaCount), phylaNames, entropies, "%s/%s" % (resultdir, sampleName))

        #        except:
        #            print "Error in", sampleName
    finally:
        print "Closing MySQL connection"
        conn.close()
