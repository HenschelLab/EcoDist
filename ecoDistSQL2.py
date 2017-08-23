
"""
Creates an ecodistribution plot for a sample (from the database!) similar to qiime barplots plus a second dimension:
for each OTU, the distribution over ecosystems is color coded/visualized
Moreover ecosystem distribution entropy for each OTU in the sample is calculated

for each otu: create a stats of Ecosystem occurrences
SQL version: OTU distribution is calculated with respect to global database (as opposed to Biom table distribution)
This way: faster, more comprehensive, but also more adaptable to otu tables with few or single samples (e.g. Sabkha, Mangroves biom table etc)

Input OTU table: Biom format
Output: ecodistribution plots for each sample in the otu table
TODO: entropy stats: Stronlgy varying entropies = mixture? all high entropies?
TODO: Fix p_i calculation to adjust for amount of samples from that ecosystem

Addition: composite ecosystems: Animal
SELECT compositeEcosystem, COUNT( * ) AS freq FROM  `OTUS_samples_unified` 
NATURAL JOIN CompositeEcosystems WHERE otu_id =  '273400'
GROUP BY compositeEcosystem
--------
compositeEcosystem	freq
Animal/Human	2
Anthropogenic|Plant|Soil	3
Anthropogenic|Soil	41
Hypersaline|Marine	1
Plant	13
Plant|Soil	30
"""
from itertools import chain, groupby
import networkx as nx
import sys,os
from envoTools import Ontology
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
#from selectEnvironments import getChildrenEnvoSQL


def ecodistributionPlot(otudata, width, colors, hatches, names, phyla, phylaNames, entropy, filename):
    from matplotlib.patches import ConnectionPatch
    """colors is an array = len(otudata[0])"""
    otudataN = np.array([row/float(row.sum()) for row in otudata]) ## normalizing row-wise
    sampleEvidence = [row.sum() for row in otudata]
    
    ind = np.hstack(([0], width.cumsum()[:-1]))  #np.arange(len(otudataN))
    left, wid = 0.1, 0.8
    fig = plt.figure(facecolor='white', figsize=(20,20))
    #figlegend = plt.figure(figsize=(12,20))
    ax0 = fig.add_axes([left, 0.89, wid, 0.10]) # sampleEvidence
    ax1 = fig.add_axes([left, 0.78, wid, 0.10], sharex=ax0) # entropy
    ax2 = fig.add_axes([left, 0.38, wid, 0.40], sharex=ax0)  # ecosystem distribution 
    ax3 = fig.add_axes([left, 0.33, wid, 0.04], sharex=ax0) # phylo distribution
    ax4 = fig.add_axes([left, 0.25, wid, 0.04]) # phylo distribution legend
    #ax5 = fig.add_axes([left, 0.0, wid, 0.04]) # phylo distribution legend
    
    bottom = np.zeros(len(otudataN))
    #legendbars = []
    ax0.bar(ind, sampleEvidence, width, color='k', log=True, linewidth=1, fill=False)
    ax0.set_xticks([])
    ax1.bar(ind, entropy, width, linewidth=1, edgecolor='k', fill=False)
    ax1.set_xticks([])
    for idx, habitat in enumerate(otudataN.T):
        color = colors[idx]
        hatch = hatches[idx]
        b = ax2.bar(ind, habitat, width, linewidth=0.1, color=color, bottom=bottom, hatch=hatch)
        #legendbars.append(b[0])
        bottom += habitat
    ax2.set_xticks([])
    #figlegend.legend(legendbars, names)
    ax2.set_ylim(0, 1)
    ind = np.hstack(([0], phyla.cumsum()[:-1]))
    #legendbarsPhyla = []
    for idx, (start, width) in enumerate(zip(ind, phyla)):
        col = phylaColorDict[phylaNames[idx]]
        rect = ax3.bar(start, 1, width, color=col, linewidth=0)
        ax4.bar(idx+0.1, 1, 0.8, color=col, linewidth=0)
        #legendbarsPhyla.append(rect[0])
        con = ConnectionPatch(xyA=(start+width/2.,0), xyB=(idx+0.5, 1), axesA=ax3, axesB=ax4, arrowstyle="->", coordsA="data", coordsB="data", shrinkB=1)
        ax3.add_artist(con)    #ax3.legend(legendbarsPhyla, phylaNames)
    ax3.set_xticks([])
    ax3.set_yticks([])
    for spine in ['right', 'top', 'left', 'bottom']:
        ax3.spines[spine].set_color('none')
        ax4.spines[spine].set_color('none')

    ax4.set_xticks(np.arange(len(phyla)) + 0.5)
    ax4.set_xticklabels(phylaNames, rotation=-90) # [p.split(";")[-1] for p in phylaNames]
    ax4.set_yticks([])    
    #figlegend.legend(legendbarsPhyla, [p.split(";")[-1] for p in phylaNames], 'center')
    fig.savefig("%s.svg"%filename)
    fig.savefig("%s.png"%filename)
    #fig.savefig("%s.pdf"%filename)
    #figlegend.savefig("%s_leg.svg"%filename)
    #figlegend.savefig("%s_leg.png"%filename)
    #figlegend.savefig("%s_leg.pdf"%filename)
    fig.clf()
    #figlegend.clf()
    plt.close(fig)
    #plt.close(figlegend)
                               
def clean(p):
    return p.strip().strip('"').split("__")[-1]

def extractPhylum(lineage, level=3):
    lineage = ";".join([clean(p) for p in lineage.split(";")])
    return ";".join(lineage.split(";")[:level])

class Sample:
    def __init__(self, sampleID='Sabkha'):
        ## this code needs to be optimized to be run on the server!
        self.sampleID = sampleID
        try:
            curs.execute('SELECT otu_id, sequence_count FROM OTUS_samples_unified WHERE sample_event_id = "%s"' % self.sampleID)   ## REMOVE LIMIT!!!
            res = curs.fetchall()
            otus = [OTU(rec["otu_id"], rec["sequence_count"]) for rec in res]
            #otusMG = [OTU(rec["otu_id"], rec["sequence_count"], includeMGRast=True) for rec in res]
        except:
            pass
            #
            pdb.set_trace()
        
        self.otus = sorted([otu for otu in otus if otu.hasData])

        if True: ## 
            if len(otus) < 5: return ## don't do anything for small samples!
            otudata = [otu.ecoDistribution for otu in self.otus]

            widths = np.array([otu.count for otu in self.otus])
            phylaNames, phylaCount = zip(*[(k, sum([e.count for e in g])) for k, g in groupby(self.otus, lambda o:o.phylum)])

            entropies =  [otu.entropy for otu in self.otus]
            wEntropies =  np.array([(otu.entropy, otu.count) for otu in self.otus])
            wEntropyAvg = np.dot(wEntropies[:,1],wEntropies[:,0])/wEntropies[:,1].sum()
            print sampleID, np.array(entropies).mean(), wEntropyAvg
            try:
                curs.execute("INSERT INTO  samples_EcoDistributionEntropy  VALUES ('%s',  '%s',  '%s')" % (sampleID, np.array(entropies).mean(), wEntropyAvg))
            except:
                pass
                #pdb.set_trace()
            ecodistributionPlot(otudata, widths, cols4plot,hatch4plot, ecosystems + ['Misc'], np.array(phylaCount), phylaNames, entropies, "%s/%s" % (resultdir, sampleID))
   
class OTU:
    def __init__(self, otuID, count, includeMGRast=False):
        self.id = otuID
        self.count = count
        self.hasData = False
        if otuID == '':
            print "Warning: empty otuID"
            return
        query = "SELECT * FROM OTUS_unified WHERE otu_id=%s" % otuID
        curs.execute(query)
        result = curs.fetchone()
        lineage = result["lineage"] if result else "?"
        self.lineage = ";".join([clean(p) for p in lineage.split(";")])
        self.phylum = ";".join(self.lineage.split(";")[:3])
        self.addEcoDistributionSQL()
        self.calculateEntropy()
        
    def addEcoDistributionSQL(self):
        #query = "SELECT compositeEcosystem, COUNT(*) AS freq FROM `OTUS_samples_unified` NATURAL JOIN CompositeEcosystems WHERE otu_id='%s' GROUP BY compositeEcosystem" % self.id
        #curs.execute(query)
        #recordsMG = curs.fetchall()
        #print "\n".join(["%-40s: %7d"%(e,i) for i,e in sorted([rec.values() for rec in recordsMG])])
        #print "---------------------------------------------------------------------------------------"
        ## TODO: normalize by # of samples in (composite) ecosystem
        queryOld = "SELECT compositeEcosystem, COUNT(*) AS freq FROM `OTUS_samples_unified` NATURAL JOIN samples_unified NATURAL JOIN CompositeEcosystems_bak WHERE otu_id='%s' AND NOT study LIKE 'MG%%' GROUP BY compositeEcosystem" % self.id ## Without MG-Rast datasets
        curs.execute(queryOld)
        self.ecoDistribution = np.zeros(len(ecosystems))
        records = curs.fetchall()
        #print "\n".join(["%-40s: %7d"%(e,i)for i,e in sorted([rec.values() for rec in records])])
        #print "#######################################################################################"

        for rec in records:            
            position = ecosystemsIndex[rec["compositeEcosystem"]]
            self.ecoDistribution[position] = rec["freq"]
        
        self.hasData = self.ecoDistribution.sum() > 0

    def calculateEntropy(self):
        ## this needs a fix!
        ## Probabilities need to be calculated according to how skewed the dataset is!
        from scipy.stats.distributions import entropy
        self.entropy = entropy(self.ecoDistribution/self.ecoDistribution.sum())
    def dumpData(self,dir):
        dump.dump((self.lineage, self.count, self.phylum, self.ecoDistribution, self.entropy), "%s/otu_%s.pcl"%(dir, self.id))        
    def __cmp__(self, o):
        return cmp(self.lineage, o.lineage)
             
if __name__ == "__main__":
    datadir = "data"
    resultdir = sys.argv[2]
    
    ## MySQL connection
    conn = MySQLdb.connect(db="EarthMicroBiome", host="research04", user="ahenschel", passwd="angi4rf")
    curs = conn.cursor(DictCursor)
    
    ## Settings
    ecosystemColors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal/Human': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Air': 'AliceBlue', 'Hypersaline':'r'}
    hatchpatterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ]*4 ## for composite ecosystems
    curs.execute("SELECT DISTINCT compositeEcosystem FROM CompositeEcosystems_bak ORDER BY  CompositeEcosystem DESC") ## CHANGE?
    ecosystems = [rec["compositeEcosystem"] for rec in curs.fetchall()]
    ecosystemsIndex  = dict([(eco,idx) for idx,eco in enumerate(ecosystems)])
    ecosystemsIndex1 = dict([(eco,idx) for idx,eco in enumerate(ecosystemColors.keys())])

    try:
        reuse = True
        if reuse:
            phylaColorDict = dump.load("%s/phylaColorDict2.pcl"%datadir)
            phylaColorDict['?'] = '#f0f0f0'
        else:
            print "Generating new color dictionary for phyla..."
            curs.execute("SELECT DISTINCT(lineage) FROM OTUS_unified")
            usedPhyla = set([extractPhylum(rec["lineage"]) for rec in curs.fetchall()])
            phylaColorDict = dict([(phylum, htmlcols[np.random.randint(len(htmlcols))]) for phylum in usedPhyla])
            dump.dump(phylaColorDict, "/home/zain/Projects/KnowYourEnv/Data/phylaColorDict2.pcl")
            
        cols4plot =  [ecosystemColors[ecosystem.split("|")[0]] for ecosystem in ecosystems] + ["b"]
        eco2hatch = lambda ecosystem: "".join([hatchpatterns[ecosystemsIndex1[eco]] for eco in ecosystem.split("|")[1:]])
        hatch4plot = map(eco2hatch, ecosystems)

        ## sample IDs
        sediments = ['Sabkha', 'Soil.Day.0', 'E2.Roots', 'E2.control,DNA']
        extremeEnvironments = ['P20.C.filt..660379', 'WPC.sed.D1.660400', 'P20.B.filt..660401', 'A23.2.sed.D1.660392', 'A23.number1.filt.D1.660399', 'WPA.filt..660391', 'IM4z.609442', 'SA4x.609398', 'P.Masambaba.SB.414876', 'P.Masambaba.SA.414862']
        sampleIDs = sediments[1:2]# + extremeEnvironments # [sediments.index('19.141761'):]
        
        for sampleID in sampleIDs:
            print " ######### %s ############### " % sampleID
            if not os.path.exists("%s/%s.pdf" %(resultdir,sampleID)):
                print "Processing", sampleID
                sample = Sample(sampleID) ## does all the work!
            else:
                print "%s exists, skipping" % sampleID

    finally:
        print "Closing MySQL connection"
        conn.close()
