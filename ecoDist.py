
"""
for each otu: create a stats of Ecosystem occurrences

"""
nomysql = False
try:
    import MySQLdb
    from MySQLdb.cursors import DictCursor
except ImportError, ie:
    nomysql = True

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
import time, pdb

def ecodistributionPlot(otudata, width, colors, names, phyla, phylaNames, entropy, filename):
    """colors is an array = len(otudata[0])"""
    otudataN = np.array([row/float(row.sum()) for row in otudata]) ## normalizing row-wise

    ind = np.hstack(([0], width.cumsum()[:-1]))  #np.arange(len(otudataN))
    left, wid = 0.1, 0.8
    fig = plt.figure(facecolor='white')
    figlegend = plt.figure()
    ax1 = fig.add_axes([left, 0.78, wid, 0.15])
    ax2 = fig.add_axes([left, 0.25, wid, 0.5])
    ax3 = fig.add_axes([left, 0.08, wid, 0.15])
    
    bottom = np.zeros(len(otudataN))
    legendbars = []
    ax1.bar(range(len(entropy)), entropy, linewidth=0)
    ax1.set_xticks([])
    for idx, habitat in enumerate(otudataN.T):
        color = colors[idx]
        b = ax2.bar(ind, habitat, width, linewidth=0, color=color, bottom=bottom)        
        #legendbars.append(b[0])
        bottom += habitat
    ax2.set_xticks([])
    #ax2.legend(legendbars, names)
    ind = np.hstack(([0], phyla.cumsum()[:-1]))
    legendbarsPhyla = []
    for idx, (start, width) in enumerate(zip(ind, phyla)):
        col = phylaColorDict[phylaNames[idx]] 
        rect = ax3.bar(start, 1, width, color=col, linewidth=0)        
        legendbarsPhyla.append(rect[0])
    #ax3.legend(legendbarsPhyla, phylaNames)
    ax3.set_xticks([])
    figlegend.legend(legendbarsPhyla, phylaNames, 'center')
    fig.savefig("%s.svg"%filename)
    fig.savefig("%s.png"%filename)
    figlegend.savefig("%s_leg.svg"%filename)
    figlegend.savefig("%s_leg.png"%filename)
    fig.clf()
    figlegend.clf()
    plt.close(fig)
    plt.close(figlegend)
                           
    
def lookupOntologyClass(sampleID): ## TODO, also implement more abstract levels
    query = "SELECT OntologyID FROM samples_EnvO_annotation_unified WHERE sample_event_ID='%s' AND ontology='envo'"%(sampleID)
    curs.execute(query)
    result = curs.fetchall()
    if not result:
        print >> sys.stderr, "Warning: no ontology annotation for sample", sampleID
        return []
    return [rec["OntologyID"] for rec in result]
    
def lookupSampleInfo(sampleID): ## TODO, also implement more abstract levels
    query = "SELECT * FROM samples_unified WHERE sample_event_ID='%s'"%sampleID    
    curs.execute(query)
    return curs.fetchone()
    
def addEnvoSampleAssociation(envoID, sampleID):
    envo2Sample[envoID].append(sampleID)
    for parent in envo.G.successors(envoID):
        addEnvoSampleAssociation(parent, sampleID)
        
def mapEcosystems2samples():
    """
    produces a dict:
    {'Biofilm': [sampleID1, sampleID2 ...], 'Soil': [sampleID1, ...]}
    """
    ecosystems2samples = defaultdict(list)
    for idx, ecosystem in enumerate(ecosystems):
        obsEnvosInEcosystem = set()
        for ecosystemCat in categories[ecosystem]:            
            if envo.synonymDict.has_key(ecosystemCat):
                envoClass = envo.synonymDict[ecosystemCat][0]
                subsumedEnvos = [en for en in envo2Sample.keys() if nx.has_path(envo.G, en, envoClass)]
                obsEnvosInEcosystem = obsEnvosInEcosystem.union(subsumedEnvos)
            else:
                print >> sys.stderr, "Warning: ignoring toplevel ecosystem branch %s " % ecosystemCat        
        ecosystems2samples[ecosystem] = list(set(chain.from_iterable([envo2Sample[e] for e in obsEnvosInEcosystem])))
    return ecosystems2samples

def clean(p):
    return p.strip().strip('"').split("__")[-1]
class OTU:
    def __init__(self, otuID):
        self.id = otuID
        #self.distribution = otuDistribution
        query = "SELECT * FROM OTUS_unified_NR WHERE otu_id=%s" % otuID
        curs.execute(query)
        result = curs.fetchone()
        lineage = result["lineage"] if result else "?"
        self.lineage = ";".join([clean(p) for p in lineage.split(";")])
        self.phylum = ";".join(self.lineage.split(";")[:3])
    def addEcoDistribution(self,ecodist):
        self.ecoDistribution = np.array(ecodist)
    def calculateEntropy(self, totalEcoDistribution):
        from scipy.stats.distributions import entropy
        self.entropy = entropy(self.ecoDistribution/self.ecoDistribution.sum())
        #normed = self.ecoDistribution/totalEcoDistribution
        #self.entropyNormalized = entropy(normed/normed.sum())
    def __cmp__(self, o):
        return cmp(self.lineage, o.lineage)
    
if __name__ == "__main__":
    t0 = time.time()
    categories = {'Biofilm': ["biofilm", "microbial mat", "biofilm material", "microbial mat material"],
                  'Anthropogenic': ['waste water', 'contamination feature', 'bioreactor', 'biofilter', 'anthropogenic environmental material',
                                    'anthropogenic habitat', 'anthropogenic abiotic mesoscopic feature', 'anthropogenic feature'],
                  'Marine': ['marine feature', 'marine sediment', 'marine biome', 'saline water habitat',
                             "archipelago", 'reef', 'seashore', 'saline hydrographic feature', 'marine snow',
                             'undersea feature', 'coast', 'saline marsh', 'mid-ocean ridge', 'hydrothermal vent', 'saline water'],
                  'Freshwater': ['fresh water', 'freshwater biome', 'freshwater habitat', 'snow field', 'glacier', 'iceberg', 'ice mass', 'glacial feature', 'lake shore', 'crater lake', 'lake sediment', 'mesotrophic water', 'aquifer'],
                  'Soil': ['beach', 'desert', 'karst', 'landslide', 'terrace', 'soil', 'clay', 'sandy sediment', 'mud',
                           'pebble sediment', 'subterrestrial habitat', 'terrestrial habitat', 'sand', 'plantation'],
                  'Animal': ['animal-associated habitat', 'bodily fluid', 'animal food product'],
                  'Plant': ['plantation', 'plant-associated habitat', 'plant food product'],
                  'Geothermal': ['volcanic feature', 'volcano'],
                  'Hypersaline': ['hypersaline']                           
                  }

    colors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Air': 'AliceBlue', 'Hypersaline':'r'}
    ecosystems = categories.keys()

    conn = MySQLdb.connect(db="GlobalMicroBiome", host="cis1-db", user="ahenschel", passwd=passwd)
    curs = conn.cursor(DictCursor)

    # Load ontology (see EnvO tools), OttoTextMining
    ontoDir = "OntologyData"
    envo = Ontology("%s/envoTerms3.pcl" % ontoDir, "%s/envo3.pcl" % ontoDir, "envo")
    
    # Load Biom table
    t0a = time.time()
    print "CheckPoint 1: %8.3f sec"%(t0a-t0)

    biom = "/data/EarthMicrobiomeProject/BetaDiversity/All_against_all_update_r2000.biom" ## Note, this is a smaller table!

    datadir = "/home/zain/Projects/KnowYourEnv/Data/"
    resultdir = "/home/zain/Projects/KnowYourEnv/Results/EcoPhylPlots"
    #pdb.set_trace()
    otuTable = parse_biom_table(open(biom, 'U')) 
    #totalOTUdistribution = otuTable.sum(axis="sample") # should not change for rarefied OTU tables  # very slow, btw!
    t1 = time.time()
    print "CheckPoint 2: %8.3f sec"%(t1-t0a)
    ## for each sample: find membership to ontology classes, including parent ontology classes (includes all ancestral classes)
    ## seems to be slow - precalc envo2Samples, ecosystem2samples
    ## Preparing dictionaries: envo2Sample       {<envoID>: [sample1, sample2 ...], ...}
    ##                         ecosystems2Samples {'marine': [sample1, sample2 ...], ...}
    precalc = True
    if not precalc:
        sampleIDs = otuTable.SampleIds
        envo2Sample = defaultdict(list)
        for sampleIDidx, (sd, sampleID, m) in enumerate(otuTable.iterSamples()):
            for envoID in lookupOntologyClass(sampleID):
                addEnvoSampleAssociation(envoID, sampleIDidx) ## can lead to redundancies! Therefore see below 2 lines
        for en, samples in envo2Sample.items():
            envo2Sample[en] = list(set(samples))
        sampleIDs = np.array(sampleIDs) ## allows for fancy indexing
        t1a = time.time()
        print "CheckPoint 2a: %8.3f sec"%(t1a-t1)

        ecosystems2samples = mapEcosystems2samples()
        dump.dump(envo2Sample, "%s/envo2samples.pcl"%datadir)
        dump.dump(ecosystems2samples, "%s/ecosystems2samples.pcl"%datadir)
    else:
        envo2samples = dump.load("%s/envo2samples.pcl"%datadir)
        ecosystems2samples = dump.load("%s/ecosystems2samples.pcl"%datadir)
        pass
    cols4plot = [colors[ecosystem] for ecosystem in ecosystems] + ['b']
    t2 = time.time()
    print "CheckPoint 3: %8.3f sec"%(t2-t1)

    if not precalc:
        usedOTUs = {}
        for otuDistribution, otuID, meta in otuTable.iterObservations():
            #if otuID not in sabkhaOTUs: continue 
            if otuDistribution.sum() < 100 and len(np.where(otuDistribution>0)[0]) < 100: continue ## change this!
            otu = OTU(otuID) ## needs cis1-db ## also get nr or Sequences        
            ecodist = [otuDistribution[ecosystems2samples[ecosystem]].sum() for ecosystem in ecosystems] ##TODO account for overlaps
            ecodist.append(max(otuDistribution.sum()-sum(ecodist),0)) ## account for Misc., this is approx. as some could be accounted for twice...
            otu.addEcoDistribution(ecodist)
            otu.calculateEntropy(0)
            usedOTUs[otuID] = otu
    else:
        usedOTUs = dump.load("%s/usedOTUs.pcl"%datadir)

    ## generate color dict for used phyla once, so to be consistent in subsequent plots
    t3 = time.time()
    print "CheckPoint 4: %8.3f sec"%(t3-t2)
    usedPhyla = set([otu.phylum for otu in usedOTUs.values()])
    phylaColorDict = dict([(phylum, htmlcols[np.random.randint(len(htmlcols))]) for phylum in usedPhyla])    

    for sampleData, sampleName, m in otuTable.iterSamples():
        #if os.path.exists("%s/%s.png" % (resultdir, sampleName)): continue
        try:
            selectedOTUs = [usedOTUs[otuTable.ObservationIds[idx]] for idx in np.where(sampleData>0)[0] if usedOTUs.has_key(otuTable.ObservationIds[idx])]
            selectedOTUs.sort()
            otudata = [otu.ecoDistribution for otu in selectedOTUs]
            widths = np.ones(len(otudata)) * 1 ## possibly stretch by occ.        
            phyla0 = [otu.phylum for otu in selectedOTUs]
            phylaNames, phylaCount = zip(*[(i,len(list(l))) for i,l in groupby(phyla0)])
            entropies =  [otu.entropy for otu in selectedOTUs]
            ecodistributionPlot(otudata, widths, cols4plot, ecosystems + ['Misc'], np.array(phylaCount), phylaNames, entropies, "%s/%s" % (resultdir, sampleName))
        except:
            print "Error in", sampleName
        
    conn.close()
    #if not precalc:
    #    dump.dump(usedOTUs, "%s/usedOTUs.pcl" % datadir)
    #    dump.dump(envo2samples, "%s/envo2samples.pcl" % datadir)
    #    dump.dump(ecosystems2samples, "%s/ecosystems2samples.pcl" % datadir)
