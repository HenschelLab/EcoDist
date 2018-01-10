
"""
Ecosystem distribution for a set of 16S rRNA profiles with OTUs called against GreenGenes 13.5.
Usage: python ecoDistSQLpost.py <biom-file> <output directory>

Creates an ecodistribution plot for a sample similar to qiime barplots plus a second dimension:
for each OTU, the distribution over ecosystems is color coded/visualized
Moreover ecosystem distribution entropy for each OTU in the sample is calculated

for each otu: create a stats of Ecosystem occurrences
SQL based version: OTU distribution is calculated with respect to global database (as opposed to Biom table distribution)
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
import time
from DBsettings import db

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
    #figlegend.legend(legendbarsPhyla, [p.split(";")[-1] for p in phylaNames], 'center')
    fig.savefig("%s.svg"%filename)
    fig.savefig("%s.png"%filename)
    figlegend.savefig("%s_leg.svg"%filename)
    figlegend.savefig("%s_leg.png"%filename)
    fig.clf()
    figlegend.clf()
    plt.close(fig)
    plt.close(figlegend)
                               
class OTUfromFile:
    def __init__(self, id, count):
        self.id = otuTable.ids(axis='observation')[idx]
        self.count = count
        filename = "%s/otu_%s.pcl" % ( otudatadir, self.id)
        if os.path.exists(filename):
            self.lineage, cnt, self.phylum, self.ecoDistribution, self.entropy = dump.load(filename)
            self.hasData = self.ecoDistribution.sum() > 0
        else:
            #print "Warning: otu %s not found" % self.id
            self.hasData = False
    def __cmp__(self, o):
        return cmp(self.lineage, o.lineage)

if __name__ == "__main__":
    ## Settings
    if len(sys.argv)!=3:
        print "Wrong number of arguments!\n"
        print __doc__
        
    ecosystemColors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal/Human': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Air': 'AliceBlue', 'Hypersaline':'r'}
    ecosystems = ecosystemColors.keys()
    ecosystemsIndex = dict([(eco,idx) for idx,eco in enumerate(ecosystems)])
    datadir = "data" #"/data/EarthMicrobiomeProject/OTUdata/"
    otudatadir = "/bioinfo-omics/Data/Projects/EcoPhyl/EcoDist/Data/OTUdata" ## should be moved!
    resultdir = sys.argv[2] 
    
    phylaColorDict = dump.load("%s/phylaColorDict.pcl" % datadir)
    phylaColorDict['?'] = '#f0f0f0'
    
    cols4plot = [ecosystemColors[ecosystem] for ecosystem in ecosystems] + ['b']

    ## Biom table
    biom = sys.argv[1]

    otuTable = parse_biom_table(open(biom, 'U'))
    print "Parsed otu table %s" % biom

    mincount = 0
    counter = 0
    for sampleData, sampleName, m in otuTable.iter(axis="sample"):
        counter += 1
        print "Creating %s OTU objects for sample %s " %(len( np.where(sampleData > mincount)[0]), sampleName)
        selectedOTUs = [OTUfromFile(idx, sampleData[idx]) for idx in np.where(sampleData > mincount)[0]]
        print "Selected %s OTUs" % len(selectedOTUs)
        selectedOTUs = [o for o in selectedOTUs if o.hasData]
        print "Selected %s OTUs" % len(selectedOTUs)
        selectedOTUs.sort()
        otudata = [otu.ecoDistribution for otu in selectedOTUs]
        widths = np.array([otu.count for otu in selectedOTUs]) #np.ones(len(otudata)) * 1 ## possibly stretch by occ.
        phylaNames, phylaCount = zip(*[(k, sum([e.count for e in g])) for k, g in groupby(selectedOTUs, lambda o:o.phylum)])
        entropies =  [otu.entropy for otu in selectedOTUs]
        ecodistributionPlot(otudata, widths, cols4plot, ecosystems + ['Misc'], np.array(phylaCount), phylaNames, entropies, "%s/%s" % (resultdir, sampleName))
