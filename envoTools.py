"""
Ontology class that provides useful functionality for:
   - finding sub-categories (children)
   - parent classes
   - keeps track of annotations (both ways: habitats (Envo-IDs) <-> samples)

Closely tied to Habitat class in habitatOntologyMapping.py (used in methods findBestMatch and

Author: Andreas Henschel


"""

import cPickle
import networkx as nx
from itertools import chain

class Ontology:
    def __init__(self, synonymFile, graphFile, color):
        self.color = color
        f = open(synonymFile)
        self.synonymDict = cPickle.load(f)
        f.close()        
        self.G = nx.read_gpickle(graphFile)
        self.synonyms = self.synonymDict.keys()
    def initPhraseSimilarity(self):
        from phraseSimilarity import PhraseSimilarity
        self.phraseSim = PhraseSimilarity(self.synonyms, stem=True)
        
    def ontoInfo(self, oboIDs):
        def info(oboID):
            if not self.G.node.has_key(oboID):
                return ""
            n = self.G.node[oboID]
            ancestry = "-".join([info(rel[1]) for rel in n['relation'] if rel[0] in ["is_a", "part_of"]]) if n.has_key("relation") else "."
            name = n["name"] if n.has_key("name") else ""
            return "%s-%s" % (name, ancestry) 
        return "/".join([info(oboID) for oboID in oboIDs]) ## should be only 1 though!
        
    def findBestMatch(self, habitat, method="jaccard"):
        ## requires initPhraseSimilarity to be run!
        if method=="jaccard":
            bestmatch, score = self.phraseSim.jaccard(habitat.bagOfWords)
        elif method == "fullmatch":
            bestmatch, score = self.phraseSim.fullMatchJaccard(habitat.bagOfWords)
        if score > 0:
            info = self.ontoInfo(self.synonymDict[bestmatch])
            habitat.annotate(bestmatch, [self.color, bestmatch, score, info, self.synonymDict[bestmatch][0]])
            self.annotate(bestmatch, habitat)
        else:
            habitat.annotate("", [self.color, "", 0, "", ""])
    def annotate(self, bestmatch, habitat):
        node = self.G.node[self.synonymDict[bestmatch][0]]
        node["habitats"] = node.get("habitats", []) + [habitat]
    def getSubsumedHabitats(self, node):
        h = list(chain.from_iterable([self.getSubsumedHabitats(cnode) for cnode in self.G.predecessors(node)])) 
        return h + self.G.node[node].get("habitats", [])
    def getChildren(self, node):
        h = list(chain.from_iterable([self.getChildren(cnode) for cnode in self.G.predecessors(node)])) 
        return h + [node]
    #def getParents(self, node):
    #    h = list(chain.from_iterable([self.getChildren(cnode) for cnode in self.G.predecessors(node)])) 
    #    return h + [node]
    def getRelevantChildren(self, node):
        h = list(chain.from_iterable([self.getRelevantChildren(cnode) for cnode in self.G.predecessors(node)]))
        if self.G.node[node].has_key("habitats"):
            return h + [node]
        else:
            return h
    def listHabitats(self, nodeName):
        ## assumes that a surrounding table is produced
        habs = []
        subsumedHabitats = self.getSubsumedHabitats(self.synonymDict[nodeName][0])
        for hab in subsumedHabitats:
            if not hab.sampleEventID in habs:
                habs.append(hab.sampleEventID)
                print hab.html

if __name__ == "__main__":
    #just testing, adjust file paths!
    ontoDir = "OntologyData" ## Might need configuration!
    envo = Ontology("%s/envoTerms3.pcl" % ontoDir, "%s/envo3.pcl" % ontoDir, "envo")
    print nx.is_directed_acyclic_graph(envo.G)
