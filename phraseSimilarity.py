import nltk
from collections import defaultdict
import dump
import numpy as np
import pdb
import time

class PhraseSimilarity:
    def __init__(self, phrases, stem=False):
        self.stem = stem
        self.phrases = phrases        
        self.tokenSets = [set(nltk.word_tokenize(s.lower())) for s in phrases] ## preserve order between phrases and tokenSets, for reusing indexes
        self.stemmer = nltk.PorterStemmer()
        if stem:
            self.tokenSets = [set([self.stemmer.stem(token) for token in tokenSet]) for tokenSet in self.tokenSets]
        self.indexing()
        self.wordFreqs = dump.load("/home/zain/Projects/NLP/Data/wordFreqs_brown.pcl")
        #self.wordFreqs = dump.load("data/wordFreqs_brown.pcl")
        
        self.D = float(sum(self.wordFreqs.values()))
        self.maxIdf = self.idf(0.5)
        tmpW = defaultdict(int)
        if stem:
            for (w,f) in self.wordFreqs.items():
                tmpW[self.stemmer.stem(w)] += f
            self.weights = dict([(w, self.idf(f)) for (w,f) in tmpW.items()])
            #self.weights = dict([(self.stemmer.stem(w), self.idf(f)) for (w,f) in self.wordFreqs.items()])
        else:
            self.weights = dict([(w, self.idf(f)) for (w,f) in self.wordFreqs.items()])
        self.phraseWeights = [self.weightByFreq(phrase) for phrase in self.tokenSets]

    def indexing(self):
        self.tokenDict = defaultdict(list)
        for idx, phrase in enumerate(self.tokenSets):
            for token in phrase:
                self.tokenDict[token].append(idx)

    def weightByFreq(self, tokens):
        return sum([self.weights.get(t, self.maxIdf) for t in tokens])
    
    def idf(self, f):
        return np.log(1 + (self.D/f))

    def fullMatch(self, phrase, tokenize=False, stem=False, idxOnly=False):
        if tokenize:
            phrase = nltk.word_tokenize(phrase)
        if stem:
            phrase = [self.stemmer.stem(token) for token in phrase]
        overlap = None
        for token in phrase:
            phrasesWithToken = set(self.tokenDict[token])
            overlap =  phrasesWithToken if overlap is None else overlap.intersection(phrasesWithToken)
        if idxOnly: return overlap
        else:
            return [self.phrases[idx] for idx in overlap]
    
    def fullMatchJaccard(self, phrase, tokenize=False):
        fullMatches = self.fullMatch(phrase, tokenize=tokenize, stem=self.stem, idxOnly=True)
        if fullMatches:
            return self.jaccard(phrase, tokenize=tokenize, phraseSelection=fullMatches)
        else:
            return "(no match)", 0
            
    def jaccard(self, phrase, tokenize=False, phraseSelection="all"):        
        def ratio(w1, w2, overlap):
            return overlap/float(w1 + w2 - overlap)
        if tokenize:
            phrase = nltk.word_tokenize(phrase)
        elif not type(phrase) == type([]):
            phrase = list(phrase)
        if self.stem:
            phrase = [self.stemmer.stem(token) for token in phrase]

        phraseWeight =  self.weightByFreq(phrase)
        phraseOverlaps = defaultdict(int)
        for token in phrase:
            for phraseIdx in self.tokenDict[token]:            
                phraseOverlaps[phraseIdx] += self.weights.get(token, self.maxIdf)
        
        
        if phraseSelection == "all":
            similarities = [(ratio(phraseWeight, self.phraseWeights[phraseIdx], overlap), phraseIdx) for phraseIdx, overlap in phraseOverlaps.items()]
        else:
            similarities = [(ratio(phraseWeight, self.phraseWeights[phraseIdx], overlap), phraseIdx) for phraseIdx, overlap in phraseOverlaps.items() if phraseIdx in phraseSelection]
        if similarities:
            score, bestMatchIdx = max(similarities)
            return self.phrases[bestMatchIdx], score
        else:
            return "(no match)", 0
    
if __name__ == "__main__":    
    t0 = time.time()
    dataDir = "/home/zain/Projects/OttoTextMining/Data"
    ontoDir = "/home/zain/Projects/OttoTextMining/OntologyData"
    synonymDict = dump.load("%s/envoTerms2.pcl" % ontoDir)
    synonyms = synonymDict.keys()
    t1 = time.time()
    ps = PhraseSimilarity(synonyms)
    t2 = time.time()
    habitats = [line.rstrip() for line in open('%s/sample_details.environments.txt' % dataDir)]
    #for hab in habitats:
    #    hab, ps.jaccard(hab, tokenize=True) 
    for hab in habitats:
        matches = ps.fullMatchJaccard(hab, True, True)
        if matches:
            print hab, matches
    t3 = time.time()
    print  "%.3f sec" % (t1-t0)
    print  "%.3f sec" % (t2-t1)
    print  "%.3f sec" % (t3-t2)
