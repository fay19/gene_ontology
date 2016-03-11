## A class to handle reading a list of stopwords and providing quick lookup
class StopwordList:
    ## Constructor
    # @param fname The filename of the stopwords to read.  File should
    # contain one word per line and already be lowercase
    def __init__(self, fname=None):
        self.words = set([])
        if fname:
            self.readFile(fname)
        

    ## Read a file that contains stopwords
    # @param fname The filename of the stopwords to read.  File should
    # contain one word per line and already be lowercase            
    def readFile(self, fname):
        words = []
        try:
            f = open(fname, 'r')
            line = f.readline()
            while line:
                word = line.strip()
                if word != "":
                    words.append(word)
                line = f.readline()
            f.close()
        except Exception, e:
            msg = "Could not read file \"%s\": %s\nWill NOT be using any stopwords"
            print msg % (fname, str(e))
        self.addWords(words)


    ## Handle queries like "word in stopwords"
    # \return True only if word is in this list
    def __contains__(self, word):
        return word in self.words


    ## Same as __contains__
    # \return True only if word is in this list
    def hasWord(self, word):
        return word in self


    ## Add a list of words to the stopwords
    # @param words A set or list of words
    def addWords(self, words):
        self.words.update(words)


    ## Get the Python list that this class operates on
    # \return a Python list of the words contained in this list
    def getList(self):
        return [word for word in self.words]


    ## Return a new StopwordList composed of the words in this one and
    # another.
    # \return A new StopwordList
    def __add__(self, other):
        sl = StopwordList()
        sl.words.update(self.words, other.words)
        return sl


    ## Determine if another stopwordlist is the same as this one
    # \return True if same, False otherwise
    def __eq__(self, other):
        return self.words == other.words


    ## Determine if another stopwordlist is the same as this one
    # \return False if same, True otherwise
    def __neq__(self, other):
        return not self == other


    ## Add a list of stopwordlists together to form one superlist
    # @param stopwordLists A list of StopwordList objects
    # \return A StopwordList that is the set that contains the union of the ones given
    @classmethod
    def join(klass, stopwordLists):
        sl = StopwordList()
        words = [swl.words for swl in stopwordLists]
        sl.words.update(*words)
        return sl

## Removes from the graph the genes that are not included in a given list of genes
# @param graph GOGeneGraph or GOGenePubmedGraph from which to remove the genes
# @param geneList The genes that should be kept if they are listed in the graph
# /return graph The graph with the unwanted genes removed
def keepGenes(graph, geneList):
    for node in graph:
        remove = list()
        for gene in graph.node[node]['data'].getPropagatedGenes():
            if gene[0] not in geneList:
                remove.append(gene)
        for gene in remove:
            graph.node[node]['data'].propGenes.remove(gene)
            
        remove = list()
        for gene in graph.node[node]['data'].getGenes():
            if gene[0] not in geneList:
                remove.append(gene)
        for gene in remove:
            graph.node[node]['data'].genes.remove(gene)

        allGenes = set(graph.geneToNode.keys())
        remove = list(allGenes.difference(set(geneList)))
        for gene in remove:
            del graph.geneToNode[gene]
        
    return graph
