from xml.sax import handler, make_parser
import re
from Document import Document


## This class parses PubmedArticleSet xml files.  When parsing is finished, the docs
# instance variable contains a dictionary of the form {'pmid': Document}
class PubmedArticleSet(handler.ContentHandler):
    def __init__(self):
        handler.feature_external_ges = "false"
        self.docs = {}
        self.doc = None
        self.chars = ""
        
               
    def startElement(self, name, attr):
        if name == 'PubmedArticle' or name == 'PubmedBookArticle':
            self.doc = Document()            
        self.chars = ""

            
    def endElement(self, name):
        if name == 'PubmedArticle':
            self.docs[self.doc.pmid] = self.doc
        if name == 'PMID' and self.doc.pmid == None:
            self.doc.pmid = self.text()
        if name == 'ArticleTitle':
            self.doc.title = self.text()
        if name == 'AbstractText':
            if self.doc.abstract == None:
                self.doc.abstract = self.text()
            else:
                self.doc.abstract += self.text()
        if name == 'DescriptorName':
            self.doc.addMeSH(self.text())
            
    
    def characters(self, data):
        self.chars += data


    def text(self):
        return self.chars.strip().encode('ascii', 'ignore')
    

    ## Method to parse a PubmedArticleSet XML file.
    # @param location The location of the xml file to parse
    # return A PubmedArticleSet object 
    @classmethod
    def parse(self, location):
        parser = make_parser()
        parser.setFeature("http://xml.org/sax/features/external-general-entities", False)
        parser.setFeature("http://xml.org/sax/features/external-parameter-entities", False)
        handler = PubmedArticleSet()
        parser.setContentHandler(handler)
        try:
            f = open(location, 'r')
            parser.parse(f)
            f.close()
        except Exception, e:
            raise RuntimeError, "Could not parse PubmedArticleSet XML file at %s" % location
        return handler

            
                                                                        



