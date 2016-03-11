from xml.etree import ElementTree
import xml.etree.cElementTree as ET

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

for i in range(41):
	with open("../pubmedAbstracts/block_%s.xml"%i, 'rt') as f:
		oldtree = ElementTree.parse(f)
	root=ET.Element("PubmedArticleSet");

	for article in oldtree.iter("MedlineCitation"):
		doc=ET.SubElement(root,"PubMedArticle");
		pmid=article.find("PMID").text
		print pmid
		title=article.find("Article").find("ArticleTitle").text
		ET.SubElement(doc, "PMID").text = pmid
		ET.SubElement(doc, "Title").text = title
		if article.find("Article").find("Abstract") is not None:
			abstract=article.find("Article").find("Abstract").find("AbstractText").text	
			ET.SubElement(doc, "Abstract").text = abstract
		elif article.find("Article").find("Abstract") is None:
			pass
	newtree = ET.ElementTree(root)
	indent(root)
	newtree.write("../parsedAbstracts/files_%s.xml"%i)


