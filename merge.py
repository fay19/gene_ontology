import MySQLdb
import sys
import prettytable
from scipy.stats import hypergeom
import collections
import cPickle as pickle
import csv
import itertools
import os, csv, sys, math, time
import xml.dom.minidom
import networkx

def checkDup(P_list):
	PTp = P_list
	PTp.sort()
	for i in (range(len(PTp)-1)):
		if PTp[i]==PTp[i+1]:
			return 1
	return 0

def Union(P_list1, P_list2):
	L_rets = []
	L_dic = {}
	for item in P_list1:
		L_rets.append(item)
		L_dic[item] = 1
	for item in P_list2:
		if not(L_dic.has_key(item)):
			L_rets.append(item)
	L_rets.sort()
	return L_rets

def Intersection(P_list1, P_list2):
	L_rets = []
	L_dic = {}
	for item in P_list1:
		L_dic[item] = 1
	for item in P_list2:
		if L_dic.has_key(item):
			L_rets.append(item)
	L_rets.sort()
	return L_rets

def Difference(P_list1, P_list2):
	#return P_list2 - P_list1
	L_rets = []
	L_dic = {}
	for item in P_list1:
		L_dic[item] = 1
	for item in P_list2:
		if not(L_dic.has_key(item)):
			L_rets.append(item)
	L_rets.sort()
	return L_rets

# # hypergeometric function for calculating probability
#     # a -- associated gene in input genesubset
#     # b -- population (here use the offical home sapiens genes in NCBI)
#     # b -- assoicated genes in population
#     # d -- input genelist size
# def hypergeom(self,a,b,c,d):
# 	prob = 1-hypergeom.cdf(a-1, b, c, d)
# 	return prob

def hyperGeo_pr(p_x,p_k,p_m,p_n):
	#Function to calculate the probability
	#input:
	#	p_x: number of white balls drawn.
	#	p_k: number of balls drawn (note p_x is not larger than p_k).
	#	p_m: number of white balls in the urn.
	#	p_n: number of black balls in the urn.
	#output:
	#	the probability that has p_x white balls in a drawing of p_k balls.

	l_nue = []
	l_deno = []
	for j in range(p_x):
		l_nue.append(p_m-j)
		l_deno.append(j+1)
	for j in range(p_k-p_x):
		l_nue.append(p_n-j)
		l_deno.append(j+1)
	for j in range(p_k):
		l_deno.append(p_m+p_n-j)
		l_nue.append(j+1)
	l_deno.sort()
	l_nue.sort()
	#l_p = 1.0
	l_p2 = 0.0

	try:
		for j in range(len(l_deno)):
			#l_p = l_p*l_nue[j]/l_deno[j]
			l_p2 += math.log(l_nue[j]) - math.log(l_deno[j])
	
	except:
		return 0
	
	#return math.exp(l_p)
	return math.exp(l_p2)
	
def hyperGeo_Pvalue(white_ball_drawn,total_ball_drawn,white_ball_in_urn,total_ball_in_urn):
	#Function to calculate the p_value
	#input:
	#	white_ball_drawn: number of white balls drawn.
	#	white_ball_in_urn: number of white balls in the urn.
	#	total_ball_in_urn: number of black balls in the urn.
	#	total_ball_drawn: number of balls drawn (note white_ball_drawn is not larger than total_ball_drawn).
	#output:
	#	1-the probability that has at least white_ball_drawn white balls in a drawing of total_ball_drawn balls.
	l_p = 0.0
	for j in range(white_ball_drawn):
		hptp = hyperGeo_pr(j,total_ball_drawn,white_ball_in_urn,total_ball_in_urn-white_ball_in_urn)
		l_p += hptp
	if l_p>1:
		l_p = 1.0
	return (1-l_p)	
	
class GoGraph(networkx.DiGraph):
	def __init__(self,P_weightGraph, host, user, password, dbname):
		networkx.DiGraph.__init__(self)		
		self.association = 'Association_proteins'
		self.association_acc = 'Association_Accumulation'
		self.mapping = 'mapping_protiens'
		self.PV = 'P_value'
		self.TotalLost = 0.0
		self.db = MySQLdb.connect(host, user, password, dbname)
		self.gene_GOterm = self.geneToGOtermAssoc()
		self.GOterm_gene = self.GOtermToGeneAssoc()
		self.NumberOfAllProtein = len(self.gene_GOterm)
		self.createDiGraph(P_weightGraph)

	def geneToGOtermAssoc(self):
		gene_GOterm = {}
		cursor = self.db.cursor()
		query = "select distinct(offsymbol) from final_symbol_term"
		cursor.execute(query)
		query_result = cursor.fetchall()
		for gene in query_result:
			query2 = "select acc from final_symbol_term where offsymbol='%s'"%(gene[0])
			cursor.execute(query2);
			query2_result = cursor.fetchall()
			assocTerms = []
			for term in query2_result:
				assocTerms.append(term[0])
			gene_GOterm[gene[0]] = assocTerms
		# print len(gene_GOterm)
		return gene_GOterm

	def GOtermToGeneAssoc(self):
		GOterm_gene = {}
		cursor = self.db.cursor()
		query = "select distinct(acc) from final_symbol_term"
		cursor.execute(query)
		query_result = cursor.fetchall()
		for term in query_result:
			query2 = "select offsymbol from final_symbol_term where acc='%s'"%(term[0])
			cursor.execute(query2);
			query2_result = cursor.fetchall()
			assocGenes = []
			for gene in query2_result:
				assocGenes.append(gene[0])
			GOterm_gene[term[0]] = assocGenes
		# print GOterm_gene
		# print len(GOterm_gene)
		return GOterm_gene

	def readWeights(self, P_filename = ''):
	#read edge-weighted go structure data
		xmldoc = xml.dom.minidom.parse(P_filename)
		xmlData = xmldoc.firstChild
		goTermList = xmlData.getElementsByTagName('term')

		goGraphRe = {}
		for term in goTermList:
				go = term.attributes['id']
				goGraphRe[go.value] = {}
				parentList = term.getElementsByTagName('parent')
				for par in parentList:
					parGO = par.attributes['id']
					goGraphRe[go.value][parGO.value] = {}
					edgeList = par.getElementsByTagName('edge')
					for edge in edgeList:
						edgeType = edge.attributes['type']
						goGraphRe[go.value][parGO.value] = float(edge.firstChild.data)
						# print go.value, parGO.value, float(edge.firstChild.data)
		return goGraphRe

	def createDiGraph(self, P_filename):
		'''
		create edge-weighted Goterm structure
		'''
		GoStructure = self.readWeights(P_filename)
		for child in GoStructure:
			for parent in GoStructure[child]:
				self.add_edge(child, parent)
				self.edge[child][parent]['weight'] = GoStructure[child][parent]
		self.node_AssociationProteinData()
		self.node_Depth()

	def propagate_AssociatedProteins(self):
		Top_sort_nodes = networkx.topological_sort(self)
		for nodeTp in Top_sort_nodes:
			ParentTp = self.successors(nodeTp)
			AssoSelf = self.node[nodeTp][self.association_acc]
			for nodePare in ParentTp:
				AssoPare = self.node[nodePare][self.association_acc]
				#print nodeTp,nodePare,AssoSelf,AssoPare
				UnionTwo = Union(AssoSelf,AssoPare)
				self.node[nodePare][self.association_acc] = UnionTwo

	def node_AssociationProteinData(self):
		for gotermTp in self.nodes():
			self.node[gotermTp]['level'] = 1000000000
			if self.GOterm_gene.has_key(gotermTp):
				self.node[gotermTp][self.association] = self.GOterm_gene[gotermTp]
				self.node[gotermTp][self.association_acc] = self.GOterm_gene[gotermTp]
			else:
				self.node[gotermTp][self.association] = []
				self.node[gotermTp][self.association_acc] = []
			# print gotermTp
			# print self.node[gotermTp][self.association]
			# print self.node[gotermTp][self.association_acc]
			self.node[gotermTp][self.mapping] = []
		self.propagate_AssociatedProteins()

	def node_Depth(self):
		Top_sort_nodes = networkx.topological_sort(self, reverse = True)
		self.node[Top_sort_nodes[0]]['level'] = 1
		for nodeTp in Top_sort_nodes:
			ChildrenTp = self.predecessors(nodeTp)
			for chd in ChildrenTp:
				if self.node[chd]['level']>self.node[nodeTp]['level']+1:
					self.node[chd]['level'] = self.node[nodeTp]['level']+1
		
	
	def node_Mapping_ProteinList_PV2(self, P_protein_InList):
	#map gene list to Go Terms, use association_acc to calculate PV
		Gene_List_Associated = []
		for gene in P_protein_InList:
			if self.gene_GOterm.has_key(gene):
				Gene_List_Associated.append(gene)
		self.NumberOfProtein_InList = len(Gene_List_Associated)

		for nodeID in self.nodes():
			P_Association = self.node[nodeID][self.association]
			P_Association2 = self.node[nodeID][self.association_acc]
			P_protein = Intersection(P_protein_InList,P_Association)
			self.node[nodeID][self.mapping] = P_protein
				
			newPV = hyperGeo_Pvalue(len(P_protein),self.NumberOfProtein_InList, len(P_Association2),self.NumberOfAllProtein)
			self.node[nodeID][self.PV] = newPV
			
			#---
			if checkDup(P_Association)==1:
				print 'Asso dup'
			if checkDup(P_protein_InList)==1:
				print 'PList dup'
			if checkDup(P_protein)==1:  
				print 'Intersection dup'
		
	def symplify_Goterm_Structure(self):
		Top_Order = networkx.topological_sort(self)
		for nodeTp in Top_Order:
			if len(self.node[nodeTp][self.mapping])==0 and len(self.predecessors(nodeTp)) == 0:
				self.remove_node(nodeTp)

	def calculate_lost(self, P_nodeSource, P_nodeTarget):
		edgeWeight = self.edge[P_nodeSource][P_nodeTarget]['weight']
		L_lost = edgeWeight*len(self.node[P_nodeSource][self.mapping])
		#print 'edgeWeight=',edgeWeight, 'sizeOfSource=',len(self.node[P_nodeSource][self.mapping]),'lost=',L_lost
		return L_lost

	def merge_nodes2(self, P_nodeSource, P_nodeTarget):
	#merge association protien in source node into target node,
	#update the P_value of the target node, and delete source node from the graph
	#use self.association_acc to calculate PV

		L_lost = self. calculate_lost(P_nodeSource, P_nodeTarget)
		P_Source = self.node[P_nodeSource][self.association]
		P_Target = self.node[P_nodeTarget][self.association]
		self.node[P_nodeTarget][self.association] = Union(P_Source, P_Target)

		#---
		if checkDup(self.node[P_nodeTarget][self.association]) == 1:
			print '---Asso dup'

		P_Source2 = self.node[P_nodeSource][self.mapping]
		P_Target2 = self.node[P_nodeTarget][self.mapping]
		self.node[P_nodeTarget][self.mapping] = Union(P_Source2, P_Target2)

		#---
		if len(self.node[P_nodeTarget][self.mapping])>self.NumberOfProtein_InList:
			print len(P_Source2), len(P_Target2)
		if checkDup(self.node[P_nodeTarget][self.mapping]) ==1:
			print '!!! mapping Dup'

		#print 'All protiens=',self.NumberOfAllProtein, 'All white protien#=',len(self.node[P_nodeTarget][self.association])
		#print 'whiteDraw=',len(self.node[P_nodeTarget][self.mapping]), 'Alldraw = ',self.NumberOfProtein_InList
		newPV = hyperGeo_Pvalue(len(self.node[P_nodeTarget][self.mapping]),self.NumberOfProtein_InList, len(self.node[P_nodeTarget][self.association_acc]),self.NumberOfAllProtein)
		self.node[P_nodeTarget][self.PV] = newPV
		self.TotalLost += L_lost
		Children_S = self.predecessors(P_nodeSource)
		Children_T = self.predecessors(P_nodeTarget)
		C_S_sub_T = Difference(Children_T,Children_S)
		#print 'Children of S not of T=',C_S_sub_T
		for proteinTp in C_S_sub_T:
			self.add_edge(proteinTp,P_nodeTarget)
			self.edge[proteinTp][P_nodeTarget]['weight'] = self.edge[proteinTp][P_nodeSource]['weight'] + self.edge[P_nodeSource][P_nodeTarget]['weight']
		self.remove_node(P_nodeSource)

	def printResult(self):
		'''
		print 'The total information lost is:', self.TotalLost
		print 'The final goterm and associated protein list:'
		for gotp in self.nodes():
			if len(self.node[gotp][self.mapping])>0:
				print gotp,':',self.node[gotp][self.mapping]
		'''
		L_Goterms = {}
		for gotp in self.nodes():
			if len(self.node[gotp][self.mapping])>0:
				L_Goterms[gotp] = [self.node[gotp][self.mapping],self.node[gotp]['level']]
		return [self.TotalLost,L_Goterms]

	def gotermSummarization(self, P_geneList, P_value, P_minProteinNo):
	#use P_GotermNumber to summarizate a given list of proteins
	#Topologically sort at first

		self.node_Mapping_ProteinList_PV2(P_geneList)
		self.symplify_Goterm_Structure()
		
		#Topologically sort all node, the first element is a leaf and the last element is the root
		New_Top_sort_nodes = networkx.topological_sort(self)

		#go through all nodes in topological order from leaves to the root.
		for nodeTp_child in New_Top_sort_nodes:
			nodePV = self.node[nodeTp_child][self.PV]
			nodeSize = len(self.node[nodeTp_child][self.mapping])

			#if the current node's pv or size is not constrained by the setting conditions, merge this node to the most close parent node.
			if nodePV>P_value or nodeSize<P_minProteinNo:
				L_parrents = self.successors(nodeTp_child)
				L_MinWeight = 100000; L_MinPa = 'none'
				for nodetp in L_parrents:
					if self.edge[nodeTp_child][nodetp]['weight']<L_MinWeight:
						L_MinWeight = self.edge[nodeTp_child][nodetp]['weight']
						L_MinPa = nodetp
				if L_MinPa == 'none':
					pass
				else:
					self.merge_nodes2(nodeTp_child,L_MinPa)

		return self.printResult()


def main():
	gograph = sys.argv[1]
	host = sys.argv[2]
	user = sys.argv[3]
	password = sys.argv[4]
	dbname = sys.argv[5]
	GeneList = ['CLEC5A','COL11A1','CRABP2','CXCL10','DCC1','DDAH2','DEFB1','DEPDC1','DKFZp762E1312','DLG7','DNA2L','DSC2','DSC3','ECT2','ENC1','EXOSC5','EYA4','EZH2','FABP5','FAM64A','FGF9','FLJ21963','GFOD1','GGH','GINS1','GLDC','GOLGA8A','GOLT1B','HMMR','HRASLS','HS3ST3A1','IGFBP2','ING4','INHBB','ISG15','ITGB8','KIAA1199','KIF11','KIF14','KIF20A','KIF23','KIF2C','KIFC1','KLK7','KPNA2','KRT6A','LMNB1','MAL','MELK','MFAP2','MGAT4A','MKI67','MRS2L','MUC16','MYO10','NCAPD2','NDC80','NEK2','NETO2','NID2','NLRP2','NMU','NRAS','NRCAM','NUSAP1','OGDHL']
	gograph = GoGraph(gograph, host, user, password, dbname)
	Result = gograph.gotermSummarization(GeneList,0.05,3)
	print 'Total information lost is', Result[0]
	for goterm in Result[1]:
		print 'GO term ID:', goterm, '--------------'
		print 'GO term level:',Result[1][goterm][1]
		print 'Gene list:',Result[1][goterm][0],'\n'

if __name__ == "__main__":
	main()	
	


'''
stTime = time.mktime(time.gmtime())

#------------------------------------------------------------------------------------------------------
#Example of using the class GoGraph:

#Data of Go Ontology structure and gene_Goterm association
weightGographData = 'newWeightedPubMedGO.xml'
geneGotermAssociationData = 'gene_association.goa_human_2012'


#Create a GoGraph object (Node: every time you use the gotermSummarization(), you need to create a new object)
G = GoGraph(weightGographData,geneGotermAssociationData)

#A list of genes need to be summarized
GeneList = ['CLEC5A','COL11A1','CRABP2','CXCL10','DCC1','DDAH2','DEFB1','DEPDC1','DKFZp762E1312','DLG7','DNA2L','DSC2','DSC3','ECT2','ENC1','EXOSC5','EYA4','EZH2','FABP5','FAM64A','FGF9','FLJ21963','GFOD1','GGH','GINS1','GLDC','GOLGA8A','GOLT1B','HMMR','HRASLS','HS3ST3A1','IGFBP2','ING4','INHBB','ISG15','ITGB8','KIAA1199','KIF11','KIF14','KIF20A','KIF23','KIF2C','KIFC1','KLK7','KPNA2','KRT6A','LMNB1','MAL','MELK','MFAP2','MGAT4A','MKI67','MRS2L','MUC16','MYO10','NCAPD2','NDC80','NEK2','NETO2','NID2','NLRP2','NMU','NRAS','NRCAM','NUSAP1','OGDHL']


#Using Go term to summarize the list of gene.
Result =  G.gotermSummarization(GeneList,0.05,3)

#0.05 is the threshold of P_Value of the Go term node in final result.
#5 is the minimum number of genes in the Go term node in final result.
#The result has the format: [value_0,{Goterm_1:[a list of genes_1,value_1],Goterm_2:[a list of genes_2,value_2],....}]
# value_0: the total information lost for the summarization.
# Goterm_1: Goterm ID that is used to summarize the given list of gene.
# a list of genes_1: a subset of genes (in the given list of gene) that are annotated by Goterm_1.
# value_1: the level of Goterm_1 on the Go Ontology. The root note is in level 1.



#print the result
print 'Total information lost is', Result[0]
for goterm in Result[1]:
	print 'GO term ID:', goterm, '--------------'
	print 'GO term level:',Result[1][goterm][1]
	print 'Gene list:',Result[1][goterm][0],'\n'

#------------------------------------------------------------------------------

print '------------------------------------------------------------'
endTime = time.mktime(time.gmtime())
print 'Total time =', ((endTime-stTime)/60), 'minutes'
'''