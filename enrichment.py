import MySQLdb
import sys
import prettytable
from scipy.stats import hypergeom
import collections
import cPickle as pickle
import csv
import itertools

"""
This class is calcaulating enrichment of given gene list

Example of how to use this class:

inputfile_path = "MLL2-MLL3.targetgenes.v9.csv". Inputfile is genelists in a csv file: every row is a list, the first column is drivers of this gene list.  
outputfile_path = "". The number of outputfiles is same with the numbers of genelists in input file. Each output file is named by the driver of each genelist. 
Output file will be contained in the givend outputfile_path directory.
threshold = 0.01

Create an enrichment object: enrichment = enrichment(host, user, password, dbname, inputfile_path, outputfile_path, threshold)
Call fuction enrichment.enrich_csv()
		
"""
class enrichment:
	
	## Class constructor. It has two fields: a database connector and numbers of gene population
	# @param host
	# @param user
	# @param password
	# @param dbname
	# @param inputfile_path
	# @param outputfile_path
	# @param threshold
	def __init__(self, host, user, password, dbname, inputfile_path, outputfile_path, threshold):
		"""
		Class constructor. It has 4 fields: a database connector, numbers of gene population, the filepath of clients' csv file, and threshold on enrichment
		@param host
		@param user
		@param password
		@param dbname
		@param inputfile_path
		@param threshold
		"""
		self.db=MySQLdb.connect(host, user, password, dbname)
		self.M = self.numbersOfPopulation()
		self.inputfile= inputfile_path
		self.outputfile = outputfile_path
		self.threshold = threshold

	## Return the number of gene population, which is the number of homo sapiens genes in NCBI gene symbol set
	def numbersOfPopulation(self):
		"""
		return the number of gene population, which is the number of homo sapiens genes in NCBI gene symbol set
		"""
		cursor = self.db.cursor()
		query = "SELECT count(*) FROM NCBIsymbol AS N where N.symbol IN (SELECT G.symbol FROM gene_product AS G, species AS S where G.species_id=S.id AND S.genus = 'Homo' AND S.species = 'sapiens') "
		cursor.execute(query)
		query_result = cursor.fetchall();
		return query_result[0]

	## Return numbers of genes that related to specific term in input gene set
	# @param term
	# @geneList
	def assocGeneSubsetNum(self, term, geneList):
		"""
		Return numbers of genes that related to specific term in input gene set
		@param term
		@geneList
		"""
		cursor=self.db.cursor()
		query="select count(distinct(offsymbol)) from final_symbol_term where acc='%s' and offsymbol in(%s)"%(term, geneList)
		cursor.execute(query)
		query_result=cursor.fetchall();
		return query_result[0][0]
		
	## Return numbers of genes that related to specific term in whole gene set
	# @param term
	def assocGeneNum(self,term):
		"""
		Return numbers of genes that related to specific term in whole gene set
		@param term
		"""
		cursor=self.db.cursor()
		query="select n from term_n where acc ='%s' "%(term)
		cursor.execute(query)
		query_result=cursor.fetchall()
		return query_result[0][0]

	## Return a list of genes associated to specific term in input gene set
	# @param term
	# @param geneList
	def assocGeneSubset(self, term, geneList):
		"""
		Return a list of genes associated to specific term in input gene set
		@param term
		@param geneList	
	    """
		cursor=self.db.cursor()
		query="select distinct(offsymbol) from final_symbol_term where acc='%s' and offsymbol in(%s)"%(term,geneList)
		cursor.execute(query)
		query_result=cursor.fetchall()
		return query_result

	## return the name of a given term
	# @param term
	def getTermName(self,term):	
		"""
		return the name of a given term
		@param term
		"""	
		cursor=self.db.cursor()
		query="select name from term where acc='%s'"%(term)
		cursor.execute(query)
		query_result=cursor.fetchall()
		return query_result[0][0]

    ## hypergeometric function for probability value 
    # @param a -- associated gene in input genesubset
    # @param b -- population (here use the offical home sapiens genes in NCBI)
    # @param b -- assoicated genes in population
    # @param d -- input genelist size
	def hypergeom(self,a,b,c,d):
		"""
		hypergeometric function for probability value 
    	@param a -- associated gene in input genesubset
    	@param b -- population (here use the offical home sapiens genes in NCBI)
    	@param b -- assoicated genes in population
   		@param d -- input genelist size
		"""
		prob = 1-hypergeom.cdf(a-1, b, c, d)
		return prob

	## Return list of list of genes in given input file
	def getGeneListsWithCSV(self):
		"""
		Return list of list of genes in given input file
		"""
		genelists=[]
		reader=csv.reader(open(self.inputfile,"rb"))
		for line in reader:
			genelist = line
			genelists.append(genelist)
		return genelists

	## Return the enrichment of given genelist
	# @param top pick the top GO term 
	def enrich_csv(self, top=None):
		"""
		Return the enrichment of given genelist
		@param top pick the top GO term 
		"""
		genelists=self.getGeneListsWithCSV()
		# genelists=self.mutualexclusive(lists)	
		cursor=self.db.cursor()
		##output is a list of lists, every list inside contains drivers and its associated term dictionary lists		
		for i in range(0,len(genelists)):
			print len(genelists[i])
			driver = genelists[i].pop(0)			
			with open(self.outputfile + '%s.csv'%(driver), 'wb') as csvfile:
				writer=csv.writer(csvfile, delimiter=',',)
				writer.writerow(["P-value", "GOTerm", "TermName", "SubGeneList"])
				#a list of term dictionaries, storing all term dictionaries assciated with specific driver
				genes = genelists[i]
				N = len(genes)
				# print genes
				term_n_dict={}
				term_gene_dict={}
				term_sigp_dict={}
				term_name_dict={}
				for item in genes:
					gene = "'"+item+"'"
					query = "SELECT distinct(acc) FROM final_symbol_term as FFS\
							where FFS.offsymbol=(%s)"%(gene)
					cursor.execute(query)
					query_result = cursor.fetchall()			
					for row in query_result:
						term=row[0]
						n=self.assocGeneNum(term)
						name=self.getTermName(term)
						term_n_dict[term]=n
						term_name_dict[term]=name;
						if term_gene_dict.has_key(term):
							term_gene_dict[term].append(item)
						else: 
							term_gene_dict[term]=[]
							term_gene_dict[term].append(item)
				# print term_gene_dict
				# keys_to_remove = [key for key, value in term_gene_dict.iteritems() if len(value)<=1]
				# for key in keys_to_remove:
				# 	del term_gene_dict[key]
		
				# print len(term_gene_dict)
				for key, value in term_gene_dict.items():
					x=len(value)
					# print x
					localn=term_n_dict[key]
					p=self.hypergeom(x,self.M,localn,N)
					if(p<self.threshold and x>1 ):
						term_sigp_dict[key]=p
				term_sigp_dict=collections.OrderedDict(sorted(term_sigp_dict.items(), key=lambda x: x[1]))
				if top is not None:
					for key in term_sigp_dict.keys()[:top]:
						output=[term_sigp_dict[key],key,term_name_dict[key],term_gene_dict[key]]
						writer.writerow(output)
				else:
					for key in term_sigp_dict.keys():
						output=[term_sigp_dict[key],key,term_name_dict[key],term_gene_dict[key]]
						writer.writerow(output)
				print "file saved"
		print "Done"		


# def main():
# 	tool = enrichment("localhost", "fanyu", "hellowork", "assocdb", "chunhui/MLL2-MLL3.targetgenes.v9.csv","chunhui/", 0.01)
# 	tool.enrich_csv();
# if __name__ == "__main__":
# 	main()	
