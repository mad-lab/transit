import sys

try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0

if hasWx:
    import wx.xrc
    from wx.lib.buttons import GenBitmapTextButton
    from pubsub import pub
    import wx.adv

import os
import time
import math
import random
import numpy
import scipy.stats
from scipy.stats import norm
import datetime
import operator

import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools




import io
from scipy.stats import hypergeom
import copy
from statsmodels.stats import multitest
# from datetime import datetime
import math

############# Description ##################

short_name = "pathway_enrichment"
long_name = "pathway_enrichment"
short_desc = "Gene set enrichment analysis"
long_desc = "Gene set enrichment analysis"
transposons = [] ##What's this for?
columns = ["[ID][descr]","Total genes","score","pval","padj","rank of genes"]

############# Analysis Method ##############

class GSEAAnalysis(base.TransitAnalysis):
	def __init__(self):		
		base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, GSEAMethod, GSEAGUI, [GSEAFile])


################## FILE ###################

class GSEAFile(base.TransitFile):

	def __init__(self):
		base.TransitFile.__init__(self, "#Example", columns)

	def getHeader(self, path):
		text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
		return text


################# GUI ##################

class GSEAGUI(base.AnalysisGUI):

	def __init__(self):
		base.AnalysisGUI.__init__(self)

########## METHOD #######################

class GSEAMethod(base.SingleConditionMethod):
	"""   
	Example
 
	"""
	def __init__(self,
				resamplingFile,				
				geneSetFile,
				output_file,p,N,M):

		base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, resamplingFile, geneSetFile, output_file, p,N,M)
		self.resamplingFile=resamplingFile
		self.geneSetFile = geneSetFile
		self.output = output_file
		self.p=p
		self.N=N
		self.M=M

####################################################################################
	###############FILES################################
	def loadD(self,fileName):
		file = open(fileName,"r")
		dict=[]
		ORFNameDict={}		
		for line in file:			
			if not(line.startswith("#")):				
				line=line.strip().split("\t") #this will return three elements				
				ORFNameDict[line[0]]=line[1]
				dict.append([line[0],float(line[10])])	
		file.close()
		return dict,ORFNameDict

		# file = open(fileName,"r")
		# line = file.readline()
		# while line.startswith("#"):
		# 	line = file.readline()		
		# dict=[]
		# line = line.split("\t")	
		# dict.append([line[0],float(line[10])])
		# ORFNameDict={line[0]:line[1]}
		# for line in file:		
		# 	line = line.strip().split("\t")
		# 	dict.append([line[0],float(line[10])])
		# 	ORFNameDict[line[0]]=line[1]
		# return dict,ORFNameDict

	def getM(self,protTable):
		file = open(protTable,"r")
		return len(file.read().splitlines())

	def loadGoTermsGoTermAsKey(self,fileName):
		dict={}
		descr={}
		file = open(fileName,"r")
		for f in file:
			if not(f.startswith("#")):				
				line = f.strip().split("\t") # It will return the id , description and list of ORFS				
				if len(line)!=3: 					
					self.output.write("Format Error in"+fileName+"\n")
				if not(line[0] in dict):
					dict[line[0]]=[]					
				dict[line[0]]+=line[2].split()
				descr[line[0]] = line[1]
		return dict,descr


	def saveInterestingPaths(self,fileName,m,n):
		inputF = open(fileName,"r")
		o = fileName.split(".")
		outputF = open(o[0]+"_"+str(m)+"_"+str(n)+"."+o[1],"w")
		for line in inputF:
			cad = line
			line = line.split(",")[1].split()
			lenLine = len(line)		
			if lenLine>m and lenLine<n:
				outputF.write(cad)
		outputF.close()
		

	################PREPROCESSING###################
	#Keeping only the value with p-value < 0.05
	#I receives a Dictionary with the gene as key
	#DEPRECATED!! D only has two columns
	def filteredPValue(self,D):	
		return {d:D[d] for d in D if D[d][4]<0.05}

	################STATISTICAL ANALYSIS############
	#M = the whole Genome
	#n = hits
	#N = sample size
	def hipergeometric(self,k,M,n,N):
		return hypergeom.sf(k,M,n,N)


	# HYPERGEOMETRIC

	def hyperGeometricTest(self,dict_k,M,dict_n,N):
		result={}
		for key in dict_n:
			n = len(dict_n[key])
			I = set(dict_k) & set(dict_n[key])
			k = len(I)
			result[key] = {"parameters":str(k)+"\t"+str(n),"p-value": hypergeom.sf(k,M,n,N), "Intersection":I}
		self.padjustForHypergeom(result)
		return result


	def padjustForHypergeom(self,result):
		keys = result.keys()
		pvalues = [result[k]["p-value"] for k in keys]
		b,adj=multitest.fdrcorrection(pvalues, alpha=0.05, method='indep')
		for i in range(len(keys)):
			result[keys[i]]["padjust"]=adj[i]



	def saveHyperGeometricTest(self,results,ORFNameDict,DESCR):
		#GoTermsDescription = loadGoTermsDescriptions("GO_terms_used_in_H37Rv.csv")
		# f = open(fileNameOut,"w")
		for k in results:
			cad = " ".join([x+"/"+ORFNameDict[x] for x in results[k]["Intersection"]])
			self.output.write(k+"\t"+DESCR[k]+"\t"+results[k]["parameters"]+"\t"+str(results[k]["p-value"])+"\t"+str(results[k]["padjust"])+"\t"+cad+"\n")

	# HYPERGEOMETRIC END

	#D is size N, it is a matrix
	#D[0]=RvName
	#D[1]=Correlation with C
	#S is a array of Genes, It might represent genes of a GoTerm, Pathway, etc
	def EnrichmentScore(self,D,S,p=1, sort=False):		
		L = copy.deepcopy(D)		
		rank = {}	
		if(sort==True):
			L = sorted(D,key = lambda x:x[1], reverse=False)	
		Ldict = {l[0]:l[1] for l in L}
		N_R=0.0000001
		N_R += sum([math.pow(abs(Ldict[l]),p) for l in Ldict if l in S]) #add abs according to the formula	
		lenL = len(L)
		miss = 1.0/(len(D)-len(S))	
		P=[[0,0]]
		ES=[]	
		#I remove abs from the formula because we're using LFCS as a metric
		if L[0][0] in S:
			P[0][0]=math.pow(abs(L[0][1]),p)/N_R
		else:
			P[0][1]=miss
		ESval = P[0][0]-P[0][1]
		ES.append(ESval)
		# print L[0][0],",",L[0][1],",",P[0][0],",",P[0][1],",",ESval
		for i in range(1,lenL):
			P.append([P[i-1][0],P[i-1][1],0])
			if L[i][0] in S:
				P[i][0]+=math.pow(abs(L[i][1]),p)/N_R
				rank[L[i][0]] = i
			else:
				P[i][1]+=miss
			ESTemp = P[i][0]-P[i][1]
			ES.append(ESTemp)
			# print L[i][0],",",L[i][1],",",P[i][0],",",P[i][1],",",ESTemp
			if ESval<ESTemp:
				ESval = ESTemp
		return ESval,N_R,rank

	#This is meant to save some time due to receiving precalculated values.
	def EnrichmentScoreOpt(self,D,S,miss,N_R,p=1):
		L = copy.deepcopy(D)	
		Ldict = {l[0]:l[1] for l in L}	
		lenL = len(L)	
		P=[[0,0]]
		ES=[]
		#Think about the abs: what is the estimator you're using	
		if L[0][0] in S:
			P[0][0]=math.pow(abs(L[0][1]),p)/N_R
		else:
			P[0][1]=miss	
		ESval = P[0][0]-P[0][1]
		ES.append(ESval)
		for i in range(1,lenL):
			P.append([P[i-1][0],P[i-1][1],0])		
			if L[i][0] in S:
				P[i][0]+=math.pow(abs(L[i][1]),p)/N_R
			else:
				P[i][1]+=miss
			ESTemp = P[i][0]-P[i][1]
			ES.append(ESTemp)
			if ESval<ESTemp:
				ESval = ESTemp
		return ESval

	def shuffleByColumnI(self,D,c):	
		lD= len(D)	
		D1 = [D[i][c] for i in range(lD)]
		random.shuffle(D1)
		for i in range(lD):
			D[i][c]=D1[i]
		return D


	def significance(self,D,S,p=1,N=1000):	
		ES_original=EnrichmentScore(D,S,p)
		pvalue=0.0
		for i in range(N):
			#Shuffle the LFCs
			D1 = shuffleByColumnI(D,0)		
			ES=EnrichmentScore(D,S,p)
			if ES_original<=ES:
				pvalue+=1.0
		return ES_original,pvalue/float(N)


	def GSEA(self,fileNameForD, fileNameForS, p=1,N=1000):
		#D is a list of 2XN
		D = self.loadD(fileNameForD)
		GoTerms,DESCR = self.loadGoTermsGoTermAsKey(fileNameForS)
		l=len(GoTerms.keys())
		Total = float(l*N)	
		gseaVal={}	
		for i in range(N):
			print "Percentage: ",(100.0*l*i)/Total,"%"
			D1 = self.shuffleByColumnI(D,0)
			D1 = sorted(D1,key = lambda x:x[1], reverse=False)
			for goTerm in GoTerms:
				if(i==0):
					gseaVal[goTerm]=[EnrichmentScore(D,GoTerms[goTerm]),0.0]
				ES=self.EnrichmentScore(D1,GoTerms[goTerm],p)
				if gseaVal[goTerm][0]<=ES:
					gseaVal[goTerm][1]+=1.0/N
		return gseaVal,GoTerms


	def padjust(self,gseaVal):
		keys = gseaVal.keys()
		pvals = [gseaVal[k][1] for k in keys]
		b,adj=multitest.fdrcorrection(pvals, alpha=0.05, method='indep')
		for i in range(len(keys)):
			gseaVal[keys[i]][2]=adj[i]


	def GSEA2(self,fileNameForD, fileNameForS, p=1,N=100):
		#D is a list of 2XN
		D,ORFNameDict = self.loadD(fileNameForD)
		GoTerms,DESCR = self.loadGoTermsGoTermAsKey(fileNameForS)
		l=len(GoTerms.keys())
		Total = float(l*N)	
		gseaVal={}
		N_R={}
		miss={}
		rank={}
		lenD=len(D)
		D1 = copy.deepcopy(D)
		self.progress_range(Total)
		for goTerm in GoTerms:
			miss[goTerm] = 1.0/(lenD-len(GoTerms[goTerm]))
			es,N_R[goTerm],rank[goTerm]=self.EnrichmentScore(D,GoTerms[goTerm],p,sort=True)				
			gseaVal[goTerm]=[es,0.0,0.0]
		for i in range(N):
			# Update Progress 
			pctg=(100.0*l*i)/Total
			text = "Running Pathway Enrichment Method... %5.1f%%" % (pctg)
			self.progress_update(text, l*i)			
			D1 = self.shuffleByColumnI(D1,0)
			D1 = sorted(D1,key = lambda x:x[1], reverse=False)		
			for goTerm in GoTerms:
				ES=self.EnrichmentScoreOpt(D1,GoTerms[goTerm],miss[goTerm],N_R[goTerm],p)
				if gseaVal[goTerm][0]<=ES:
					gseaVal[goTerm][1]+=1.0/N
		self.padjust(gseaVal)
		return gseaVal,GoTerms,rank,ORFNameDict,DESCR

	def saveExit(self,GSEADict, PathDict,rank,ORFNameDict,DESCR):		
		for gsea in GSEADict:
			# cad=gsea+":\n \t"+" ".join(PathDict[gsea])+"\n\tEnrichment Score: "+str(GSEADict[gsea][0])+" P-value:"+str(GSEADict[gsea][1])+" P-Adjust:"+str(GSEADict[gsea][2])+"\n"			
			d = rank[gsea]			
			sorted_d = sorted(d.items(), key=operator.itemgetter(1))
			rankCad = " ".join([str(k[1])+":"+k[0]+"/"+ORFNameDict[k[0]] for k in sorted_d])
			# cad="["+"][ ".join(gsea.split("-"))+"],"+str(len(PathDict[gsea]))+","+str(GSEADict[gsea][0])+","+str(GSEADict[gsea][1])+","+str(GSEADict[gsea][2])+","+rankCad+"\n"

			cad=gsea+"\t"+DESCR[gsea]+"\t"+str(len(PathDict[gsea]))+"\t"+str(GSEADict[gsea][0])+"\t"+str(GSEADict[gsea][1])+"\t"+str(GSEADict[gsea][2])+"\t"+rankCad+"\n"
			self.output.write(cad)
		self.output.close()


#################Z Test#################################
# D is the whole set with p values <=0.05
# S is the Sanger Categories with their genes.
	
	def Ztest(self,D,S):		
		auxD = [d[0] for d in D]		
		DictD={d[0]:d[1] for d in D}
		ES={}
		rank={}		
		for key in S: #GoTerms or Category in S			
			r=[DictD[k] for k in S[key] if k in DictD]
			lenr=len(r)
			if(lenr>1):				
				meanR = numpy.mean(r)
				es = math.sqrt(lenr)*meanR
				ES[key]= [lenr,es,norm.sf(es)]				
				rank[key]={i: auxD.index(i) for i in S[key] if i in auxD}
		self.padjustTest(ES)		
		return ES, rank

	def ChiSquareTest(self,D,S):		
		auxD = [d[0] for d in D]
		DictD={d[0]:d[1] for d in D}
		ES={}
		rank={}

		for key in S: #GoTerms or Category in S
			r=[DictD[k] for k in S[key] if k in DictD]
			lenr=len(r)
			if(lenr>19):
				meanR = numpy.mean(r)			
				lenSkey = lenr-1
				es=(sum([(ri-meanR)**2 for ri in r]) - (lenSkey)) / (2*lenSkey)
				ES[key]= [lenr,es, norm.sf(es)]
				rank[key]={i: auxD.index(i) for i in S[key] if i in auxD}
		self.padjustTest(ES)
		return ES,rank


	def printTest(self,test,ORFNameDict,rank,DESCR):		
		for key in test:
			cad = " ".join([str(rank[key][g])+":"+g+"/"+ORFNameDict[g] for g in rank[key]])
			self.output.write(key+"\t"+DESCR[key]+"\t"+str(test[key][0])+"\t"+str(test[key][1])+"\t"+str(test[key][2])+"\t"+str(test[key][3])+"\t"+cad+"\n")
		self.output.close()	

	def padjustTest(self,test):
		keys = test.keys()
		pvals = [test[k][2] for k in keys]
		b,adj=multitest.fdrcorrection(pvals, alpha=0.05, method='indep')
		for i in range(len(keys)):
			test[keys[i]]+=[adj[i]]

	def t_ishEstimators(self,D):
		lenD = len(D)
		for i in range(lenD):
			if D[i][1] == 0:
				D[i][1]=-4.0
			elif D[i][1] == 1:
				D[i][1]=3.0
			else:
				D[i][1]=norm.ppf(D[i][1])
####################################################################################
	@classmethod
	def fromGUI(self, wxobj):
		pass

	@classmethod
	def fromargs(self, rawargs): 
		(args, kwargs) = transit_tools.cleanargs(rawargs)

		p=int(kwargs.get("p", 1))
		N=int(kwargs.get("S", 1000))
		M = kwargs.get("M", "GSEA")
		resamplingFile = args[0]		
		geneSetFile = args[1]
		outpath = args[2]
		output = open(outpath, "w") #The outputfile is opened here!!!

		
		return self(resamplingFile,
				geneSetFile,
				output,
				p,N,M)

	def Run(self):
		self.transit_message("Starting Pathway Enrichment Method")
		start_time = time.time()

		if self.M =="GSEA":
			method = "GSEA"
			gseaVal,PathDict,rank,ORFNameDict,DESCR=self.GSEA2(self.resamplingFile,self.geneSetFile,self.p,self.N)
		elif self.M =="HYPE":
			method = "HYPERGEOMETRIC"
			D,ORFNameDict = self.loadD(self.resamplingFile)
			GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)
			DE = [d[0] for d in D if d[1]<=0.05]
			S = len(DE) #Sample Size in hyperGeometric	
			W = len(D) #Whole Genes			
			# k is the number of genes in DE that are in GoTermsWithRv
			results=self.hyperGeometricTest(DE,W,GoTermsWithRV,S)
		elif self.M=="GSEA-Z":
			method = "GSEA-Z"
			D,ORFNameDict = self.loadD(self.resamplingFile)
			GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)

			self.t_ishEstimators(D)
			results,rank = self.Ztest(D,GoTermsWithRV)
		elif self.M=="GSEA-CHI":
			method = "GSEA-CHI"
			D,ORFNameDict = self.loadD(self.resamplingFile)
			GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)
			self.t_ishEstimators(D)
			results,rank= self.ChiSquareTest(D,GoTermsWithRV)
		else:
			method = "Not a valid option"
			self.progress_update("Not a valid option", 100)

		self.output.write("#Pathway Enrichment\n")
		if self.wxobj:
			members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
			memberstr = ""
			for m in members:
				memberstr += "%s = %s, " % (m, getattr(self, m))
			self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata).encode('utf-8'), self.annotation_path.encode('utf-8'), self.output.name.encode('utf-8')))
		else:
			self.output.write("#Console: python %s\n" % " ".join(sys.argv))

		self.output.write("#Data: %s\n" % self.resamplingFile) 
		self.output.write("#Annotation path: %s\n" % self.annotation_path.encode('utf-8')) 
		self.output.write("#Time: %s\n" % (time.time() - start_time))
		self.output.write("#Methodology: %s\n"%method)
		if self.M =="GSEA":
			columns = ["[ID][descr]","Total genes","score","pval","padj","rank of genes"]
			self.output.write("#%s\n" % "\t".join(columns))
			self.saveExit(gseaVal,PathDict,rank,ORFNameDict,DESCR)
		elif self.M =="HYPE":
			columns=["cat id descr","Total genes","Total in intersection","pval","padj","genes in intersection"]
			self.output.write("#%s\n" % "\t".join(columns))
			self.saveHyperGeometricTest(results,ORFNameDict,DESCR)
		elif self.M =="GSEA-Z":
			columns=["#ID-Description","Total Genes","Score","P-Value","P-Adjust","genes"]
			self.output.write("#%s\n" % "\t".join(columns))
			self.printTest(results,ORFNameDict,rank,DESCR)
		elif self.M =="GSEA-CHI":			
			columns=["#ID-Description","Total Genes","Score","P-Value","P-Adjust","genes"]
			self.output.write("#%s\n" % "\t".join(columns))
			self.printTest(results,ORFNameDict,rank,DESCR)
		self.transit_message("") # Printing empty line to flush stdout 
		self.transit_message("Adding File: %s" % (self.output.name))
		self.add_file(filetype="Pathway Enrichment")
		self.finish()
		self.transit_message("Finished Pathway Enrichment Method") 


	@classmethod
	def usage_string(self):
		return """python %s pathway_enrichment <resampling files> <annotation file> <output file> [-p .-S -M < GSEA, HYPE, Z, CHI >] """ % (sys.argv[0])


if __name__ == "__main__":

	(args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

	print "ARGS:", args
	print "KWARGS:", kwargs

	G = GSEAMethod.fromargs(sys.argv[1:])

	print G
	G.console_message("Printing the member variables:")   
	G.print_members()

	print ""
	print "Running:"

	G.Run()



