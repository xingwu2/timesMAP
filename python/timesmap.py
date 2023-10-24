#import modules
import time
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import re

import multiplicative_gibbs as multi
import additive_gibbs as add

parser = argparse.ArgumentParser()
parser.add_argument('-g',type = str, action = 'store', dest = 'geno')
parser.add_argument('-c',type = str, action = 'store', dest = 'covar')
parser.add_argument('-y',type = str, action = 'store', dest = 'pheno')
parser.add_argument('-m',type = int, action = 'store', dest = 'mode',default = "1",help="1: multiplicative;	2: additive")
parser.add_argument('-o',type = str, action = 'store', dest = 'output',help = "the prefix of output files")
args = parser.parse_args()


X = np.loadtxt(args.geno,delimiter="\t")

y = []
with open(args.pheno) as f:
	for line in f:
		line = line.strip("\n")
		y.append(float(line))

y = np.asarray(y)

C =  np.array(pd.read_csv(args.covar,sep="\t",header=None)) 

if args.mode == 1:
	trace,alpha_trace,beta_trace,gamma_trace = multi.sampling(y,C,X,1,1,0.001,12000,args.output)
	trace.to_csv(args.output+"_multiplicative_trace.txt",sep="\t",header=False,index=False)
	
	alpha_trace_avg = np.mean(alpha_trace,axis=0)
	alpha_trace_sd = np.std(alpha_trace,axis = 0)
	OUTPUT_ALPHA = open(args.output+"_multiplicative_alpha.txt","w")
	for i in range(C.shape[1]):
		print("%f\t%f" %(alpha_trace_avg[i],alpha_trace_sd[i]),file = OUTPUT_ALPHA)

	beta_trace_avg = np.mean(beta_trace,axis=0)
	beta_trace_sd = np.std(beta_trace,axis = 0)
	OUTPUT_BETA = open(args.output+"_multiplicative_beta.txt","w")
	for i in range(len(beta_trace_avg)):
		print("%s\t%f\t%f" %(i,beta_trace_avg[i],beta_trace_sd[i]),file = OUTPUT_BETA)

	pip = np.mean(gamma_trace,axis = 0)
	OUTPUT_PIP = open(args.output+"_multiplicative_pip.txt","w")
	for i in range(len(pip)):
		print("%s\t%f" %(i,pip[i]),file = OUTPUT_PIP)

elif args.mode == 2:
	trace,alpha_trace,beta_trace,gamma_trace = add.sampling(y,C,X,1,1,0.001,12000,args.output)
	trace.to_csv(args.output+"_additive_trace.txt",sep="\t",header=False,index=False)

	alpha_trace_avg = np.mean(alpha_trace,axis=0)
	alpha_trace_sd = np.std(alpha_trace,axis = 0)
	OUTPUT_ALPHA = open(args.output+"_additive_alpha.txt","w")
	for i in range(C.shape[1]):
		print("%f\t%f" %(alpha_trace_avg[i],alpha_trace_sd[i]),file = OUTPUT_ALPHA)

	beta_trace_avg = np.mean(beta_trace,axis=0)
	beta_trace_sd = np.std(beta_trace,axis = 0)
	OUTPUT_BETA = open(args.output+"_additive_beta.txt","w")
	for i in range(len(beta_trace_avg)):
		print("%s\t%f\t%f" %(i,beta_trace_avg[i],beta_trace_sd[i]),file = OUTPUT_BETA)

	pip = np.mean(gamma_trace,axis = 0)
	OUTPUT_PIP = open(args.output+"_additive_pip.txt","w")
	for i in range(len(pip)):
		print("%s\t%f" %(i,pip[i]),file = OUTPUT_PIP)


