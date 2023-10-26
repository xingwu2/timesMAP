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
	trace_avg,trace_sd,alpha_avg,alpha_sd,beta_avg,beta_sd,pip = multi.sampling(y,C,X,1,1,0.001,12000,args.output)
	
	OUTPUT_TRACE = open(args.output+"_multiplicative_param.txt","w")
	for i in range(len(trace_avg)):
		print("%f\t%f" %(trace_avg[i],trace_sd[i]),file = OUTPUT_TRACE)
	
	OUTPUT_ALPHA = open(args.output+"_multiplicative_alpha.txt","w")
	for i in range(len(alpha_avg)):
		print("%f\t%f" %(alpha_avg[i],alpha_sd[i]),file = OUTPUT_ALPHA)

	OUTPUT_BETA = open(args.output+"_multiplicative_beta.txt","w")
	for i in range(len(beta_avg)):
		print("%s\t%f\t%f" %(i,beta_avg[i],beta_sd[i]),file = OUTPUT_BETA)

	OUTPUT_PIP = open(args.output+"_multiplicative_pip.txt","w")
	for i in range(len(pip)):
		print("%s\t%f" %(i,pip[i]),file = OUTPUT_PIP)

elif args.mode == 2:
	trace_avg,trace_sd,alpha_avg,alpha_sd,beta_avg,beta_sd,pip = add.sampling(y,C,X,1,1,0.001,12000,args.output)
	
	OUTPUT_TRACE = open(args.output+"_additive_param.txt","w")
	for i in range(len(trace_avg)):
		print("%f\t%f" %(trace_avg[i],trace_sd[i]),file = OUTPUT_TRACE)

	OUTPUT_ALPHA = open(args.output+"_additive_alpha.txt","w")
	for i in range(len(alpha_avg)):
		print("%f\t%f" %(alpha_avg[i],alpha_sd[i]),file = OUTPUT_ALPHA)


	OUTPUT_BETA = open(args.output+"_additive_beta.txt","w")
	for i in range(len(beta_avg)):
		print("%s\t%f\t%f" %(i,beta_avg[i],beta_sd[i]),file = OUTPUT_BETA)

	OUTPUT_PIP = open(args.output+"_additive_pip.txt","w")
	for i in range(len(pip)):
		print("%s\t%f" %(i,pip[i]),file = OUTPUT_PIP)


