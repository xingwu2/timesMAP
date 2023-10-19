#import modules
import time
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import re

import multiplicative_gibbs as multi

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
alpha_trace.to_csv(args.output+"_multiplicative_alpha_trace.txt",sep="\t",header=False,index=False)
beta_trace.to_csv(args.output+"_multiplicative_beta_trace.txt",sep="\t",header=False,index=False)
gamma_trace.to_csv(args.output+"_multiplicative_gamma_trace.txt",sep="\t",header=False,index=False)


