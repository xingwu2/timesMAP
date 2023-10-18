import numpy as np
import scipy as sp
import math
from sklearn import preprocessing
import pandas as pd
import time
import sys
#import pybind11
#import pymc3 as pm3

#import multiplicative

def circle_product_matrix(X, beta):
    # Element-wise multiplication of each row of X with beta.
    x_beta = np.multiply(X, beta)
    # Compute the product along the second axis (columns) to get the product for each row.
    y = np.prod(x_beta + 1, axis=1)
    return(y)

def circle_prodcut_vector(x,beta):
	y = x * beta + 1
	return(y)

def sample_pie(gamma,pie_a,pie_b):
	a_new = np.sum(gamma)+pie_a
	b_new = np.sum(1-gamma)+pie_b
	pie_new = np.random.beta(a_new,b_new)
	return(pie_new)

def sample_sigma_1(beta,gamma,a_sigma,b_sigma):
	a_new = 0.5*np.sum(gamma)+a_sigma
	b_new = 0.5*np.sum(np.multiply(np.square(beta),gamma))+b_sigma
	sigma_1_neg2 =np.random.gamma(a_new,1.0/b_new)
	sigma_1_new = math.sqrt(1/sigma_1_neg2)
	return(sigma_1_new)

def sample_sigma_e(y,H,beta,C,alpha,a_e,b_e):
	n = len(y)
	a_new = float(n)/2+a_e
	resid = y - circle_product_matrix(H,beta) - np.matmul(C,alpha)
	b_new = np.sum(np.square(resid))/2+b_e
	sigma_e_neg2 =np.random.gamma(a_new,1.0/b_new)
	sigma_e_new = np.sqrt(1/sigma_e_neg2)
	return(sigma_e_new)


def sample_alpha(y,H,beta,C,alpha,sigma_e):

	r,c = C.shape
	H_beta = circle_product_matrix(H,beta)
	if c == 1:
		new_variance = 1/(np.sum(C[:,0]**2) * sigma_e**-2)
		new_mean = new_variance*np.dot((y-H_beta),C[:,0])*sigma_e**-2
		new_alpha = np.random.normal(new_mean,math.sqrt(new_variance))
	else:
		C_alpha = np.matmul(C,alpha)
		for i in range(c):
			new_variance = 1/(np.sum(C[:,i]**2) * sigma_e**-2)
			C_alpha_negi = C_alpha - C[:,i] * alpha[i]
			new_mean = new_variance*np.dot(y-C_alpha_negi-H_beta,C[:,i])*sigma_e**-2
			alpha[i] = np.random.normal(new_mean,math.sqrt(new_variance))
			C_alpha = C_alpha_negi + C[:,i] * alpha[i]
	return(alpha)


def sample_gamma(y,C,alpha,H,beta,pie,sigma_1,sigma_e,gamma):
	C_alpha = np.matmul(C,alpha)
	sigma_e_neg2 = sigma_e**-2
	sigma_1_neg2 = sigma_1**-2
	H_beta = circle_product_matrix(H,beta)
	H_beta_neg = np.transpose(np.ones((len(beta),len(y))) * H_beta) / (H*beta + 1)
	H_beta_neg_H =  np.multiply(H ,H_beta_neg)
	variance = 1/(np.sum((H_beta_neg_H)**2,axis=0)*sigma_e_neg2+sigma_1_neg2)
	f = np.sqrt(1/(np.sum(H_beta_neg_H**2,axis=0) * sigma_1**2 * sigma_e_neg2 + 1))
	gamma_0_pie = np.zeros(len(beta))
	for i in range(len(beta)):
		residual = y - C_alpha - H_beta_neg[:,i]
		mean = variance[i]*np.dot(residual,H_beta_neg_H[:,i])*sigma_e_neg2
		A = f[i] * np.exp(0.5*mean**2/variance[i])
		gamma_0_pie[i] = (1-pie) / ((1-pie)+pie*A)
	gamma = np.random.binomial(1,1-gamma_0_pie)
	return(gamma)

def sample_beta(y,C,alpha,H,beta,gamma,sigma_1,sigma_e):

	C_alpha = np.matmul(C,alpha)
	sigma_e_neg2 = sigma_e**-2
	sigma_1_neg2 = sigma_1**-2
	H_beta = circle_product_matrix(H,beta)

	for i in range(len(beta)):
		if gamma[i] == 0:
			H_beta = np.divide(H_beta,circle_prodcut_vector(H[:,i], beta[i]))
			beta[i] = 0
		else:
			H_beta_negi = np.divide(H_beta,circle_prodcut_vector(H[:,i], beta[i]))
			new_variance = 1/(np.sum((H_beta_negi*H[:,i])**2) * sigma_e_neg2+sigma_1_neg2)
			residual = y - C_alpha -  H_beta_negi
			new_mean = new_variance*np.dot(H_beta_negi*H[:,i],residual)*sigma_e_neg2
			beta[i] = np.random.normal(new_mean,math.sqrt(new_variance))
			H_beta = H_beta_negi * circle_prodcut_vector(H[:,i], beta[i])
	return(beta)


def sampling(y,C,HapDM,sig1_initiate,sige_initiate,pie_initiate,step_size,iters,prefix):


	LOG = open(prefix+".log","w")
	##specify hyper parameters
	pie_a = 1
	pie_b = 1
	a_sigma = 1
	b_sigma = 1
	a_e = 1
	b_e = 1

	sigma_1 = sig1_initiate
	sigma_e = sige_initiate
	pie = pie_initiate


	print("initiate:",sigma_1,sigma_e,pie,file = LOG)

	#initiate beta,gamma and H matrix
	C_r, C_c = C.shape

	H = np.array(HapDM)

	#H = preprocessing.scale(H,with_std=False)
	print(H[:5,:5])
	H_r,H_c = H.shape

	#initiate alpha, alpha_trace, beta_trace and gamma_trace

	it = 0
	burn_in_iter = 2000
	trace = np.empty((iters-2000,7))
	alpha_trace = np.empty((iters-2000,C_c))
	gamma_trace = np.empty((iters-2000,H_c))
	beta_trace = np.empty((iters-2000,H_c))
	top5_beta_trace = np.empty((iters-2000,5))

	alpha = np.random.random(size = C_c)
	gamma = np.random.binomial(1,pie,H_c)
	beta = np.array(np.zeros(H_c))

	for i in range(H_c):
		if gamma[i] == 0:
			beta[i] = 0
		else:
			beta[i] = np.random.normal(0,sigma_1)

	#start sampling

	while it < iters:
		before = time.time()
		beta_pre = np.array(beta)
		gamma_pre = np.array(gamma)
		alpha_pre = np.array(alpha)
		sigma_1_pre = sigma_1
		sigma_e_pre = sigma_e
		pie_pre = pie

		sigma_1_update = sample_sigma_1(beta_pre,gamma_pre,a_sigma,b_sigma)

		if sigma_1_update < 0.05:
			sigma_1_update = 0.05
			pie_update = 0
		else:
			pie_update = sample_pie(gamma_pre,pie_a,pie_b)
		sigma_e_update = sample_sigma_e(y,H,beta_pre,C,alpha_pre,a_e,b_e)
		gamma_update = sample_gamma(y,C,alpha,H,beta,pie_update,sigma_1_update,sigma_e_update,gamma)
		alpha_update = sample_alpha(y,H,beta_pre,C,alpha_pre,sigma_e_update)
		beta_update = sample_beta(y,C,alpha_update,H,beta_pre,gamma_update,sigma_1_update,sigma_e_update)
		after = time.time()

		genetic_var = np.var(circle_product_matrix(H,beta_update))
		pheno_var = np.var(y - np.matmul(C,alpha_update))
		large_beta = np.absolute(beta_update) > 0.3
		large_beta_ratio = np.sum(large_beta) / len(beta_update)
		large_beta_heritability = np.var(circle_product_matrix(H[:,large_beta],beta_update[large_beta])) / pheno_var
		total_heritability = genetic_var / pheno_var

		if it > 100 and large_beta_heritability > 1 and large_beta_heritability > total_heritability:
			print("unrealistic beta sample",it,genetic_var,pheno_var,large_beta_heritability,total_heritability)
			continue

		# elif it > 100 and large_beta_heritability > total_heritability:
		#	print("unrealistic beta sample",it,genetic_var,pheno_var,large_beta_heritability,total_heritability)
		# 	continue

		else:
			beta = np.array(beta_update)
			gamma = np.array(gamma_update)
			alpha = np.array(alpha_update)
			sigma_1 = sigma_1_update
			sigma_e = sigma_e_update
			pie = pie_update

			print(it,str(after - before),pie,large_beta_ratio,large_beta_heritability,total_heritability)

			if it >= burn_in_iter:
				trace[it-burn_in_iter,:] = [sigma_1,sigma_e,large_beta_ratio,large_beta_heritability,total_heritability,pie,it]
				gamma_trace[it-burn_in_iter,:] = gamma
				beta_trace[it-burn_in_iter,:] = beta
				alpha_trace[it-burn_in_iter,:] = alpha
				top5_beta_trace[it-burn_in_iter,:] = np.sort(np.absolute(beta))[::-1][:5]

			if it >= burn_in_iter + 9999: # after burn-in iterations, test convergence

				max_z = []

				for a in range(C_c):
					after_burnin_alpha = alpha_trace[:,a]
					alpha_zscores = pm3.geweke(after_burnin_alpha)[:,1]
					max_z.append(np.amax(np.absolute(alpha_zscores)))

				for b in range(5):
					after_burnin_beta = top5_beta_trace[:,b]
					beta_zscores = pm3.geweke(after_burnin_beta)[:,1]
					max_z.append(np.amax(np.absolute(beta_zscores)))

				#convergence for pie
				after_burnin_pie = trace[:,2]
				pie_zscores = pm3.geweke(after_burnin_pie)[:,1]
				max_z.append(np.amax(np.absolute(pie_zscores)))

				#convergence for large_heritability
				after_burnin_var = trace[:,3]
				var_zscores = pm3.geweke(after_burnin_var)[:,1]
				max_z.append(np.amax(np.absolute(var_zscores)))

				#convergence for sigma_1
				after_burnin_sigma1 = trace[:,0]
				sigma1_zscores = pm3.geweke(after_burnin_sigma1)[:,1]
				max_z.append(np.amax(np.absolute(sigma1_zscores)))

				#convergence for sigma_e
				after_burnin_sigmae = trace[:,1]
				sigmae_zscores = pm3.geweke(after_burnin_sigmae)[:,1]
				max_z.append(np.amax(np.absolute(sigmae_zscores)))

				if  np.amax(max_z) < 1.5:
					print("convergence has been reached at %i iterations." %(it),file=LOG)
					break

				else:
					trace_ = np.empty((1000,7))
					gamma_trace_ = np.empty((1000,H_c))
					beta_trace_ = np.empty((1000,H_c))
					alpha_trace_ = np.empty((1000,C_c))
					top5_beta_trace_ = np.empty((1000,5))

					trace = np.concatenate((trace[-(iters - burn_in_iter-1000):,:],trace_),axis=0)
					gamma_trace = np.concatenate((gamma_trace[-(iters - burn_in_iter-1000):,:],gamma_trace_),axis=0)
					beta_trace = np.concatenate((beta_trace[-(iters - burn_in_iter-1000):,:],beta_trace_),axis=0)
					alpha_trace = np.concatenate((alpha_trace[-(iters - burn_in_iter-1000):,:],alpha_trace_),axis=0)
					top5_beta_trace = np.concatenate((top5_beta_trace[-(iters - burn_in_iter-1000):,:],top5_beta_trace_),axis = 0)

					burn_in_iter += 1000
					iters += 1000

			if (it - burn_in_iter) >= 0 and (it - burn_in_iter ) % 1000 == 0:
				print("%i iterations have sampled" %(it), str(after - before),trace[it-burn_in_iter,:],file=LOG)

			it += 1

	LOG.close()
	trace = pd.DataFrame(trace)
	alpha_trace = pd.DataFrame(alpha_trace)
	beta_trace = pd.DataFrame(beta_trace)
	gamma_trace = pd.DataFrame(gamma_trace)

	return(trace,alpha_trace,beta_trace,gamma_trace)



HapDM = np.loadtxt("large_X.txt",delimiter="\t")
print(HapDM.shape)

y = []
with open("large_multiplicative_y.txt","r") as f:
	for line in f:
		line = line.strip("\n")
		y.append(float(line))

y = np.asarray(y)

C =  np.array(pd.read_csv("large_C.txt",sep="\t",header=None))
print(C.shape)




trace,alpha_trace,beta_trace,gamma_trace = sampling(y,C,HapDM,0.5,0.5,0.5,700,12000,"test")
trace.to_csv("test_multiplicative_trace.txt",sep="\t",header=False,index=False)
alpha_trace.to_csv("test_multiplicative_alpha_trace.txt",sep="\t",header=False,index=False)
beta_trace.to_csv("test_multiplicative_beta_trace.txt",sep="\t",header=False,index=False)
gamma_trace.to_csv("test_multiplicative_gamma_trace.txt",sep="\t",header=False,index=False)
