import HMM
import numpy as np
from numpy import linalg as LNG 
import math
import math
from numpy import linalg as LNG 

import useful as usfl

    

N=2 
# each window is covered by cover[t] parameter which in [-0.001...0.999]. Cover_cut remove windows with low covering and transforms method to Skov. 
cover_cut = 0.8 

#transitions and emissions in case of gaps in observable chromosome (centromeres/telomeres), not archaic gaps. 
a_gaps=np.identity(N)
b_gaps=np.ones((N,1))    
    
#EM-common+gaps
 # forward-algo
def alpha_scaled_opt_gaps(a,b, b_mas, o, p, gaps, cover):   
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, len(o)))
    if cover[0]>cover_cut:
        alpha[:,0] = b_mas[int(cover[0]*10),: , o[0][0],o[0][1]] * p
    else:
    
        alpha[:, 0] = b[:, o[0][0]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    

    for t in range(1, len(o)):   
        if usfl.point_in_set(t, gaps)==True:

            for i in range(0, N):
                
                alpha[i, t] = np.dot(alpha[:, t-1],a_gaps[:,i]) * b_gaps[i][0]     

        else:
            if cover[t] > cover_cut:
        
                for i in range(0, N):            
                    alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b_mas[int(cover[t]*10),i, o[t][0],o[t][1]] 
            else:
                for i in range(0, N):            
                    alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt_gaps(a,b, b_mas, o, scaling_factors, gaps, cover):
    
    beta = np.zeros((N, len(o)))
    
    length = len(o)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(len(o)-2,-1,-1):        
        if usfl.point_in_set(t, gaps)==True:
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a_gaps[i, l] * b_gaps[l][0] * beta[l, t+1]       
        else:
            if cover[t+1] > cover_cut:
               
                for i in range(0, N):             
                    for l in range(0, N):
                        beta[i, t] += a[i, l] * b_mas[int(cover[t+1]*10),l, o[t+1][0],o[t+1][1]] * beta[l, t+1]
                 
            else:
                for i in range(0, N):             
                    for l in range(0, N):
                        beta[i, t] += a[i, l] * b[l, o[t+1][0]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma_gaps(alpha, beta, gaps):
        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])

        for i in range(0,N):        
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi_gaps( a, b, o, alpha, beta, gaps):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        if point_in_set(t, gaps)==True:
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a_gaps[i, j] * b_gaps[j][0] * beta[j, t+1]        
            for i in range(0, N):
                for j in range(0, N):       
                    ksi[i, j, t] = (alpha[i, t]*a_gaps[i, j]*b_gaps[j,0]* beta[j, t+1]) / denom
        else:
        
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]
                
        
            for i in range(0, N):
                for j in range(0, N):
                    ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]) / denom
    
    return ksi

    
def new_lambda_i_gaps(O, Gamma, gaps, cover):
    nom, denom = 0, 0
    i=0
    
    for o in O:
        gamma=Gamma[i]
        i+=1
        
        for t in range(0, len(o), 1):
            if usfl.point_in_set(t, gaps)==False:
                if cover[t] > cover_cut:   
                    nom += o[t, 1] *  gamma[1, t]
                    denom += (int(cover[t]*10)/10+0.1)*gamma[1,t]


    if denom==0: # only in case of no archaic covering
        return 0
    return nom/denom
   
def new_lambda_n_gaps(O, Gamma, gaps, cover):
    nom, denom = 0, 0
    i=0
    for o in O:    
        gamma=Gamma[i]
        i+=1        

        for t in range(0, len(o), 1):   
            if usfl.point_in_set(t, gaps)==False:
                if cover[t] > cover_cut:
                    nom += o[t, 1] *  gamma[0, t]  + o[t, 0]  * gamma[1, t]
                    denom += (int(cover[t]*10)/10+0.1) * gamma[0, t] + gamma[1, t]   
                else:
                    nom +=  o[t, 0]  * gamma[1, t]
                    denom += gamma[1, t]                    
    
    return nom/ denom
          
def new_lambda_af_gaps(O, Gamma, gaps):
    nom, denom = 0, 0
    i=0
    
    for o in O:    
        gamma=Gamma[i]
        i+=1   
        for t in range(0, len(o), 1):   
            if usfl.point_in_set(t, gaps)==False:
                nom += o[t, 0] * gamma[0, t]
                denom += gamma[0, t]  
    
    return nom/ denom
          
    
    
    
    
def E_step_gaps(cut,  p, O, n_states, mu,rr, lambda_old, gaps, cover):

   
    b_Skov = HMM.initBwN(cut, lambda_old[0:3], n_states)
    a = HMM.initA(cut,rr, lambda_old[4]/(mu*cut), lambda_old[3])
    b_our_mas = np.array([HMM.initB_arch_cover(mu,cut, lambda_old, n_states, 0.1+i*0.1) for i in range(10)])
    

    
    GAMMA=[]
    for o in O:
        alpha, sc_factors = alpha_scaled_opt_gaps(a,b_Skov, b_our_mas, o, p, gaps, cover)
        beta = beta_scaled_opt_gaps(a, b_Skov, b_our_mas, o, sc_factors, gaps, cover)    
        gamma = def_gamma_gaps(alpha, beta, gaps)    

        GAMMA.append(gamma)       

    return new_lambda_i_gaps(O, GAMMA, gaps, cover), new_lambda_n_gaps(O, GAMMA, gaps, cover), new_lambda_af_gaps(O, GAMMA, gaps), lambda_old[3], lambda_old[4]   
    
    
    
    
def EM_algorithm_gaps(p,o_mas, n_states, mu, rr, lambda_0, epsilon, cut, n_em_steps, gaps, cover ):
    
    lmbd = np.array(lambda_0)
    
    em_steps = 0
    for i in range(n_em_steps):
        


        lmbd_new = np.array(E_step_gaps(cut,  p, o_mas, n_states, mu, rr, lmbd,gaps, cover))
        if lmbd_new[4]>0.1:
            print('Oops. something went wrong with sample')
            return(lambda_0)
        em_steps += 1

        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new
        

    print('Число шагов в EM -алгоритме', em_steps )
    return lmbd_new
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
