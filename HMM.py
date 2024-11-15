import numpy as np
import pandas as pd

import math


#number of states
N=2 
cover_cut=0.8


#Ti: Introgression of Nd
def initA(cut,RR,Ti, a) -> np.array:
    A = np.zeros((2,2))
    
    A[0][1]=Ti*RR*cut*a
    A[0][0]=1-A[0][1]
 
    A[1][0]=Ti*RR*cut*(1-a)
    A[1][1]=1-A[1][0]
    
    return A

#Ti: Introgression of Nd
#Taf: Time out of Africa
#Tn: Time of Split between Nd and Sapiens

def initB(m,cut, lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(2,n_st,n_st))
    meani = lmbd[0]
    meann = lmbd[1]
    meanaf = lmbd[2]
    
    Pi = np.empty(n_st)
    Paf=np.empty(n_st)
    Pn=np.empty(n_st)

    
    Pi[0]=np.exp(-meani)
    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    
    sumi=0
    sumaf=0
    sumn=0
    
    for i in range(1,n_st):
        Pi[i]=Pi[i-1]*meani/i
        Paf[i]=Paf[i-1]*meanaf/i
        Pn[i]=Pn[i-1]*meann/i
        
        sumi=sumi+Pi[i]
        sumaf=sumaf+Paf[i]
        sumn=sumn+Pn[i]

    Pi[0]=1-sumi
    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st): 
        for j in range(n_st):
            B[0][i][j]=Paf[i]*Pn[j]
            B[1][i][j]=Pn[i]*Pi[j]
                  
    return B
    
def initB_arch_cover(m,cut, lmbd, n_st, cover):
    
    B = np.empty(shape=(2,n_st,n_st))
    meani, meann, meanaf, meani2 = lmbd[0], lmbd[1]*cover, lmbd[2], lmbd[0] * cover
    Pi, Paf, Pn, Pi2 = [np.empty(n_st) for i in range(4)]
    Pi[0], Paf[0], Pn[0], Pi2[0] = np.exp(-meani), np.exp(-meanaf), np.exp(-meann), np.exp(-meani2)
    sumi, sumaf, sumn, sumi2 = [0 for i in range(4)]

    
    for i in range(1,n_st):
        Pi[i], Paf[i], Pn[i], Pi2[i] = Pi[i-1]*meani/i, Paf[i-1]*meanaf/i, Pn[i-1]*meann/i, Pi2[i-1]*meani2/i

        
        sumi += Pi[i]
        sumaf += Paf[i]
        sumn += Pn[i]
        sumi2 += Pi2[i]

    Pi[0], Paf[0], Pn[0], Pi2[0] = 1-sumi, 1-sumaf, 1-sumn, 1-sumi2
    
    for i in range(n_st): 
        for j in range(n_st):
            B[0][i][j]=Paf[i]*Pn[j]
            B[1][i][j]=Pn[i]*Pi2[j]
                  
    return B
        
    
    

# Emission probabilities matrix (Skov'18) for the segmets without Neanderthal cover

#Ti: Introgression of Nd
#Taf: Time out of Africa
#Tn: Time of Split between Nd and Sapiens

def initBwN(cut, lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(2,n_st))
    meann = lmbd[1]
    meanaf = lmbd[2]


    Paf=np.empty(n_st)
    Pn=np.empty(n_st)

    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    

    sumaf=0
    sumn=0
    
    for i in range(1,n_st):
        Paf[i]=Paf[i-1]*meanaf/i
        Pn[i]=Pn[i-1]*meann/i
        
        sumaf=sumaf+Paf[i]
        sumn=sumn+Pn[i]

    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st): 
        B[0][i]=Paf[i]
        B[1][i]=Pn[i]
               

    return B
    
    
    
def viterbi_modified(V, initial_distribution, a, b_our_mas, b_Skov, archaic_cover):
    
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))

    if archaic_cover[0]>cover_cut:
        omega[0, :] = np.log(initial_distribution * b_our_mas[int(archaic_cover[0]*10),:, V[0][0],V[0][1]])

    else:
        
        omega[0, :] = np.log(initial_distribution * b_Skov[:, V[0][0]])        
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            if archaic_cover[t]>cover_cut:
                probability = omega[t - 1] + np.log(a[:, j]) + np.log(b_our_mas[int(archaic_cover[t]*10),j, V[t][0], V[t][1]])
                
            else:
                probability = omega[t - 1] + np.log(a[:, j]) + np.log(b_Skov[j, V[t][0]])
                
               
 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)

    return result
    




def get_HMM_tracts(seq):
    migrating_tracts = []
    for i in range(N):
        migrating_tracts.append([])
    start=0
    for i in range(1,len(seq)):
        if seq[i]!=seq[i-1]:
            migrating_tracts[seq[i-1]].append([start,i-1])
            start=i
    migrating_tracts[seq[len(seq)-1]].append([start,len(seq)-1])
    return migrating_tracts








def viterbi_modified2(V, initial_distribution, a, b_our_mas, b_Skov, archaic_cover):
    
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))


        
    omega[0, :] = np.log(initial_distribution * b_Skov[:, V[0][0]])        
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability

            probability = omega[t - 1] + np.log(a[:, j]) + np.log(b_Skov[j, V[t][0]])

               
 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)

    return result

import EM

def posterior(o, initial_distribution, a, b_our_mas, b_Skov, archaic_cover, gaps, cut_off):
    
    alpha, sc_factors = EM.alpha_scaled_opt_gaps(a,b_Skov, b_our_mas, o, initial_distribution, gaps, archaic_cover)
    beta = EM.beta_scaled_opt_gaps(a, b_Skov, b_our_mas, o, sc_factors, gaps, archaic_cover)    
    gamma = np.array(EM.def_gamma_gaps(alpha, beta, gaps)) 


    result=[]
    for j in range(len(gamma[0])):
        if gamma[:,j].argmax()==1 and gamma[:,j].max()>=cut_off:
            result.append(1)
        else:
            result.append(0)


    return result








