import numpy as np

from functools import reduce
from operator import xor
from itertools import chain
import pandas as pd


# two set of intervals intersection 
def intersections(a,b):
    if b==[[]]:
        return a
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges
    
    


def exclude_gaps(set1, set2): #set2 is set of gaps [[],[],[]]

    l = sorted((reduce(xor, map(set, chain(set1 , set2)))))
    XOR=[l[i:i + 2] for i in range(0, len(l), 2)]

    return intersections(XOR, set1)
    
    

def point_in_set(p, m_a):
    if m_a==[[]]:
        return False
    def point_in(p, a):
        if p>=a[0] and p<=a[1]:
            return True
        else:
            return False
    f=False
    for j in m_a:
        f=point_in(p,j)
        if f==True:
            return f
    return f



def read_bed(f):


    with open(f,'r') as f:
        domain=f.readlines()
    domain=[domain[i].replace('\n','').split('\t')[1:3] for i in range(len(domain))]    

    for i in  range(len(domain)):
        domain[i][0]=int(domain[i][0])
        domain[i][1] = int(domain[i][1])
    return domain
    
    
def read_aa_file(f_aa): #read file with ancestral alleles
    with open(f_aa,'r') as f:
        l=f.readlines()
    dct={}
    for  i in l:    
        m=i.replace('\n','').split('\t')
        dct[int(m[1])]=m[2]
    return dct




def main_read(file, f_aa): # return the dictionary with REF,ALT, Outgroup and Archaic allels. And Ancestral allles (if possible)
    dct_AA=read_aa_file(f_aa)
    with open(file,'r') as f :
        lns=f.readlines()
        
    nested_dict = {}
        
    for i in range(1,len(lns)):
        m=lns[i].replace(' \n','').split('\t') 
        out = [int(m[3][k])  for k in range(len(m[3])) if m[3][k]!='.']      
        arch = [int(m[4][k])  for k in range(len(m[4])) if m[4][k] != '.']    
        o=m[5].split(' ')
        o = [int(o[j]) for j in range(len(o))]
   
        nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'Outgroup': out, 'Archaic': arch, 'Obs': o  }

    for j in dct_AA.keys():
        AA=dct_AA[j]
        if nested_dict.get(j) is not None:
            st=nested_dict[j]['REF']+nested_dict[j]['ALT'].replace(',','')
        
        
            for i in range(len(st)):
                if st[i]==AA:
                    k=i        
        
            nested_dict[j]['AA']=k
    for i, j in nested_dict.items():
        if j.get('AA') is None:
            nested_dict[i]['AA']='.'

    return nested_dict

def main_read2(file2): # return the dictionary with REF,ALT, Outgroup and Archaic allels. And Ancestral allles (if possible)
    
    with open(file2,'r') as f :
        lns=f.readlines()
    
    nested_dict = {}
        
    for i in range(1,len(lns)):
        
        m=lns[i].replace(' \n','').split('\t') 

        out=m[4].replace(',','')               
        if out!='':
            out = [int(out[k])  for k in range(len(out))]  
        else:
            out = []
        
        arch= m[5].replace('.','')
        if len(arch)!=0:
            arch=[int(k) for k in arch.split(',')]
        else:
            arch=[]

        o=m[6].replace('\n','').split(' ')
        o=[int(k) for k in o ]
        


        if m[3]=='.': 
            s=''
            nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'Outgroup': out, 'Archaic': arch, 'Obs': o  }
            
        else:
            s=m[3]   
            nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'Outgroup': out, 'Archaic': arch, 'Obs': o, 'AA': s  }


    return nested_dict

      
def make_obs_ref(nested_dict, domain, ind, L,  ref):


        
    start, end = domain[0][0], domain[-1][1]
    if float(end-start) % L != 0:
        T = int((end-start)/ L) + 1
    obs_ref = np.zeros(T, int)  
    
    for k in nested_dict.keys():
        j=int((k-start)/L) #window-positions
        
        if nested_dict[k].get('AA') is None:
        
            if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0:
                obs_ref[j]+=1
        else:
            if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0 and nested_dict[k]['Obs'][ind]!=nested_dict[k]['AA']:
                obs_ref[j]+=1
    return obs_ref
    
    
    
def read_par_HMM(file):
    with open(file,'r') as f:

        GEN_time = float(f.readline())
        MU = float(f.readline())
        RR = float(f.readline())
        L = int(f.readline())



        Lambda_0=np.zeros(5)
        Lambda_0[1] = float(f.readline())/GEN_time*MU*L
        Lambda_0[2] = float(f.readline())/GEN_time*MU*L
        Lambda_0[0] = float(f.readline())/GEN_time*MU*L
        Lambda_0[4] = float(f.readline())/GEN_time*MU* L
        Lambda_0[3] = float(f.readline())

    return GEN_time, MU, RR, L, Lambda_0
   
   
def ne_set_interval(set_intervals):
    
    ne_set = []
    for j in range(len(set_intervals)):
        
        if j==0:
            ne_set.append([0,set_intervals[j][0]-1])
        if j== len(set_intervals)-1:
            ne_set.append([set_intervals[j][1]+1, seq_length-1])
        if j!=0 and j!= len(set_intervals)-1:
            ne_set.append([set_intervals[j-1][1]+1, set_intervals[j][0]-1])
    return ne_set

def len_tracts(set_intervals):
    if len(set_intervals)==0:
        return 0
    else:
        s=0
        for j in range(len(set_intervals)):
            s+= set_intervals[j][1]-set_intervals[j][0]+1
        return s
    


def read_out(file):
    with open(file, 'r') as f:
        l=f.readlines()

    l=[l[i].split('\t') for i in range(len(l))]
    l=[l[i][2] for i in range(len(l))]

    nd_HMM_tracts=[]
    for j in range(len(l)):
        m1=[]
        l2=l[j].replace('[[','').replace(']]','').replace('\n','').split('], [')
        for j1 in range(len(l2)):
            m2=l2[j1].split(', ')
            m2=[int(m2[j2]) for j2 in range(len(m2))]
            m1.append(m2)
        nd_HMM_tracts.append(m1)
    return nd_HMM_tracts
    
    
# return European tracts with input=Neanderthal tracts
def tracts_eu(tr_nd, seq_length):
    result = []

    if tr_nd[0][0] > 0:
        result.append([0,tr_nd[0][0]-1])
        
    for i in range(len(tr_nd)-1):
        result.append([tr_nd[i][1]+1, tr_nd[i+1][0]-1])
        
    if tr_nd[-1][1]!=seq_length-1:
        result.append([tr_nd[-1][1]+1,seq_length-1])
      
    return result    

def noND_txt(ts_tree, n_ref_pop, n_eu):
    f_names = ['./Skov'+str(i)+'.txt' for i in range(n_eu)]


    for i in range(n_eu):
        with open(f_names[i], "w") as file:
            file.write("chrom  pos     ancestral_base  genotype\n")
            for v in ts_tree.variants():       
                flag = False
                inx=[k+n_eu for k in range(n_ref_pop)]
                for j in inx:
                    if abs(v.genotypes[j]-v.genotypes[i])==0:
                        flag = True
                        break                

                if flag == False and v.genotypes[i] == 1:              
                    file.write('chr1   '+str(int(v.position))+'    0               1\n')


def read_noND(cutoff, n_eu):

    f2_names = ['./Skov/Skov.result.'+str(i)+'.hap1.txt' for i in range(n_eu)]
    tracts_skov_neand = []
    for j in range(n_eu):
        
        with open(f2_names[j],'r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            lines[i] = lines[i].split('\t')
        lines2=[]    
        for i in range(len(lines)):
           
            if 'Archaic' in lines[i] and float(lines[i][5])>cutoff:
                lines2.append(lines[i])
    
        tracts_skov1 = []
        for i in range(len(lines2)):
            tracts_skov1.append([int(lines2[i][1]), int(lines2[i][2])])
        
        tracts_skov_neand.append(tracts_skov1)
    return tracts_skov_neand
