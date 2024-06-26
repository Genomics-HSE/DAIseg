import numpy as np
import sys
import useful as usfl

CHR=sys.argv[1]

f_yri = sys.argv[2]
f_neand = sys.argv[3]
f_obs = sys.argv[4]
f_aa = sys.argv[5]



 
def make_obs(lines, lines_ref, L, ind,dc):

    start, end = lines_ref[0,0], lines_ref[-1,0]
    if float(end-start) % L != 0:
        T = int((end-start)/ L) + 1
    obs_ref = np.zeros(T, int)
    
    obs = lines[:,ind]
    
    for i in range(len(lines)):
        j=int((lines_ref[i][0]-start)/L)

        if obs[i]==0 and dc[lines_ref[i][0]]!=obs[i]:
            if lines_ref[i][1]==-1:
                obs_ref[j]+=1
                
        if obs[i]==1 and dc[lines_ref[i][0]]!=obs[i]:
            if lines_ref[i][2]==-1:
                obs_ref[j]+=1
                
    return obs_ref   
    
    
with open(f_aa,'r') as f:
    l=f.readlines()
 
dct={}
for  i in l:
    
    m=i.replace(',','').replace('[','').replace(']','').replace('\n','').split(' ')

    dct[int(m[0])]=int(m[1])


with open(f_obs,'r') as f:
    
    
    lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].strip(' \n').split(' ')

        for j in range(len(lines[i])):
            lines[i][j]=int(lines[i][j])
lines = np.array(lines)  

    
    
   
with open(f_neand,'r') as f:
    lines_neand = f.readlines()
with open(f_yri,'r') as f:
    lines_yri = f.readlines()
    
for i in range(len(lines_yri)):
    lines_yri[i] = lines_yri[i].strip('\n').split('\t')
    for j in range(len(lines_yri[i])):
        lines_yri[i][j]=int(lines_yri[i][j])
        
for i in range(len(lines_neand)):
    lines_neand[i] = lines_neand[i].strip('\n').split('\t')
    for j in range(len(lines_neand[i])):
        lines_neand[i][j]=int(lines_neand[i][j])    
        
lines_yri = np.array(lines_yri)
lines_neand = np.array(lines_neand)





start, end = lines_yri[0,0], lines_yri[-1,0]

n_eu = len(lines[0])
#n_eu=3
SEQ=[]
N_ST=[]
L=1000




        

for ind in range(n_eu):
    sq=np.vstack([make_obs(lines, lines_yri, L, ind, dct ),make_obs(lines, lines_neand, L, ind, dct)])
    sq=sq.transpose()
    n_st = sq.max()+1
    SEQ.append(sq)
    N_ST.append(n_st)
SEQ=np.array(SEQ)


with open('pos.chr'+str(CHR)+'.txt','w') as f:
    f.write(str(lines_yri[0,0])+' ' +str(lines_yri[-1,0])+'\n')


with open('obs.outgroup.chr'+str(CHR)+'.txt', "w") as file1,  open('obs.archaic.chr'+str(CHR)+'.txt', "w") as file2:
    for j in range(len(SEQ[0])):
        s1, s2 = '',''
        for i in range(n_eu):
            s1 += str(SEQ[:,j,0][i])+' '
            s2 += str(SEQ[:,j,1][i])+' '    
        file1.write(s1[:-1]+'\n')
        file2.write(s2[:-1]+'\n')
        
        
