import numpy as np
import useful as usfl

import sys

f_bed=sys.argv[1]
f_pos=sys.argv[2]
CHR=sys.argv[3]

genome=sys.argv[4]

with open(f_bed,'r') as f:
    d=f.readlines()
d=[d[i].replace('\n','').split('\t')[1:3] for i in range(len(d))]    

for i in  range(len(d)):
    d[i][0]=int(d[i][0])
    d[i][1] = int(d[i][1])
domain =d


with open(f_pos,'r') as f:
    lines=f.readlines()

lines=[int(lines[i].replace('\n', '')) for i in range(len(lines))]
L=1000
n_windows=(domain[-1][1]-domain[0][0])//L + 1

windows_cover=-np.ones(n_windows)
for i in lines:
    if usfl.point_in_set(i, domain) is True:
        windows_cover[(i-domain[0][0])//L]+=1
windows_cover=windows_cover/L

with open(str(genome)+'.arch.covering.chr'+str(CHR)+'.txt','w') as f:
    for j in windows_cover:
        f.write(str(j)+'\n')

