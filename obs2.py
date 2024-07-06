import numpy as np
import sys
import useful as usfl

CHR=sys.argv[1]
f=sys.argv[2]
f_aa = sys.argv[3]
f_bed = sys.argv[4]


domain = usfl.read_bed(f_bed)
dct_all = usfl.main_read(f, f_aa)



L=1000

print(f)
print(f[5:])
with open(f[5:],'w') as f1:
    f1.write('#POSITIONS\t#REF\t#ALT\tANCESTRAL\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n')
    for i  in dct_all.keys():
        j=dct_all[i]        
        s1=str(j['Outgroup']).replace('[','').replace(']','').replace(' ','')
        s2=str(j['Archaic']).replace('[','').replace(']', '').replace(' ','')
        if s2=='':
            s2+='.'
        s3=str(j['Obs']).replace('[','').replace(']','').replace(',','')
        f1.write(str(i)+'\t'+str(j['REF'])+'\t'+str(j['ALT'])+'\t'+str(j['AA'])+'\t'+s1+'\t'+s2+'\t'+s3+'\n')


n_eu = len(dct_all[list(dct_all.keys())[0]]['Obs'])

SEQ=[]
N_ST=[]


for ind in range(n_eu):
    sq=np.vstack([usfl.make_obs_ref(dct_all, domain, ind, L,  'Outgroup'), usfl.make_obs_ref(dct_all, domain, ind, L,  'Archaic')])
    sq=sq.transpose()
    n_st = sq.max()+1
    SEQ.append(sq)
    N_ST.append(n_st)
SEQ=np.array(SEQ)








