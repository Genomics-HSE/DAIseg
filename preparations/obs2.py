import numpy as np
import sys
import useful as usfl

CHR=sys.argv[1]
f=sys.argv[2]
f_aa = sys.argv[3]
f_bed = sys.argv[4]
f_out = sys.argv[5]

domain = usfl.read_bed(f_bed)
dct_all = usfl.main_read(f, f_aa)

print(len(dct_all.keys()))

L=1000
k=0
with open(f_out,'w') as f1:
    f1.write('#POSITIONS\t#REF\t#ALT\tANCESTRAL\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n')
    for i  in dct_all.keys():
        j=dct_all[i]        
        s1=str(j['Outgroup']).replace('[','').replace(']','').replace(' ','')
        s2=str(j['Archaic']).replace('[','').replace(']', '').replace(' ','')
        if s2=='':
            s2+='.'
        s3=str(j['Obs']).replace('[','').replace(']','').replace(',','')
        
        flag=[]
        if j['Archaic']==[]:
            for o in set(j['Obs']):
                if o in set(j['Outgroup']):
                    flag.append(True)
            if flag==[True for o in set(j['Obs'])]:
                k+=1
   
            else:
                f1.write(str(i)+'\t'+str(j['REF'])+'\t'+str(j['ALT'])+'\t'+str(j['AA'])+'\t'+s1+'\t'+s2+'\t'+s3+'\n')
        else:
            for o in set(j['Obs']):
                if (o in set(j['Outgroup'])) and (o in set(j['Archaic'])):
                    flag.append(True)
            if flag==[True for o in set(j['Obs'])]:
                k+=1
   
            else:
                f1.write(str(i)+'\t'+str(j['REF'])+'\t'+str(j['ALT'])+'\t'+str(j['AA'])+'\t'+s1+'\t'+s2+'\t'+s3+'\n')        
                











