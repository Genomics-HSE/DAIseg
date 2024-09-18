import sys


CHR=int(sys.argv[1])



for cr in range(CHR,CHR-1, -1):
    with open('./Ancestral.Alleles/POS.REF.ALT.INFO.chr'+str(cr)+'.txt','r') as f:
        l=f.readlines()
    m_AA=[]
    m_AA_bi=[]
    kk=0
    flag=False
    for i in l:
        kk+=1
        m=i.split(' ')
        
        k=m[3].find('AA=')
        if k!=-1:
            s=m[3][k+3] 
            flag=True

        if m[2] in ['A', 'C', 'G', 'T'] and m[1] in ['A', 'C', 'G', 'T'] and flag:


	    
            
            
            
            if s in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']:
                s=s.upper()
                if s=='a':
                    m_AA_bi.append([int(m[0]), 'A'])
                if s=='c':
                    m_AA_bi.append([int(m[0]), 'C'])
                if s=='g':
                    m_AA_bi.append([int(m[0]), 'G'])
                if s=='t':
                    m_AA_bi.append([int(m[0]), 'T'])    
                else:                                  
                    m_AA_bi.append([int(m[0]), s])

        flag=False

            
        


    with open('./Ancestral.Alleles/grch38.AA.chr'+str(cr)+'.txt','w') as f:
        for i in m_AA_bi:
            f.write('chr'+str(CHR)+'\t'+str(i[0])+'\t'+str(i[1])+'\n')

