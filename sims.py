# --- Demographic model ---


import msprime 
import numpy as np
import useful

def history_archaic(gen_time, len_seq, rr,mu, n_e, t,  n, rand_sd, n_neand, t_neand_samples, n_eu, n_eu_growth, t_eu_growth, n_eu_bottleneck, gr_rt, p_admix,ploid ):
   
    n_ANC, n_ND, n_AMH, n_OOA,n_AF, n_EU = n_e
    t_NEAND_migration, t_NEAND_AMH, t_OOF_AF = t

    
    demography = msprime.Demography()
    
    
    demography.add_population(name="AF", initial_size=n_AF)
    demography.add_population(name="EU", initial_size=n_EU)
    demography.add_population(name="AMH", initial_size=n_AMH)
    demography.add_population(name="NEAND", initial_size=n_ND)
    demography.add_population(name="ANCES", initial_size=n_ANC)  #common population for Neanderthal and AMH
    demography.add_population(name="OOA", initial_size = n_OOA)
    
    
    demography.add_population_parameters_change(time=0, initial_size=n_EU, population=1, growth_rate=gr_rt)
    demography.add_population_parameters_change(time=t_eu_growth, initial_size=n_eu_growth, population=1, growth_rate=0)   
    
    demography.add_admixture(time=t_NEAND_migration, derived="EU", ancestral=["OOA", "NEAND"], 
                             proportions=[1-p_admix, p_admix])


    demography.add_population_split(time = t_OOF_AF, derived=["AF", "OOA"], ancestral="AMH")
    demography.add_population_split(time = t_NEAND_AMH, derived=["AMH", "NEAND"], ancestral="ANCES")

#    print(demography.debug())

    ts = msprime.sim_ancestry(
        samples=        [       
                msprime.SampleSet(n_eu, ploidy=ploid, population='EU'), 
                msprime.SampleSet(n, ploidy=ploid, population='AF'),
                msprime.SampleSet(n_neand, ploidy=ploid, population='NEAND', time = t_neand_samples)          
               
        ],
    
        ploidy=ploid,    
        sequence_length=len_seq,
        recombination_rate=rr, 
        demography=demography,
        record_migrations=True   
    )
   
    
    ts = msprime.sim_mutations(ts, rate=mu)    
    return ts


#несколько вспомогательных функций
def connected(m):
    for i in range(len(m)-1):
        if m[i][1] == m[i+1][0]:
            return True
    return False
        
def remove_one(m):
    mas = m
    while connected(mas) == True:
        for i in range(len(mas)-1):
            if mas[i][1] == mas[i+1][0]:
                mas[i][1] = mas[i+1][1]
                mas.pop(i+1)
                break
    return mas


#Вход: ts, название популяции, индивид(которого мы препарируем), время предка
def get_migrating_tracts_ind(ts, pop_name, ind, T_anc):
    
    pop=-1
    for i in ts.populations():
        if i.metadata['name']==pop_name:
            pop=i.id
    
    mig = ts.tables.migrations
    migration_int = []

    for tree in ts.trees():  #перебираем все деревья. Как известно, каждому дереву отвечает участок днк  
        anc_node = ind #chose observable node
        while tree.time( tree.parent(anc_node) ) <= T_anc : #идем в прошлое до вершины anc_node по предкам нашего мексиканца, пока не наткнемся на миграцию 
            anc_node = tree.parent(anc_node)
        migs = np.where(mig.node == anc_node)[0] #выбирем все строки, соответствующие заданному узлу

        #идем по таблице миграций с anc_node и проверяем, чтобы миграции попадали в тот самый участок днк
        for i in migs:

            stroka = mig[i]
            if stroka.time == T_anc and stroka.dest == pop and tree.interval.left >= stroka.left and tree.interval.right <= stroka.right:
                migration_int.append([tree.interval.left, tree.interval.right])

    migration_int2 = []
    for i in range(len(migration_int)):
        if migration_int[i][0] != migration_int[i][1]:
            migration_int2.append(migration_int[i])
    migration_int = migration_int2
    
    mi = remove_one(migration_int)
    mi.sort()  

    return mi 


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



#from Skov 
def print_neand_dosages(ts):
    
    seq_len = ts.get_sequence_length()

    ME_ids = ts.get_samples(1)
    de_seg = {i: [] for i in ME_ids}
    ar_seg = {i: [] for i in ME_ids}
    for mr in ts.migrations():
        if mr.source == 1 and mr.dest == 3:
            for tree in ts.trees(leaf_lists=True):
                if mr.left > tree.get_interval()[0]:
                    continue
                if mr.right <= tree.get_interval()[0]:
                    break
                for l in tree.leaves(mr.node):
                    if l in ME_ids:
                        # print(l)
                        de_seg[l].append(tree.get_interval())

    def combine_segs(segs, get_segs = False):
        merged = np.empty([0, 2])
        if len(segs) == 0:
            if get_segs:
                return([])
            else:
                return(0)
        sorted_segs = segs[np.argsort(segs[:, 0]), :]
        for higher in sorted_segs:
            if len(merged) == 0:
                merged = np.vstack([merged, higher])            
            else:
                lower = merged[-1, :]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1, :] = (lower[0], upper_bound) 
                else:
                    merged = np.vstack([merged, higher])
        if get_segs:
            return(merged)
        else:
            return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)

    true_de_prop = [combine_segs(np.array(de_seg[i])) for i in sorted(de_seg.keys())]
    true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]
    print("Neand ancestry: ", true_de_prop)
      

def create_obs_txt(file, ploidy, n_diplo, ts, n_ref_pop, n_neanderthal, n):
    N_neanderthal, n_eu, n, N_ref_pop =ploidy*n_neanderthal, ploidy * n_diplo, ploidy*n, ploidy* n_ref_pop
    
    with open(file, 'w') as f:
        f.write('#POSITIONS\t#REF\t#ALT\tANCESTRAL\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n')
        for v in ts.variants():
            outgroup= str(list(set(v.genotypes[n_eu :( n_eu+N_ref_pop)]))).replace('[','').replace(']','').replace(' ','')
            archaic= str(list(set(v.genotypes[n_eu+n :( n_eu+n+N_neanderthal)]))).replace('[','').replace(']','').replace(' ','')
    
            obs=''
            for i in v.genotypes[0 :n_eu]:
                obs+=str(i)+' '
    
            f.write(str(int(v.site.position))+'\t'+str(v.alleles[0])+'\t'+
                    str(v.alleles[1]) + '\t'+ str(0)+'\t' + outgroup+'\t'+archaic+'\t'+str(obs)+'\n')    



def create_bed_smpls_arch_cov(sample_file, bed_file, arch_cover_file,len_sequence, cover, CHR, n_eu_diplo, ploidy):
    N_eu=n_eu_diplo*ploidy
    L=1000
    
    #create bed file
    with open(bed_file,'w') as f:
        f.write('1\t0\t'+str(int(len_sequence)-1)+'\n')
    
    domain=useful.read_bed(bed_file)
    
    n_windows=(domain[-1][1]-domain[0][0])//L + 1
    windows_cover=np.ones(n_windows)*cover
    
    #create archaic covering file. 
    
    with open(arch_cover_file,'w') as f:
        for j in windows_cover:
            f.write(str(j)+'\n')
    
    
    #create file with sample's names
    with open(sample_file,'w') as f:
        for i in range(N_eu):
            f.write('eu'+str(i)+'\n')

def real_nd_tracts(ts, n_eu_diplo, ploidy, T):
    
    ND_true_tracts = []
    for idx in range(0, ploidy*n_eu_diplo): 
        if (idx % 20) ==0 and n_eu_diplo>20:
            print('Done', idx)
        ND_true_tracts.append( get_migrating_tracts_ind(ts, 'NEAND', idx, T[0]))      
       
    s=0
    for i in range(ploidy*n_eu_diplo):
        for j in ND_true_tracts[i]:
            s+=j[1]-j[0]    

    return ND_true_tracts
    print('средняя доля неандертальца',s/(n_eu * len_sequence))
