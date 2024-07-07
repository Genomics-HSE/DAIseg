# Files's summary
*  __outgroup.txt__(Africa), __archaic.txt__(Neanderthals)  and __obs.samples.txt__(European), are .txt files which consist of the samples' ids of reference Africans, Neanderthals and observable Europeans written in a column
   ```note
   NA18484
   NA18489
   GM19129
   ```





*  __all.chr22.vcf.gz{.tbi}__ files containing all reference genomes (Outgroup and Archaic) and observable samples with snps only (excluding indels, deletions etc.). The main reason of it is to avoid inconsistencies.
  
 * __output.txt__ is a  file 
    ```note
    HG01510    0    [[t_1,t_2], [t_3,t_4], [t_5,t_6]]
    HG01510    1    [[t'_1,t'_2], [t'_3,t'_4], [t'_5,t'_6]]
    ...
    ...
    ```
    where each two lines correspond to the one diploid sample from obs.samples.txt.


*  __hg19.AA.chr22.txt__  file with information about ancestral allels.
     ```note
     position_0 A
     position_1 C
     ...
     ```
   The link on the [ancestral alles files based on hg19][4] 

* __region.bed__ is file with desired regions
  ```note
  22	16050000	16697850
  22	16847850	20509431
  22	20609431	50364777
  22	50414777	51244566
  ```
  
*  __par.file.txt__
   ```note
   29 # years per generation
   1.25e-08    #mutation rate μ
   1e-08    #recombination rate
   1000    #window size
   t_arch^c    #Coalescent time of AMH and Neanderthals
   t_split^c    #Coalescent time out of Africa 
   t_intr^c    #coalescent time of archaic segments in modern genome with neanderthal samples
   t_intr #introgression time 
   0.025    #admixture proportion of archaic introgression
     ```


     By default, the  time values are  550.000, 70.000, 55.000, 55.000 are used to make  initiall guess for the EM algorithm on Step 2. These values are good to find archqic segments but using EM algorithm allows to find short segments.

* __allels.ref.and.obs.chr22.txt__ is a file with all needed informations such as REF/ALT alleles, Ancestral Allele, Outgroup and Archaic Alleles, and Observations (haplotypes are written in columns )
     ```note
     #POSITIONS	#REF	#ALT	ANCESTRAL	#OUTGROUP	#ARCHAIC	#OBSERVATIONS
     16050075	A	G	.	0	.	0 0
     16050115	G	A	.	1,0	.	0 0
     ...
     21954214	C	T	0	0	0	0 0
     ```

* __arch.covering.chr22.txt__ file with window covering by neanderthals (-0.001 means that there is no information in this window, 1.0 means the there is some information in each positions of the window from at least one archaic sample )
  ```note
  -0.001
  -0.001
  ...
  0.022
  0.001
  ...
  0.999
  ```
Where the first window starts with the fist position in the region.bed file. The length of window is 1000.
