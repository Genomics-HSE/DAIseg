
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 


# Run DAI.seg without EM algorithm
```bash
python3 dai.seg.2.py --obs_samples obs.samples.txt --bed test.bed   --HMM_par par.file.txt --EM no --prepared_file allels.ref.and.obs.chr22.txt --o out.chr22.txt --arch_cover arch.covering.chr22.txt
```

where the examples par.file.txt could be found in the main directory

# Run DAI.seg using EM algorithm

```bash
python3 dai.seg.2.py --obs_samples obs.samples.txt --bed test.bed   --HMM_par par.file.txt --EM yes --EM_steps 20  --prepared_file allels.ref.and.obs.chr22.txt --o out.EM.txt --arch_cover arch.covering.chr22.txt
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.

if you have the files obs.samples.txt, region.bed, par.file.txt, arch.covering.chr22.txt, allels.ref.and.obs.chr22.txt (the detailed description see in __File's summary__ section) 



     





## Merging 1000GP  and Archaic genomes


Download [1000GP panel][1] and  archaic samples  [Link1][2] and [Link2][3]. Make .txt files with samples' names  obs.samples.txt, outgroup.txt, archaic.txt

Add full path to files  of 1000GP,  Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables NAME1000 and n1, n2, n3 in  panel.preparation.*.sh and run 

```bash
./new.panel.preparation.Linux.sh 22 Outgroup.txt obs.samples.txt test.bed all.chr22.vcf.gz
```
 
The resulting vcf.gz file is all.chr22.vcf.gz{.tbi}

## Ancestral alleles

If you working with hg19 the list of acestral allels could be extract from vcf [1000GP panel][1]. Run
```bash
./ancestral.alleles.sh 22 way.to.1000GP.chr22
```

to make file hg19.AA.chr22.txt in created directory Ancestral.Allels.




## Archaic covering

Add full path to files  of   Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables n1, n2, n3 in  archaic.covering.sh and run 
```bash
./archaic.covering.sh 22 test.bed
```
to obtain the window-covering by archaic samples.





## Make observations 

You need  vcf file, lists of samples obs.samples.txt, outgroup.txt, archaic.txt and file with [ancestral alleles positions][4]
 to run  

```bash
./new.make.obs.sh 22 all.chr22.vcf.gz obs.samples.txt Outgroup.txt archaic.txt  ./Ancestral.Alleles/hg19.AA.chr22.txt test.bed
```

and to make observation files obs.neand.chr22.txt, obs.outgroup.chr22.txt



# Files's summary
*  __outgroup.txt__(Africa), __archaic.txt__(Neanderthals)  and __obs.samples.txt__(European), are .txt files which consist of the samples' ids of reference Africans, Neanderthals and observable Europeans written in a column
   ```note
   NA18484
   NA18489
   GM19129
   ```



   ```

     By default, the  time values are  550.000, 70.000, 55.000, 55.000 are used to make  initiall guess for the EM algorithm on Step 2. These values are good to find archqic segments but using EM algorithm allows to find short segments.


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




[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
[4]: https://drive.google.com/file/d/1Vw-QEG9uu1trkbGHpDVXhMlbGt-RQhbN/view?usp=sharing
