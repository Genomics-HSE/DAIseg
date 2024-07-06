
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 


# Files's summary
*  __outgroup.txt__(Africa), __archaic.txt__(Neanderthals)  and __obs.samples.txt__(European), are .txt files which consist of the samples' ids of reference Africans, Neanderthals and observable Europeans written in a column
   ```note
   NA18484
   NA18489
   GM19129
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


*  __all.chr22.vcf.gz{.tbi}__ files containing all reference genomes (Outgroup and Archaic) and observable samples with snps only (excluding indels, deletions etc.). The main reason of it is to avoid inconsistencies.
  
 * __output.txt__ is a  file 
    ```note
    HG01510    0    [[t_1,t_2], [t_3,t_4], [t_5,t_6]]
    HG01510    1    [[t'_1,t'_2], [t'_3,t'_4], [t'_5,t'_6]]
    ...
    ...
    ```
    where each two lines correspond to the one diploid sample from obs.samples.txt.

*  __pos.chr22.txt__ first-last positions of desired region
    ```note
    start_chr end_chr
    ```


*  __POS.AA.chr22.txt__  file with information about ancestral allels ("-1"=="no information").
     ```note
     position1 -1
     position2 A
     position3 C
     ...
     ```
   The link on the [ancestral alles files based on hg19][4] 


* __gaps.by.pos.chr.22.txt__ is file with list of gaps
  ```note
  [[a,b], [c,d]. [e,f]]
  ```

## Merging 1000GP  and Archaic genomes
Download [1000GP panel][1] and  archaic samples  [Link1][2] and [Link2][3]. Make .txt files with samples' names  obs.samples.txt, outgroup.txt, archaic.txt

Add full path to files  of 1000GP,  Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables NAME1000 and n1, n2, n3 in  panel.preparation.*.sh and run 

```bash
./panel.preparation.Linux.sh 22 obs.samples.txt outgroup.txt
```
 
The resulting vcf.gz file is all.chr22.vcf.gz{.tbi}


## Make observations 

You need  vcf file, lists of samples obs.samples.txt, outgroup.txt, archaic.txt and file with [ancestral alleles positions][4]
 to run  

```bash
./make.obs.sh 22 all.chr22.vcf.gz obs.samples.txt Outgroup.txt archaic.txt  POS.AA.chr22.txt
```

and to make observation files obs.neand.chr22.txt, obs.outgroup.chr22.txt




## Run DAI.seg without EM algorithm
```bash
python3 dai.seg.py --obs_samples obs.samples.txt --location pos.chr22.txt --gaps gaps.by.pos.chr22.txt --HMM_par par.file.txt --EM no --obs_af obs.outgroup.chr22.txt --obs_archaic obs.archaic.chr22.txt --o out.chr22.txt
```

where the examples par.file.txt, pos.chr22.txt, POS.AA.chr22.txt could be found in the main directory

## Run DAI.seg using EM algorithm

```bash
python3 dai.seg.py --obs_samples obs.samples.txt --location pos.chr22.txt --gaps gaps.by.pos.chr.22.txt --HMM_par par.file.txt --EM yes --obs_af obs.outgroup.chr22.txt --obs_archaic obs.archaic.chr22.txt --o out.chr22.txt
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.


[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
[4]: https://drive.google.com/file/d/1Vw-QEG9uu1trkbGHpDVXhMlbGt-RQhbN/view?usp=sharing
