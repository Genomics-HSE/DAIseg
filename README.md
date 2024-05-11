
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 

__Input__: .vcf.gz{.tbi} file where merging all Neanderthal, Outgroupand  ingroup observable samples together and three .txt files with ids of each group.

__Output__: .txt file where each line corresponds to the array of tracts with modern and archaic ancestries.


# Pipeline briefly
0. (optionally) Run __panel.preparation.Linux.sh obs.samples.txt outgroup.txt__(or __panel.preparation.MacOS.sh__) with samples' name files to merge 1000GP, neanderthal samples and obtain .vcf.gz file.
1. Using .vcf.gz{.tbi} and files with samples's names to run __./make.obs.sh outgroup.txt archaic.txt obs.samples.txt__ to make observation files.
2. Run __dai.seg.py__ to obtain archaic tracts of samples from  __observations.txt__  with the posssibility of using EM algorithm.



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
start end    #position of the first SNP in .vcf file
t_arch^c    #Coalescent time of AMH and Neanderthals
t_split^c    #Coalescent time out of Africa 
t_intr^c    #coalescent time of archaic segments in modern genome with neanderthal samples
t_intr #introgression time 
0.025    #admixture proportion of archaic introgression
```

By default, the  time values are  550.000, 70.000, 55.000, 55.000 are used to make  initiall guess for the EM algorithm on Step 2. These values are good to find archqic segments but using EM algorithm allows to find short segments.


*  __all.chr22.vcf.gz{.tbi}__ files containing all reference genomes (Outgroup and Archaic) and observable samples with snps only (excluding indels, deletions etc.). The main reason of it is to avoid inconsistencies.
  
* __output.tracts.txt__ is a file 
```note
[[[t_1,t_2], [t_3,t_4], [t_5,t_6]], [[t_2+1, t_3-1], [t_4+1, t_5-1]]]
[[[t'_1,t'_2], [t'_3,t'_4], [t'_5,t'_6]], [[t'_2+1, t'_3-1], [t'_4+1, t'_5-1]]]
...
...
```
where each two lines correspond to the one diploid sample from obs.samples.txt.





# Step 0. Merging 1000GP  and Archaic genomes (~5 min)
1. Download data:
- 1000GP panel:
>http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 

- Archaic samples:
>http://cdna.eva.mpg.de/neandertal/Vindija/VCF/

>http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/ (split by chromosomes!)

2. Make .txt files with the corresponding samples' names:
- __outgroup.txt__
- __obs.samples.txt__
- __archaic.txt__

3. Edit __panel.preparation.sh__:
- `$NAME1000` = Add full path to 1000GP data
- `$n1, $n2, $n3` = Add full path to Archaic data
- `$CHR` = chromosome of interest
  
4. Run command:

Linux: 
```bash
./panel.preparation.Linux.sh obs.samples.txt outgroup.txt
```

Mac:
```bash
./panel.preparation.MacOS.sh obs.samples.txt outgroup.txt
```
 
5. Result:
   
- Merged vcf: __all.chr${CHR}.vcf.gz__
- Index file for vcf: __all.chr${CHR}.vcf.gz.tbi__

# Step 1. Make observations (~ 3 min for two samples)

1. Run command:

```bash
./make.obs.sh outgroup.txt archaic.txt obs.samples.txt
```

For that you need files __all.chr${CHR}.vcf.gz{.tbi}__,  __outgroup.txt__, __observations.txt__, __archaic.txt__ 

2. Result (see the File's summary paragraph):

- Files with observations: __obs.neand.txt__, __obs.outgroup.txt__
- File with default parameters: __par.file.txt__
- File with start-end positions: **positions.chr{$CHR}.txt**

# Step 2. Run DAI.seg 
## Without EM algorithm

1. Run command:
```python
python3 dai.seg.py  --EM no --HMM_par par.file.txt  --o output.tracts.txt
```

2. Result (see the File's summary paragraph):

- __output.tracts.txt__

## With EM algorithm

- par.file.txt obtained on the Step 1 could be used as the initial guess for EM algorithm.

1. Run command:

To obtain estimations only for coalescent times: `--EM_est coal`

```python
python3 dai.seg.py --EM yes --EM_est coal --HMM_par par.file.txt --o output.tracts.txt
```
To make estimations of coalescent times, split times and admixture proportion: `--EM_est all`

```python
python3 dai.seg.py --EM yes --EM_est all --HMM_par par.file.txt --o output.tracts.txt
```

2. Result (see the File's summary paragraph):

- __output.tracts.txt__
