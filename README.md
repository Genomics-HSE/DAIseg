
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes.  



The description of the used files is [here][1]









There are two options without EM-algorithm and with EM algorithm. 


# Run DAI.seg without EM algorithm
To run DAIseg you need
1) 'obs.txt' with  the following structure
```bash
#POSITIONS	#REF	#ALT	ANCESTRAL	#OUTGROUP	#ARCHAIC	#OBSERVATIONS
48	G	A	.	0	0	0 0 0 0 
67	A	T	0	0	0,1	0 0 0 0 
143	T	A	0	0,1	.	0 0 0 0
...
```
which consists of the rows corresponding to biallelic SNPs with columns of SNS's position, REF and ALT alleles, Ancestral allele(if possible), Outgroup(African) and Neanderthal variants(if it is possible) and 2*number_of_samples haplotypes. 

2) samples.txt = list of samples
```bash
name1
name2
...
```

3) .bed file
```bash
CHR	po1	pos2
CHR	po3	pos4
...
```
5) 


```bash
python3 dai.seg.py --obs_samples path.to/obserables.list --bed path.to/file.bed   --HMM_par par.file.txt --EM no --prepared_file ./preparations/hg19.all.chr22.txt --o out.chr22 --arch_cover ./preparations/hg19.arch.covering.chr22.txt --decoding posterior/viterbi --cut_off 0.9
```


# Run DAI.seg using EM algorithm

```bash
python3 dai.seg.py --obs_samples path.to/obserables.list --bed path.to/file.bed   --HMM_par par.file.txt --EM yes --EM_steps 20  --EM_samples 5 --prepared_file ./preparations/hg19.all.chr22.txt --o out.EM.chr22 --arch_cover ./preparations/hg19.arch.covering.chr22.txt  --decoding posterior/viterbi --cut_off 0.9
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.

# Preparations of files

Download file with ancestrall alleles from [here][5] or go to the directory ./preparations/ to run script and get ancestral alleles from 1000GP hg19.

To read more details for files preparation see [readme][2]. To avoid details use script 
```bash
 ./full.preparation.Linux(MacOS).sh hg19 22 path.to/file.bed n1 n2 n3 1000GP path.to/obserables.list path.to/outgroup.list path.to/archaic.list name.out_vcf name.out_txt path.to/ancestral.alleles
```
where 1000GP, n1, n2, n3 are [vcf with 1000 Genome project][3], [neanderthal vcfs][4] and [chagyrskaya][5] (needed to be splited). 

Write full path to the list of samples path.to/outgroup.list,  path.to/obserables.list. path.to/archaic.list

name.out_vcf name.out_txt are the names of the resulting files(be saved in the hg19 directory).




[1]: https://github.com/Genomics-HSE/DAIseg/blob/main/File.types.md
[2]: https://github.com/Genomics-HSE/DAIseg/blob/main/hg19/README.md

[5]: https://drive.google.com/drive/folders/1_zE9eaV3psFPRdFatkq-R1yGluvjgiX6?usp=sharing





[3]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[4]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[5]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/

