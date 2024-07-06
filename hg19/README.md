
# Preparation connected with hg19

Choose this directory.

## Archaic covering 
The goal is to create __arch.covering.chr22.txt__ file with the window-covering by archaic samples. 
Add full path to files  of   Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables n1, n2, n3 and run 
```bash
./archaic.covering.sh 22 path.to/file.bed n1 n2 n3
```


## Ancestral alleles
The goal is to create __hg19.AA.chr22.txt__ file  in created directory Ancestral.Allels with known ancestral alleles in positions.

If you working with hg19 the list of acestral allels could be extract from vcf [1000GP panel][1]. Run
```bash
./ancestral.alleles.sh 22 path.to/1000GP.chr22
```




## Merging 1000GP  and Archaic genomes

If you would like to work with 1000GP and archaic samples only we propose you pipeline briefly. 


Download [1000GP panel][1] and  archaic samples  [Link1][2] and [Link2][3]. Make .txt files with samples' names  obs.samples.txt, outgroup.txt, archaic.txt

Add full path to files  of 1000GP,  Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to parameters 1000GP and n1, n2, n3  and run 

```bash
./new.panel.preparation.Linux.sh 22 .path.to/outgroup.list path.to/observables.list path.to/file.bed 1000GP n1 n2 n3 all.chr22.vcf.gz
```
 
The resulting vcf.gz file is all.chr22.vcf.gz{.tbi}













## Make summary file 

You need  vcf file, lists of samples obs.samples.txt, outgroup.txt, archaic.txt and file with ancestral alleles positions
 to run  

```bash
./new.make.obs.sh 22 all.chr22.vcf.gz path.to/observables.list path.to/outgroup.list path.to/archaic.list ./Ancestral.Alleles/hg19.AA.chr22.txt path.to/file.bed
```





[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/

