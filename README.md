
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 


# Run DAI.seg without EM algorithm

If you have the files obs.samples.txt, region.bed, par.file.txt, arch.covering.chr22.txt, allels.ref.and.obs.chr22.txt (the detailed description see in __File's summary__ section) 

```bash
python3 dai.seg.2.py --obs_samples obs.samples.txt --bed test.bed   --HMM_par par.file.txt --EM no --prepared_file allels.ref.and.obs.chr22.txt --o out.chr22.txt --arch_cover arch.covering.chr22.txt
```

where the examples par.file.txt could be found in the main directory

# Run DAI.seg using EM algorithm

```bash
python3 dai.seg.2.py --obs_samples obs.samples.txt --bed test.bed   --HMM_par par.file.txt --EM yes --EM_steps 20  --prepared_file allels.ref.and.obs.chr22.txt --o out.EM.txt --arch_cover arch.covering.chr22.txt
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.





     







# Whole pipeline

## Archaic covering 
The goal is to create __arch.covering.chr22.txt__ file with the window-covering by archaic samples. 
Add full path to files  of   Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables n1, n2, n3 and run 
```bash
./archaic.covering.sh 22 test.bed n1 n2 n3
```


## Ancestral alleles
The goal is to create __hg19.AA.chr22.txt__ file  in created directory Ancestral.Allels with known ancestral alleles in positions.

If you working with hg19 the list of acestral allels could be extract from vcf [1000GP panel][1]. Run
```bash
./ancestral.alleles.sh 22 way.to.1000GP.chr22
```




## Merging 1000GP  and Archaic genomes

If you would like to work with 1000GP and archaic samples only we propose you pipeline briefly. 


Download [1000GP panel][1] and  archaic samples  [Link1][2] and [Link2][3]. Make .txt files with samples' names  obs.samples.txt, outgroup.txt, archaic.txt

Add full path to files  of 1000GP,  Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to parameters 1000GP and n1, n2, n3  and run 

```bash
./new.panel.preparation.Linux.sh 22 Outgroup.txt obs.samples.txt test.bed 1000GP n1 n2 n3 all.chr22.vcf.gz
```
 
The resulting vcf.gz file is all.chr22.vcf.gz{.tbi}













## Make observations 

You need  vcf file, lists of samples obs.samples.txt, outgroup.txt, archaic.txt and file with [ancestral alleles positions][4]
 to run  

```bash
./new.make.obs.sh 22 all.chr22.vcf.gz obs.samples.txt Outgroup.txt archaic.txt  ./Ancestral.Alleles/hg19.AA.chr22.txt test.bed
```

and to make observation files obs.neand.chr22.txt, obs.outgroup.chr22.txt








[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
[4]: https://drive.google.com/file/d/1Vw-QEG9uu1trkbGHpDVXhMlbGt-RQhbN/view?usp=sharing
