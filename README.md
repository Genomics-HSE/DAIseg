
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes.  



The description of the used files is [here][1]




# Preparations of files

To read more details for files preparation see [readme][2]. To avoid details use script 
```bash
 ./full.preparation.sh hg19 22 file.bed n1 n2 n3 1000GP path.to.outgroup.list  path.to.obserables.list path.to.archaic.list name.out_vcf name.out_txt
```
where 1000GP, n1, n2, n3 are [vcf with 1000 Genome project][3], [neanderthal vcfs][4] and [chagyrskaya][5] (needed to be splited). 


There are two options without EM-algorithm and with EM algorithm. 


# Run DAI.seg without EM algorithm



```bash
python3 dai.seg.2.py --obs_samples ./samples/obs.samples.txt --bed ./hg19/regions/chr22.hg19.bed   --HMM_par par.file.txt --EM no --prepared_file ./hg19/allels.ref.and.obs.chr22.txt --o out.chr22.txt --arch_cover ./hg19/arch.covering.chr22.txt
```


# Run DAI.seg using EM algorithm

```bash
python3 dai.seg.2.py --obs_samples ./samples/obs.samples.txt --bed ./hg19/regions/chr22.hg19.bed   --HMM_par par.file.txt --EM yes --EM_steps 20  --prepared_file ./hg19/allels.ref.and.obs.chr22.txt --o out.EM.txt --arch_cover ./hg19/arch.covering.chr22.txt
```
to obtain estimations of the  coalescent times and run DAIseg. Here par.file.txt is used as the initial guess for EM algorithm.





[1]: https://github.com/Genomics-HSE/DAIseg/blob/main/File.types.md
[2]: https://github.com/Genomics-HSE/DAIseg/blob/main/hg19/README.md













[3]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[4]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[5]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/

