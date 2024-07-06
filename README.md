
# DAIseg
DAIseg method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes.  

The full description of the used files is in __File.description.md__ 
Run DAIseg if you have the files obs.samples.txt, region.bed, par.file.txt, arch.covering.chr22.txt, allels.ref.and.obs.chr22.txt. The premaded par.file.txt is in main directory which is corresponds to the scenario of Neanderthal introgression into the Europeans. 

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





     















[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
[4]: https://drive.google.com/file/d/1Vw-QEG9uu1trkbGHpDVXhMlbGt-RQhbN/view?usp=sharing
