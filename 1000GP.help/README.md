these files 

1.1 obtain list of samples by population for 1000GP.grch38

```bash
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel


./make.pop.list.unrelated.grch38.sh IBS 
```


1.2  obtain list of unrelated samples for 1000GP.grch37
```bash
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped



./make.pop.list.unrelated.grch37.sh IBS 
```

  
2. generate jsons
...
