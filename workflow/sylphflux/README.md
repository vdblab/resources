
### Getting data and building Sylph db

```
wget https://github.com/opencobra/COBRA.papers/raw/refs/heads/master/2021_demeter/input/AGORA2_infoFile.xlsx
snakemake --snakefile ./workflow/get_data.smk --config input_manifest=$PWD/AGORA2_infoFile.xlsx   --directory $PWD/db/
```

## Note on Failures

Several of the genomes failed to download (eg Ruminiclostridium_thermocellum_DSM_2360, Yersinia_pseudotuberculosis_PB1).  those were found to be in "suppressed" status according to NCBI. There are other entries that lack an assembly accession; this workflow does not attempt to fetch those raw data and assemble.  In our hands we can fetch 6289 of the 6541 entries.
