
### Getting data and building Sylph db


The AGORA2 metadata list ~7k genomes, which are a mixure of genbank, refseq and unassambled reads. the "Genome Link" column has some accessions with full FTP links, wheras others just have the link to the assembly webpage.  The script `./scripts/split_accs_and_sras.R` fetches and cleans up the file, outputting two files: `agora2_accs` which contains the assembly accesions, and `agora2_sras` which contains the SRA accessions.

```
 snakemake --snakefile ./workflow/get_data_full.smk --config agora2_sras=$PWD/agora2_sras agora2_accs=$PWD/agora2_accs    --directory /data1/collab004/resources/dbs/sylph/agora2_clean_full/ --rerun-triggers mtime  --rerun-incomplete
 ```

## Note on Failures

<!-- Several of the genomes failed to download (eg Ruminiclostridium_thermocellum_DSM_2360, Yersinia_pseudotuberculosis_PB1).  those were found to be in "suppressed" status according to NCBI. There are other entries that lack an assembly accession; this workflow does not attempt to fetch those raw data and assemble.  In our hands we can fetch 6289 of the 6541 entries. -->
this successfully downloads 7182 genomes; some may have been removed from NCBI, causing the descrepency between that number and whats listed in the papers supplementary
