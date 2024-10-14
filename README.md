# vdblab databases
This snakemake workflow is used to create the database/reference/resource infrastructure to support the vdblab-pipelines.


## Usage
We recommend running in `-n` or `--dryrun` mode first to ensure you are only making needed changes; similarly, use `--rerun-triggers mtime` to prevent formatting changes or additional rules from automatically triggering rerunning of successful jobs.
```sh
snakemake  --directory /data/brinkvd/resources/ --rerun-triggers mtime -n
```


## Structure and naming
The main structure of the resulting folder has 4 subdirectories: **taxonomy**, **dbs**, **references**, and **indexes**. These represent the *types* of static resources used. The naming of the hierarchy should look like the following:


```
taxonomy/<name>/<version>/
# or
dbs/<name>/<version>/
# or
references/<organism>/<name>/<version>
# or
indexes/<organism>/<name>/<version>/<tool>
```

### Taxomony
This is for the files that define a taxonomic hierarchy.  Currenly we have the NCBI taxonomy downloaded here, but we could also add GTDB, SILVA, etc.
### dbs
Database archives should be listed in the `workflow/downloads.txt` file to be fetched. For anything more complex than unpacking, a new rule may need be created for that database.
### references
Reference sequences are arranged by organism, name, and version. When possible, retain the original file names from NCBI, so we can programatically access sequences and annotations with the same prefix. References should be listed in the `workflow/downloads.txt` file to be fetched, listing the appropriate fields for organism, name, version, and URL.
### indexes
The indexes are stored in a similar hierarchy to the reference sequences, with the addition of a tool-specific level.


# manual steps
Phanta must be downloaded/uploaded manually, as Dropbox is not accessible via the cluster.
