# vdblab databases
This snakemake workflow is used to create the database/reference/resource infrastructure to support the vdblab-pipelines.


## Usage
```sh
snakemake --profile /home/watersn/GitHub/vdblab-pipelines/msk-lsf --directory /path/to/new/resources/
``


## Structure and naming
The main structure of the resulting folder has 3 subdirectories: **taxonomy**, **dbs**, and **references**.  These represent the *types* of static resources used. The naming of the hierarchy should look like the following:


```
taxonomy/<name>/<version>/
# or
dbs/<name>/<version>/
# or
references/<organism>/<name>/<version>
```

### taxomony
### dbs
### references
