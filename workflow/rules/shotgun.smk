import os
import sys
from pathlib import Path


refs = Path("references")
dbs = Path("dbs")
tax = Path("taxonomy")



metaphlan =  multiext("mpa_vJan21_CHOCOPhlAnSGB_202103", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")


default_container = "docker://ghcr.io/vdblab/utility:0b"

rule all:
    input:
        metaphlan,

rule metaphlan:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command to
    continue the rest of the program
    cd /data/brinkvd/resources/biobakery_workflows_dbs/metaphlanV4/ &&
    metaphlan --install --bowtie2db /data/brinkvd/resources/biobakery_workflows_dbs/metaphlanV4/

    """
    threads: 8
    output:
        db=multiext(dbs / "metaphlan" / "202103" / "mpa_vJan21_CHOCOPhlAnSGB_202103", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
    resources:
        mem_mb=32*1024,
        runtime= "6:00",
    params:
        prefix=dbs / "metaphlan" / "202103"
    shell: """
    cd {params.prefix}
    # you can get this path by running metaphlan --install and watching the logs.
    wget -c  --read-timeout=5 http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar
    metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $PWD"
    """


rule antismash:
    """
    This creates many files, but we add just a few for the output
    """
    threads: 2
    output:
        clusters=dbs / "antismash" / "6.0" / "clusterblast" / "clusters.txt",
        proteins=dbs / "antismash" / "6.0" / "clusterblast" / "proteins.fasta",
        resfam = directory(dbs / "antismash" / "6.0" / "resfam"),
        tigrfam = directory(dbs / "antismash" / "6.0" / "tigrfam"),
        pfam = directory(dbs / "antismash" / "6.0" / "pfam"),
        clustercompare = directory(dbs / "antismash" / "6.0" / "clustercompare"),
    resources:
        mem_mb=2*1024,
        runtime= "12:00",
    container:
        "docker://antismash/standalone:6.0.0"
    params:
        prefix=dbs / "antismash" / "6.0"
    shell: """
    download-antismash-databases --database-dir {params.prefix}
    """

rule metaerg:
    threads: 2
    output:
        blast = directory(dbs / "metaerg" / "2022" / "blast"),
        diamond = directory(dbs / "metaerg" / "2022" / "diamond"),
        hmm = directory(dbs / "metaerg" / "2022" / "hmm"),
        sqlite3 = directory(dbs / "metaerg" / "2022" / "sqlite3"),
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    container:
        "docker://antismash/standalone:6.0.0"
    params:
        prefix=dbs / "metaerg" / "2022"
    shell: """
    wget http://ebg.ucalgary.ca/metaerg/db.tar.gz -P {params.prefix}
    cd {params.prefix}
    tar xzf db.tar.gz
    """



# ### Metaphlan

# ```
# mkdir /data/brinkvd/resources/biobakery_workflows_dbs/v31
# # wget -P /data/brinkvd/resources/biobakery_workflows_dbs/v31 http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v31_CHOCOPhlAn_201901_marker_info.txt.bz2

# singularity run --bind /data/brinkvd/ /data/brinkvd/.singularity/d11b613c13947f4c48596feeca9d18cc.simg metaphlan  --install --index mpa_v31_CHOCOPhlAn_201901 --bowtie2db /data/brinkvd/resources/biobakery_workflows_dbs/v31/

# $ tree /data/brinkvd/resources/biobakery_workflows_dbs/v31/
# /data/brinkvd/resources/biobakery_workflows_dbs/v31/
# ├── mpa_v31_CHOCOPhlAn_201901.1.bt2
# ├── mpa_v31_CHOCOPhlAn_201901.2.bt2
# ├── mpa_v31_CHOCOPhlAn_201901.3.bt2
# ├── mpa_v31_CHOCOPhlAn_201901.4.bt2
# ├── mpa_v31_CHOCOPhlAn_201901.fna.bz2
# ├── mpa_v31_CHOCOPhlAn_201901_marker_info.txt.bz2
# ├── mpa_v31_CHOCOPhlAn_201901.md5
# ├── mpa_v31_CHOCOPhlAn_201901.pkl
# ├── mpa_v31_CHOCOPhlAn_201901.rev.1.bt2
# ├── mpa_v31_CHOCOPhlAn_201901.rev.2.bt2
# └── mpa_v31_CHOCOPhlAn_201901.tar

# ```


# """
# #### Version 4
# ### Humann
# # This is uneeded (see https://forum.biobakery.org/t/announcing-metaphlan-3-1-and-humann-3-1/3881)
# #singularity run --bind /data/brinkvd/ /data/brinkvd/.singularity/d11b613c13947f4c48596feeca9d18cc.simg humann_databases --download uniref uniref90_diamond /data/brinkvd/resources/biobakery_workflows_dbs/v31/uniref90_diamond/

# # this stalls becuase the biobakery servers are super unreliable
# # singularity run --bind /data/brinkvd/ /data/brinkvd/.singularity/d11b613c13947f4c48596feeca9d18cc.simg humann_databases --download chocophlan full /data/brinkvd/resources/biobakery_workflows_dbs/v31/choco/
# # have to use wget instead
# wget -c --read-timeout=.1 http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz -P /data/brinkvd/resources/biobakery_workflows_dbs/v31/choco/

# cd /data/brinkvd/resources/biobakery_workflows_dbs/v31/choco/ && tar xzf full_chocophlan.v201901_v31.tar.gz && cd ../ && mv choco choco_v201901_v31


# ### Human

# """





# ## kaiju

# ```
# cd /data/brinkvd/resources/kaiju
# conda create -n kaiju kaiju
# conda activate kaiju
# bsub -n 10 -e dload_kaiju_db.e -o dload_kaiju_db.o "kaiju-makedb -t 10 -s nr_euk"



# # that died partway, retstarted manually the commands found here https://github.com/bioinformatics-centre/kaiju/blob/591302244a95e8410e5f8df2b95bc21392611e0d/util/kaiju-makedb#L230
# export DB=nr_euk
# wget -c ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
# wget -c  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

# gunzip  $DB/prot.accession2taxid.gz

# # submit the job to format the AA file
# bsub -n 12 -R "rusage[mem=4GB]" -e convertnr.e -o convertnr.o "gunzip -c $DB/nr.gz | kaiju-convertNR  -t nodes.dmp -g $DB/prot.accession2taxid -e /home/watersn/miniconda3/envs/kaiju/bin/kaiju-excluded-accessions.txt -a -o $DB/kaiju_db_$DB.faa -l /home/watersn/miniconda3/envs/kaiju/bin/kaiju-taxonlistEuk.tsv"

# # Submit the job to make the BWT index
# bsub -n 12 -R "rusage[mem=16GB]" -W12:00 -e makebwt.e -o makebwt.o "kaiju-mkbwt -e 5 -n $12 -a ACDEFGHIKLMNPQRSTVWY -o $DB/kaiju_db_$DB $DB/kaiju_db_$DB.faa"

# # And finally submit the job to build the FM index file
# bsub -n 12 -R "rusage[mem=16GB]" -W12:00 -e makefm.e -o makefm.o "kaiju-mkfmi $DB/kaiju_db_$DB"


# ```


rule download_and_format_CARD:
    threads: 2
    output:
        blast = directory(dbs / "metaerg" / "2022" / "blast"),
        diamond = directory(dbs / "metaerg" / "2022" / "diamond"),
        hmm = directory(dbs / "metaerg" / "2022" / "hmm"),
        sqlite3 = directory(dbs / "metaerg" / "2022" / "sqlite3"),
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    container:
        "docker://ghcr.io/vdblab/rgi:6.0.0"
    params:
        prefix=dbs / "CARD" / "v3.2.5"
    shell: """
    cd {params.prefix}
    echo "Downloading CARD database; this happens the first time RGI is run in a new conda env"
    wget https://card.mcmaster.ca/download/0/broadstreet-v3.2.5.tar.bz2 -O data #https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    rgi load --card_json ./card.json
    rgi card_annotation -i card.json > card_annotation.log 2>&1
    CARD_VER=$(rgi database --version)
    rgi load -i card.json --card_annotation card_database_v${CARD_VER}.fasta
    echo "Downloading variants"
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/download/6/prevalence-v4.0.0.tar.bz2 #https://card.mcmaster.ca/latest/variants
    mkdir -p wildcard
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    echo "unpacking variants"
    gunzip wildcard/*.gz

    # these get run on execution currently; not the most efficient but I can't
    # seem to point rgi to an existing DB that isn't in the default search path
    #rgi wildcard_annotation -i wildcard --card_json card.json  -v $CARD_VER > wildcard_annotation.log 2>&1

    #rgi load --wildcard_annotation wildcard_database_v${CARD_VER}.fasta \
        --wildcard_index wildcard/index-for-model-sequences.txt \
	--card_annotation card_database_v${CARD_VER}.fasta --local
    """



## Bowtie indexes

```
mkdir bowtie_indexes && cd
wget  https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip
wget https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
unzip chm13.draft_v1.0_plusY.zip
unzip GRCm39.zip

```

## snap indexes
```

mkdir snap_indexes
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip chm13v2.0.fa.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz

(pull docker://ghcr.io/vdblab/snap-aligner:2.0.1; it saved as /data/brinkvd/.singularity/5b0edf82245696a61e21418a09f2e963.simg)
bsub -n 16 -R "rusage[mem=4]" -e buildchm13v2.e -o buildchm13v2.o  "module load singularity/3.7.1 && singularity run -B $PWD /data/brinkvd/.singularity/5b0edf82245696a61e21418a09f2e963.simg snap-aligner index chm13v2.0.fa chm13v2.0/ -t16"
bsub -n 16 -R "rusage[mem=4]" -e buildGRCm39.e -o buildGRCm39.o  "module load singularity/3.7.1 && singularity run -B $PWD  /data/brinkvd/.singularity/5b0edf82245696a61e21418a09f2e963.simg snap-aligner index GCF_000001635.27_GRCm39_genomic.fna GRCm39/ -t16"


```
