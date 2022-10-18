import os
import sys
from pathlib import Path


refs = Path("references")
dbs = Path("dbs")
tax = Path("taxonomy")



all_metaphlan =  multiext(str(dbs / "metaphlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103" / "mpa_vJan21_CHOCOPhlAnSGB_202103"), ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
all_metaerg = dbs / "metaerg" / "2022" / "db" / "blast" / "silva_LSURef.fasta"

all_card =    dbs / "CARD" / "v3.4.5/" / "card.json"

# minomer: includes prebuilt kraken :/
all_humann_dbs =  expand(f"{dbs}/{{name}}/{{version}}/{{base}}.db_ready", zip,
                         org=DB_MANIFEST["org"], name=DB_MANIFEST["name"], base=DB_MANIFEST["base"],
                         version=DB_MANIFEST["version"]
                         )


default_container = "docker://ghcr.io/vdblab/utility:0b"


rule metaphlan4:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command to
    metaphlan --install does more than just downloading the db, it also reformats fasta files, extracts files, runs bowtie to index, and more
    https://github.com/biobakery/MetaPhlAn/blob/3e818b755d2f084d4df31004e83c7cdd5c97a093/metaphlan/__init__.py#L87

    """
    output:
      #  outdir = directory(dbs / "metaphlan_chocophlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103"  ),
        all_metaphlan,
        #fasta = dbs / "metaphlan_chocophlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103" / "mpa_v30_CHOCOPhlAn_201901.fna.bz2",
    threads: 32
    resources:
        mem_mb=32*1024,
        runtime="6:00",
    params:
        dbdir = os.path.dirname(all_metaphlan[0]),#lambda wildcards, output: os.path.dirname(oup
    container: "docker://ghcr.io/vdblab/biobakery-profiler:20221001",
    shell:"""
    cd {params.dbdir}
    wget  --continue --read-timeout=1  http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar
    metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $PWD
    rm mpa_vJan21_CHOCOPhlAnSGB_202103.tar
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

rule metaerg_download:
    threads: 2
    output:
        dbs / "metaerg" / "2022" / "db.tar.gz"
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    container: default_container
    params:
        prefix=dbs / "metaerg" / "2022"
    shell: """
    wget http://ebg.ucalgary.ca/metaerg/db.tar.gz -P {params.prefix}
    """

rule metaerg_unpack:
    """why a whole separate rule? The download is ~23 gb, and unpacking takes ages
    """
    threads: 2
    input:
        dbs / "metaerg" / "2022" / "db.tar.gz"
    output:
        blast = directory(dbs / "metaerg" / "2022" / "db" /  "blast"),
        blast_silva = dbs / "metaerg" / "2022" / "db" /  "blast" / "silva_LSURef.fasta",
        diamond = directory(dbs / "metaerg" / "2022" /  "db" /"diamond"),
        hmm = directory(dbs / "metaerg" / "2022" /  "db" /"hmm"),
        sqlite3 = directory(dbs / "metaerg" / "2022" /  "db" /"sqlite3"),
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    params:
        prefix=dbs / "metaerg" / "2022"
    container: default_container
    shell: """
    cd {params.prefix}
    tar xzf db.tar.gz
    rm db.tar.gz
    """

rule prepare_humann_db:
    """
    THe server hosting these files is, again, rather unreliable and slow.
    We use wget instead of their humann_database --download tool which just downloads and extracts the gzipped tars
    See: https://github.com/biobakery/humann/blob/b1b674c122b5a538200bc3ba6a8efdade15ff496/humann/utilities.py#L510

    """
    input:
        #expand(f"{dbs}/{{name}}/{{version}}/{{base}}", zip, name=DB_MANIFEST["name"], version=DB_MANIFEST["version"], base=DB_MANIFEST["base"])
        f"{dbs}/{{name}}/{{version}}/{{base}}"
    params:
        basename=lambda wildcards, input: os.path.basename(input[0]),
        dirname=lambda wildcards, input: os.path.dirname(input[0]),
    output:
        f"{dbs}/{{name}}/{{version}}/{{base}}.db_ready",
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    container: default_container
    shell:"""
    tar xzf {input[0]} --directory {params.dirname}
    echo "{wildcards.name} db version {wildcards.version} prepared on $(date) by $(whoami)" > {output[0]}
    rm {input[0]}
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
# #

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
        raw = dbs / "CARD" / "v3.4.5/" / "card.json" ,
    resources:
        mem_mb=2*1024,
        runtime= "6:00",
    container:
        "docker://ghcr.io/vdblab/rgi:6.0.0"
    params:
        dbdir=lambda wildcards, output: os.path.dirname(output.raw),
        wildcards_db_version="v4.0.0",
        CARD_version = "v3.2.5",
    shell: """
    set -eux
    cd {params.dbdir}
    echo "Downloading CARD database; this happens the first time RGI is run in a new conda env"
    wget https://card.mcmaster.ca/download/0/broadstreet-{params.CARD_version}.tar.bz2 -O data #https://card.mcmaster.ca/latest/data
    echo "uncompressing data"
    tar -xvf data ./card.json
    rgi card_annotation -i card.json

    echo "Downloading variants"
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/download/6/prevalence-{params.wildcards_db_version}.tar.bz2 #https://card.mcmaster.ca/latest/variants
    mkdir -p wildcard
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    echo "unpacking variants"
    gunzip wildcard/*.gz

    rgi wildcard_annotation -i wildcard --card_json card.json -v {params.wildcards_db_version}

    """
