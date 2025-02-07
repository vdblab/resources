import os
import sys
from datetime import date
from pathlib import Path


refs = Path("references")
dbs = Path("dbs")
tax = Path("taxonomy")

localrules:
    get_mkspike_genomes

# changed from all_refseq to _refseq to prevent this from being run every single time the workflow is executed
today = date.today()
_refseq = f"{refs}/refseq/refseq/refseq/{today}-refseq.done"

all_metaphlan_21 = multiext(str(dbs / "metaphlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103" / "mpa_vJan21_CHOCOPhlAnSGB_202103"),
                            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
all_metaphlan_24 = multiext(str(dbs / "metaphlan" / "mpa_vJun23_CHOCOPhlAnSGB_202403" / "mpa_vJun23_CHOCOPhlAnSGB_202403"),
                            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")

all_metaerg = dbs / "metaerg" / "2022" / "db" / "blast" / "silva_LSURef.fasta"

all_card =    dbs / "CARD" / "v3.2.5/" / "card.json"

all_shortbred = dbs / "shortbred" / "colibactin" / "colibactin_markers_20240715.fa"
all_sylph =  dbs / "sylph" / "20250205-3kingdom/" / "3kingdom_metadata.tsv"

# formerly called all_humann_dbs
TAR_DB_MANIFEST =  DB_MANIFEST[DB_MANIFEST.url.str.contains("tar|tgz", regex=True)]
NOTAR_DB_MANIFEST =  DB_MANIFEST[~DB_MANIFEST.url.str.contains("tar|tgz", regex=True)]

# print(DB_MANIFEST.shape)
# print(TAR_DB_MANIFEST.shape)
# print(NOTAR_DB_MANIFEST.shape)

all_tarred_dbs =  expand(f"{dbs}/{{name}}/{{version}}/{{base}}.db_ready", zip,
                         org=TAR_DB_MANIFEST["org"], name=TAR_DB_MANIFEST["name"], base=TAR_DB_MANIFEST["base"],
                         version=TAR_DB_MANIFEST["version"]
                         )

# these are to trigger those files that just need to be downloaded and not unpacked
all_untarred_dbs =  expand(f"{dbs}/{{name}}/{{version}}/{{base}}", zip,
                         org=NOTAR_DB_MANIFEST["org"], name=NOTAR_DB_MANIFEST["name"], base=NOTAR_DB_MANIFEST["base"],
                         version=NOTAR_DB_MANIFEST["version"]
                         )

all_silva_db =  multiext(f"{dbs}/SILVA/138.1_SSURef_NR99/SILVA_138.1_SSURef_NR99_tax_silva", ".nsq", ".nin",  ".nhr")

all_spikes =  refs /  "mkspikeseq" / "mkspikeseq" / "0.0.0" / "spikes.done"

default_container = "docker://ghcr.io/vdblab/utility:0b"


rule metaphlan4_download:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command
    metaphlan --install does more than just downloading the db, it also reformats fasta files, extracts files, runs bowtie to index, and more
    https://github.com/biobakery/MetaPhlAn/blob/3e818b755d2f084d4df31004e83c7cdd5c97a093/metaphlan/__init__.py#L87

    """
    output:
        fasta = temp(dbs / "metaphlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103" / "mpa_vJan21_CHOCOPhlAnSGB_202103.tar"),
    threads: 1
    resources:
        runtime=12*60
    params:
        dbdir = lambda wc, output: os.path.dirname(output[0]),
    container: "docker://ghcr.io/vdblab/biobakery-profiler:20221001",
    shell:"""
    cd {params.dbdir}
    wget  --continue --read-timeout=1 http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar
    """

rule metaphlan4_download_2024:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command
    metaphlan --install does more than just downloading the db, it also reformats fasta files, extracts files, runs bowtie to index, and more
    https://github.com/biobakery/MetaPhlAn/blob/3e818b755d2f084d4df31004e83c7cdd5c97a093/metaphlan/__init__.py#L87

    """
    output:
        tar = temp(dbs / "metaphlan" / "mpa_vJun23_CHOCOPhlAnSGB_202403" / "mpa_vJun23_CHOCOPhlAnSGB_202403.tar"),
    threads: 1
    resources:
        runtime=12*60
    params:
        dbdir = lambda wc, output: os.path.dirname(output[0]),
    container: "docker://ghcr.io/vdblab/biobakery-profiler:20240513a",
    shell:"""
    cd {params.dbdir}
    wget  --continue --read-timeout=1  http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403.tar
    """

rule metaphlan4_process:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command to
    metaphlan --install does more than just downloading the db, it also reformats fasta files, extracts files, runs bowtie to index, and more
    https://github.com/biobakery/MetaPhlAn/blob/3e818b755d2f084d4df31004e83c7cdd5c97a093/metaphlan/__init__.py#L87

    """
    input:
        fasta = dbs / "metaphlan" / "mpa_vJan21_CHOCOPhlAnSGB_202103" / "mpa_vJan21_CHOCOPhlAnSGB_202103.tar",
    output:
        all_metaphlan_21,
    threads: 32
    resources:
        mem_mb=32*1024,
        runtime=12*60,
    params:
        dbdir = os.path.dirname(all_metaphlan_21[0]),#lambda wildcards, output: os.path.dirname(oup
    container: "docker://ghcr.io/vdblab/biobakery-profiler:20221001",
    shell:"""
    cd {params.dbdir}
    metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $PWD
    rm mpa_vJan21_CHOCOPhlAnSGB_202103.tar
    """

rule metaphlan4_process_2024:
    """
    This downloads the chocophlan db.  Their install code downloads very slowly
    and often fails so we get it with wget then rerun the install command to
    metaphlan --install does more than just downloading the db, it also reformats fasta files, extracts files, runs bowtie to index, and more
    https://github.com/biobakery/MetaPhlAn/blob/3e818b755d2f084d4df31004e83c7cdd5c97a093/metaphlan/__init__.py#L87

    """
    input:
        rules.metaphlan4_download_2024.output.tar
    output:
        all_metaphlan_24,
    threads: 32
    resources:
        mem_mb=32*1024,
        runtime=12*60,
    params:
        dbdir = os.path.dirname(all_metaphlan_24[0]),#lambda wildcards, output: os.path.dirname(oup
        indexname = os.path.basename(os.path.dirname(all_metaphlan_24[0]))
    container: "docker://ghcr.io/vdblab/biobakery-profiler:20240513a",
    shell:"""
    cd {params.dbdir}
    metaphlan --install --index {params.indexname} --bowtie2db $PWD
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
        runtime= 12*60,
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
        runtime= 12*60,
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
        runtime= 6*60,
    params:
        prefix=dbs / "metaerg" / "2022"
    container: default_container
    shell: """
    cd {params.prefix}
    tar xzf db.tar.gz
    rm db.tar.gz
    """

rule prepare_tarred_db:
    """
    This was originally for just the humann dbs but is used for iphop and kraken prebuilts as well
    THe server hosting these files is, again, rather unreliable and slow.
    We use wget instead of their humann_database --download tool which just downloads and extracts the gzipped tars
    See: https://github.com/biobakery/humann/blob/b1b674c122b5a538200bc3ba6a8efdade15ff496/humann/utilities.py#L510

    """
    input:
        f"{dbs}/{{name}}/{{version}}/{{base}}"
    params:
        basename=lambda wildcards, input: os.path.basename(input[0]),
        dirname=lambda wildcards, input: os.path.dirname(input[0]),
    output:
        f"{dbs}/{{name}}/{{version}}/{{base}}.db_ready",
    resources:
        mem_mb=2*1024,
        runtime= 12* 60,
    container: default_container
    shell:"""
    tar xzf {input[0]} --directory {params.dirname}
    echo "{wildcards.name} db version {wildcards.version} prepared on $(date) by $(whoami)" > {output[0]}
    rm {input[0]}
    """







# Keeping the following as a starting point should anyone attempt to make a kaiju db for here
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
        raw = dbs / "CARD" / "v3.2.5/" / "card.json" ,
    resources:
        mem_mb=2*1024,
        runtime=6*60,
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


rule get_quast_silva:
    output:
        silva=f"{dbs}/SILVA/138.1_SSURef_NR99/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
    shell:"""
    wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz -O {output.silva}
    wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz.md5  -O {output.silva}.md5
    cd $(dirname {output.silva})
    md5sum -c $(basename {output.silva}).md5
    """


rule format_quast_silva_blast_db:
    """ Quast wants the long-outdated old blast format that results in a nhr, nin, and nsq file
    """
    container:
        "docker://ncbi/blast:2.7.1"
        #"docker://nickp60/micro16s-blast",
    input:
        raw=f"{dbs}/SILVA/138.1_SSURef_NR99/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
    output:
        multiext(f"{dbs}/SILVA/138.1_SSURef_NR99/SILVA_138.1_SSURef_NR99_tax_silva", ".nsq", ".nin",  ".nhr")
    resources:
        mem_mb=8 * 1024,
        runtime=4*60,
    params:
        taxdir=tax / "NCBI" / "taxonkit",
        dbpath = lambda wildcards, input: input.raw.replace(".fasta.gz", ""),
        dbname = lambda wildcards, input: os.path.basename(input.raw.replace(".fasta.gz", "")),
    shell:"""
    set -eux
    gunzip -c {input.raw} | tr -d "_" | \
    makeblastdb  -in - -dbtype nucl -out {params.dbpath} -title {params.dbname}
    find dbs
    """

rule get_checkm:
    """ Depreciated; moved to a docker container with this built in
    """
    output:
        f"{dbs}/CHECKM/20150116/checkm_data_2015_01_16.tar.gz",
    shell:"""
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -O {output}
    cd $(dirname {output})
    tar xzf $(basename {output})
    find .
    """

rule get_refseq_genomes:
    # taken directly from https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/
    output:
        _refseq
    params:
        outdir=os.path.dirname(_refseq),
    resources:
        runtime=8*60,
    shell:"""
    cd {params.outdir}
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    awk -F "\t" '$12=="Complete Genome" &&  $11=="latest"{{print $20}}' assembly_summary.txt > ftpdirpaths
    awk 'BEGIN{{FS=OFS="/";filesuffix="genomic.fna.gz"}}{{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}}' ftpdirpaths > ftpfilepaths
    # get those not already here. note you will need to force rerunning if needed
    cat ftpfilepaths | while read f;
    do
    echo $f
    base=$(basename $f)
    if [ ! -f "$base" ]; then
    wget $f
    fi
    done
    ls . > $(basename {output})
    """

all_cazi_outputs =  dbs / "dbCAN2" / "v11/" / "stp.hmm.h3f"

rule cazi_db:
    container: "docker://haidyi/run_dbcan:3.0.1"
    resources:
        mem_mb=32*1024,
        runtime=12*60,
    output:
        trig = dbs / "dbCAN2" / "v11/" / "stp.hmm.h3f",
    params:
        db=lambda wildcards, output: os.path.dirname(output.trig)
    threads: 32
    shell:"""
    mkdir -p {params.db}
    cd {params.db}
    # don't ask why some are http and others https; I took this from the docs
    wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa && diamond makedb --in CAZyDB.08062022.fa -d CAZy --threads {threads}
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt && mv dbCAN-HMMdb-V11.txt dbCAN.txt && hmmpress dbCAN.txt
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb --threads {threads}
    wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm
    wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm
    """


rule hecatomb:
    """
    These files are pulled from the hecatomb database config file:
    https://github.com/shandley/hecatomb/blob/main/hecatomb/snakemake/config/dbFiles.yaml
    """
    threads: 2
    output: """dbs/hecatomb/v0/contaminants/nebnext_adapters.fa
       dbs/hecatomb/v0/contaminants/primerB.fa
       dbs/hecatomb/v0/contaminants/rc_primerB_ad6.fa
       dbs/hecatomb/v0/contaminants/truseq.fa
       dbs/hecatomb/v0/contaminants/vector_contaminants.fa
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB.dbtype
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_delnodes.dmp
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_h
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_h.dbtype
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_h.index
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB.index
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB.lookup
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_mapping
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_merged.dmp
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_names.dmp
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB_nodes.dmp
       dbs/hecatomb/v0/aa/virus_primary_aa/sequenceDB.source
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB.dbtype
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_delnodes.dmp
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_h
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_h.dbtype
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_h.index
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB.index
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB.lookup
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_mapping
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_merged.dmp
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_names.dmp
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB_nodes.dmp
       dbs/hecatomb/v0/aa/virus_secondary_aa/sequenceDB.source
       dbs/hecatomb/v0/host/bat/masked_ref.fa.gz
       dbs/hecatomb/v0/host/camel/masked_ref.fa.gz
       dbs/hecatomb/v0/host/cat/masked_ref.fa.gz
       dbs/hecatomb/v0/host/celegans/masked_ref.fa.gz
       dbs/hecatomb/v0/host/cow/masked_ref.fa.gz
       dbs/hecatomb/v0/host/dog/masked_ref.fa.gz
       dbs/hecatomb/v0/host/human/masked_ref.fa.gz
       dbs/hecatomb/v0/host/macaque/masked_ref.fa.gz
       dbs/hecatomb/v0/host/mosquito/masked_ref.fa.gz
       dbs/hecatomb/v0/host/mouse/masked_ref.fa.gz
       dbs/hecatomb/v0/host/rat/masked_ref.fa.gz
       dbs/hecatomb/v0/host/tick/masked_ref.fa.gz
       dbs/hecatomb/v0/host/virus_shred.fasta.gz
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB.dbtype
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB_h
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB_h.dbtype
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB_h.index
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB.index
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB.lookup
       dbs/hecatomb/v0/nt/virus_primary_nt/sequenceDB.source
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB.dbtype
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB_h
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB_h.dbtype
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB_h.index
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB.index
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB.lookup
       dbs/hecatomb/v0/nt/virus_secondary_nt/sequenceDB.source
       dbs/hecatomb/v0/tables/2020_07_27_Viral_classification_table_ICTV2019.txt
       dbs/hecatomb/v0/tables/phage_taxonomic_lineages_12102020.txt""".split()
    resources:
        mem_mb=2*1024,
        runtime= 12*60,
    params:
        prefix=dbs / "hecatomb" / "v0",
        url="https://hecatombdatabases.s3.us-west-2.amazonaws.com/databases/"
    container: default_container
    shell: """
    for path in {output}
    do
        urlbase=$(echo $path | sed "s|dbs/hecatomb/v0/||g")
        echo "fetching $path"
        wget  --continue --read-timeout=1 {params.url}${{urlbase}} -O $path
    done
    """

rule uncompress_uniref90:
    input:
        uniref_ref=dbs / "uniref90" / "20240715" / "uniref90.fasta.gz",
    output:
        uniref_ref=temp(dbs / "uniref90" / "20240715" / "uniref90.fasta"),
    resources:
        mem_mb= 64 * 1024,
        runtime= 8 * 60,
    shell: """    gunzip -c $PWD/{input.uniref_ref}  > {output.uniref_ref} """

rule build_shortbred_uniref90_blastdb:
    input:
        uniref_ref=dbs / "uniref90" / "20240715" / "uniref90.fasta",
    output:
        uniref_ref=dbs / "uniref90" / "20240715" / "uniref90.pal",
    params:
        outname = lambda wc, input: input.uniref_ref.replace(".fasta", "")
    container: "docker://biobakery/shortbred:0.9.5"
    resources:
        mem_mb= 64 * 1024,
        runtime= 8 * 60,
    shell: """
    makeblastdb -in {input.uniref_ref} -out {params.outname} -dbtype prot -logfile {params.outname}.log
    """

rule build_colibactin_shortbred_markers:
    input:
        target_fa=os.path.join(workflow.basedir, "../data/colibactin_genes.fasta"),
        uniref_ref=dbs / "uniref90" / "20240715" / "uniref90.pal",
#        usearch_path=config["usearch_path"]
    output:
        marker=dbs / "shortbred" / "colibactin" / "colibactin_markers_20240715.fa",
        tmp_dir=temp(directory("results/tmp_marker"))
    container: "docker://biobakery/shortbred:0.9.5"
    params:
        dbname = lambda wc, input: input.uniref_ref.replace(".pal", "")
    threads: 16
    resources:
        mem_mb= 32 * 1024,
        runtime= 2 * 60,
    shell:
        '''
        shortbred_identify.py \
            --goi {input.target_fa} \
            --refdb {params.dbname} \
            --markers {output.marker} \
            --threads {threads} \
            --tmp {output.tmp_dir}
        '''

spikes = ["Salinibacter_ruber", "Trichoderma_reesei", "Haloarcula_hispanica"]

rule get_mkspike_genomes:
    output:
        expand(f"{refs}/{{org}}/{{name}}/{{version}}/{{base}}.fa",
               org="mkspikeseq", name="mkspikeseq",
               base=spikes, version="0.0.0"),
        expand(f"{refs}/{{org}}/{{name}}/{{version}}/{{base}}_full.bed",
               org="mkspikeseq", name="mkspikeseq",
               base=spikes, version="0.0.0"),
        refs /  "mkspikeseq" / "mkspikeseq" / "0.0.0" / "spikes.done",
        table=refs /  "mkspikeseq" / "mkspikeseq" / "0.0.0" / "mkspike.csv",
    params:
        outdir = lambda wc, output: os.path.dirname(output[0]),
        abspath = os.path.abspath(os.getcwd()),
    # get a container for this, as it requires curl,wget, samtools
    shell:"""
    set -euo
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.fna.gz | gunzip > {params.outdir}/Salinibacter_ruber.fa
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip  > {params.outdir}/Trichoderma_reesei.fa
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/905/GCF_000223905.1_ASM22390v1/GCF_000223905.1_ASM22390v1_genomic.fna.gz | gunzip  > {params.outdir}/Haloarcula_hispanica.fa
    for fa in {params.outdir}/*fa;
    do
        echo $fa
        samtools faidx $fa
        sort ${{fa}}.fai  > ${{fa}}.fai.sorted
        thisbase=$(echo $fa | sed "s|.fa$||g")
        awk 'BEGIN {{FS="\t"}}; {{print $1 FS "0" FS $2}}' ${{fa}}.fai.sorted | sort -k 1,1 -k2,2n > ${{thisbase}}_full.bed
        echo -e "$(basename $thisbase),{params.abspath}/${{fa}},{params.abspath}/${{thisbase}}_full.bed" >> {output.table}
    done
    touch {params.outdir}/spikes.done
    """



rule download_sylph_dbs:
   output:
       dbs / "sylph" / "20250205-3kingdom/" / "fungi-refseq-2024-07-25-c200-v0.3.syldb",
       dbs / "sylph" / "20250205-3kingdom/" / "gtdb-r220-c200-dbv1.syldb",
       dbs / "sylph" / "20250205-3kingdom/" / "imgvr_c200_v0.3.0.syldb",
       dbs / "sylph" / "20250205-3kingdom/" / "fungi_refseq_2024-07-25_metadata.tsv.gz",
       dbs / "sylph" / "20250205-3kingdom/" / "gtdb_r220_metadata.tsv.gz",
       dbs / "sylph" / "20250205-3kingdom/" / "IMGVR_4.1_metadata.tsv.gz",
       dbs / "sylph" / "20250205-3kingdom/" / "3kingdom_metadata.tsv",
   params:
       outdir=lambda wc, output: os.path.dirname(output[0]),
   shell:
    """
    cd {params.outdir}
    wget http://faust.compbio.cs.cmu.edu/sylph-stuff/fungi-refseq-2024-07-25-c200-v0.3.syldb
    wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb
    wget http://faust.compbio.cs.cmu.edu/sylph-stuff/imgvr_c200_v0.3.0.syldb
    wget https://zenodo.org/records/14320496/files/fungi_refseq_2024-07-25_metadata.tsv.gz
    wget https://zenodo.org/records/14320496/files/gtdb_r220_metadata.tsv.gz
    wget https://zenodo.org/records/14320496/files/IMGVR_4.1_metadata.tsv.gz
    zcat fungi_refseq_2024-07-25_metadata.tsv.gz gtdb_r220_metadata.tsv.gz IMGVR_4.1_metadata.tsv.gz > 3kingdom_metadata.tsv
    """
