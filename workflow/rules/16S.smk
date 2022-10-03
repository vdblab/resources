import os
import sys
from pathlib import Path


refs = Path("references")
dbs = Path("dbs")
tax = Path("taxonomy")


db16S = multiext(f"{dbs}/ncbi16S/2022/16S_ribosomal_RNA", ".fna", ".ndb", ".nhr", ".nin", ".nnd", ".nni", ".nog", ".nos",  ".nsq",
                 ".ntf", ".nto",  "_id_and_taxonomy.txt")



default_container = "docker://ghcr.io/vdblab/utility:0b"

rule all:
    input:
        db16S,

rule install_taxonkit_taxonomy:
    """ I wish this wasn't needed, but taxonkit via singularity doesn't seem
    to be able to access where taxonkit installs the nodes db
    Further, taxonkit needs its OWN copy, they can't just
     use the unpacked taxdump.  Other tools might use the taxdump (eg kraken)
    so we leave it here
    """
    container: default_container
    output:
        dir= directory(tax / "NCBI"),
        tkdir= directory(tax / "NCBI" / "taxonkit" ),
        nodes = tax / "NCBI" / "nodes.dmp",
        tknodes = tax / "NCBI" / "taxonkit" / "nodes.dmp",

    shell:"""
    cd {output.dir}
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz
    cp  names.dmp nodes.dmp delnodes.dmp merged.dmp taxonkit
    rm taxdump.tar.gz
    """

rule ncbi_16S_get:
    """
    Where do the 16S database files come from?

    NCBI releases a prebuilt BLAST database for 16s RNA, but they don't include the fasta file.
    <https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz>

    So we create it from the released Refseq targetted loci releases

    """
    container:
        "docker://nickp60/micro16s-blast",
    output:
        f"{dbs}/ncbi16S/2022/16S_ribosomal_RNA.tar.gz",
    params:
        taxdir=tax / "NCBI" / "taxonkit",
    shell:"""
#    wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz
#    wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
#    zcat archaea.16SrRNA.fna.gz bacteria.16SrRNA.fna.gz > {output}
    wget  https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz -O {output}
    """



rule ncbi_16S_format:
    """

    that is then post-processed to get the labels file using `blastdbcmd` (from the blast commandline package) to get the accessions and taxids

    Then, because `blastdbcmd` doesn't give the whole taxonomy hierarchy, use taxonkit to extract the data in a way that works with the database
    """
    container:
        "docker://nickp60/micro16s-blast",
    input:
        raw=rules.ncbi_16S_get.output[0],
        taxdb_nodes = "taxonomy/NCBI/nodes.dmp",
    output:
        # it doens't work to just give db16S as the output; complains about a missing comma?
        multiext(f"{dbs}/ncbi16S/2022/16S_ribosomal_RNA", ".fna", ".ndb", ".nhr", ".nin", ".nnd", ".nni", ".nog", ".nos",  ".nsq",
                 ".ntf", ".nto",  "_id_and_taxonomy.txt")
    resources:
        mem_mb=8 * 1024,
        runtime="4:00",
    params:
        taxdir=tax / "NCBI" / "taxonkit",
        dbpath = lambda wildcards, input: input.raw.replace(".tar.gz", ""),
        dbname = lambda wildcards, input: os.path.basename(input.raw.replace(".tar.gz", "")),
    shell:"""
    # extra logging given the number of bash variables etc
    set -eux
    HERE=$PWD
    cd $(dirname {input.raw})
    tar xzf $(basename {input.raw})
    cd $HERE

    # taxonkit needs to know where to find taxdump from NCBI
    export TAXONKIT_DB=${{PWD}}/{params.taxdir}/
    export BLASTDB=$PWD/$(dirname {input.raw})

    # create a nice fasta from the database's sequences
    blastdbcmd -entry all -db {params.dbname} -out {output[0]}
    # extract an accession/taxid table
    blastdbcmd -db {params.dbname} -outfmt "%a %T" -entry "all" > acc_taxid

    # create the fasta/taxonomy labels file
    cut -f 2 -d" "  acc_taxid  | taxonkit lineage  | awk '$2!=""'  | taxonkit reformat -P | cut -f 3 | paste acc_taxid - > {params.dbpath}_id_and_taxonomy.txt
    rm acc_taxid
    find .
    """
