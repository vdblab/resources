import os
import pandas as pd

accsdf = pd.read_csv(config["agora2_accs"], sep="\t").set_index("MicrobeID", drop=False)
srasdf = pd.read_csv(config["agora2_sras"], sep="\t").set_index("MicrobeID", drop=False)

# accsdf["x"].to_csv("accs", header=False, index=False)
# accsdf.to_csv("agora2_accs"),
df = pd.concat([accsdf, srasdf])

localrules: all, sylphdb_tax_fmt, download


envvars:
    "NCBI_API_KEY"

rule all:
    input:
        "agora2.syldb",
        "agora2.metadata.tsv",
        expand("genomes/skesa/{id}.fna.gz", id=srasdf.MicrobeID.values),


rule download:
    input:
        inp=config["agora2_accs"]
    output:
        fnas = expand("genomes/raw/{id}.fna.gz", id=accsdf.MicrobeID.values),
   #     outdir = temp(directory("tmp"))
    resources:
        runtime = 1 * 60
    container: "docker://staphb/ncbi-datasets:18.9.0"
    shell: """
    #https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/
    cut -f 9 {input.inp} |tail -n+2 > accs
    datasets download genome accession --inputfile accs --dehydrated --filename tmp.zip
    unzip tmp.zip -d tmp
    mkdir -p genomes/raw/
    datasets rehydrate --gzip --directory tmp/
    cut -f1,9 {input.inp} | tail -n+2 | while read id accession ; do target="tmp/ncbi_dataset/data/${{accession}}/" ; if [ -d "$target" ] ; then   cp ${{target}}*.fna.gz    genomes/raw/${{id}}.fna.gz ; else touch  genomes/raw/${{id}}.fna.gz ; fi ; done
    """
# rule rename:
#     input:
#         fa=lambda wildcards: "genomes/{id}.fna.gz".format(id= accsdf.loc[(wildcards.id, wildcards.assembler), "R1s"].split(","))
#     	expand(, id=accsdf.x.values),
#     output:
#     	expand("genomes_renamed/{id}.fna.gz", id=accsdf[accsdf.MicrobeID == id),
#     shell:
#         """
#         cp {input} {output}
#         """

rule skesa:
    output:
    	"genomes/skesa/{id}.fna.gz"
    params:
        sra = lambda wildcards: srasdf.loc[wildcards.id, "x"]
    container: "docker://staphb/skesa:latest"
    resources:
        mem_mb = 32 * 1024,
        runtime = 2 * 60
    threads: 8

    shell:"""
    skesa --sra_run {params.sra} --cores {threads} --memory 30 | gzip > {output}

    """



rule sylphdb:
    input:
        expand("genomes/raw/{id}.fna.gz", id=accsdf.MicrobeID.values),
        expand("genomes/skesa/{id}.fna.gz", id=srasdf.MicrobeID.values)
    output:
        "agora2.syldb"
    container: "docker://ghcr.io/vdblab/sylph:0.6.1a"
    threads: 48
    shell: """
    sylph sketch genomes/*/*.fna.gz -t 48 -c 200 --out-name-db agora2
    """

rule sylphdb_tax_fmt:
    output:
        "agora2.metadata.tsv"
    run:
       df["taxout"] = "d__" + df["Kingdom"] + ";p__" + df["Phylum"] + ";c__" + df["Class"] + ";o__" + df["Order"]+ ";f__" + df["Family"]+ ";g__"+ df["Genus"] + ";s__" + df["Species"]
       df["pathid"] = df["MicrobeID"].astype(str) + ".fna.gz"
       df[["pathid", "taxout"]].to_csv(output[0], sep="\t", header=False, index=False)
