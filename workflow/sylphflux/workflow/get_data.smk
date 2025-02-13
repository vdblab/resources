import os
import pandas as pd

df = pd.read_excel(config["input_manifest"]).set_index("MicrobeID", drop=False)
print(df.shape)
df = df[~df["Genome link"].str.contains("fastq")].fillna("") # necessary for bugs like Prevotellaceae_bacterium_Marseille_P2826
print(df.shape)

localrules: all, sylphdb_tax_fmt

rule all:
    input:
        "agora2.syldb",
        "agora2.metadata.tsv",

def get_ftp_path(wc):
    pathbase = df.loc[wc.id, "Genome link"]
    return(
    pathbase + "/" + os.path.basename(pathbase) + "_genomic.fna.gz",
    pathbase + "/" + "assembly_status.txt"
    )

rule download:
    output:
        outf="genomes/{id}.fna.gz"
    params:
        ftppath = lambda wc: get_ftp_path(wc)[0],
        statuspath = lambda wc: get_ftp_path(wc)[1],
    shell: """
    if wget --spider {params.statuspath}  2>/dev/null ; then
        status=$(curl {params.statuspath} | head -n 1)
    	if [ "$status" == "status=suppressed" ]
    	then
	     echo "supressed record"
	     touch {output.outf}
    	else
             wget --continue --read-timeout=1  {params.ftppath} -O {output}
        fi
    else
        touch {output.outf}
	echo "bad url"
    fi


    """

rule sylphdb:
    input:
    	expand("genomes/{id}.fna.gz", id=df.MicrobeID.values),
    output:
        "agora2.syldb"
    container: "docker://ghcr.io/vdblab/sylph:0.6.1a"
    threads: 48
    shell: """
    sylph sketch genomes/*.fna.gz -t 48 -c 200 --out-name-db agora2
    """

rule sylphdb_tax_fmt:
    output:
        "agora2.metadata.tsv"
    run:
       df["taxout"] = "d__" + df["Kingdom"] + ";p__" + df["Phylum"] + ";c__" + df["Class"] + ";o__" + df["Order"]+ ";f__" + df["Family"]+ ";g__"+ df["Genus"] + ";s__" + df["Species"]
       df["pathid"] = df["MicrobeID"].astype(str) + ".fna.gz"
       df[["pathid", "taxout"]].to_csv(output[0], sep="\t", header=False, index=False)
