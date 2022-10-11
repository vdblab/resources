import os
import sys
from pathlib import Path


refs = Path("references")
idxs = Path("indexes")



rule minimap:
    """
    we use Minimap to check for host contamination.
    """
    input:
        ref=f"{refs}/{{org}}/{{name}}/{{thisref}}",
    output:
        ref=f"{idxs}/{{org}}/{{name}}/minimap2/{{thisref}}.mmi",
    container:
        "docker://evolbioinfo/minimap2:v2.24"
    resources:
        mem_mb=64 * 1024,
    threads: 16
    shell:"""
    minimap2  -x sr --sr -u both --secondary=no -N 30 -c -d {output.ref} {input.ref}
    """

bowtie_prebuilts = {
    "GRCm39": {"org": "mouse", "url": "https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip", "dest": f"{idxs}/human/GRCm39/bowtie/GRCm39/"},
    "CHM13":  {"org": "human", "url": "https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip", "dest": f"{idxs}/human/CHM13/bowtie/CHM13/"}
}

rule get_prebuild_bowties:
    output:
        ref=directory(expand("{dest}", dest = [v["dest"] for k, v in bowtie_prebuilts.items()])),
    params:
        url=expand("{url}", url = [v["url"] for k, v in bowtie_prebuilts.items()]),
    shell:"""
    cd {output.ref}
    wget {params.url}
    unzip $(basename {params.url}
    """


rule snap:
    container:""
    input:
        ref=f"{refs}/{{org}}/{{name}}/{{thisref}}",
    output:
        ref=f"{idxs}/{{org}}/{{name}}/snap/{{thisref}}/Genome",
    threads: 16
    resources:
        mem_mb=64 * 1024
    script:"""
    zcat {input.ref}  > tmp.fa
    snap-aligner index tmp.fa chm13v2.0/ -t{threads}
    rm tmp.fa
    """
