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


# bowtie_prebuilts = {
#     "GRCm39": {"org": "mouse",
#                "url": "https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip",
#                "dest": f"{idxs}/mouse/GRCm39/bowtie/GRCm39/"},
#     "CHM13":  {"org": "human",
#                "url": "https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip",
#                "dest": f"{idxs}/human/CHM13/bowtie/CHM13/"}
# }

# all_bowtie_prebuit_triggers = [
#     f"{idxs}/human/CHM13/v1.0/bowtie/chm13.draft_v1.0_plusY.1.bt2",
#     f"{idxs}/mouse/GRCm39/GCA_000001635.27/bowtie/GRCm39.1.bt2"]


# # I was trying to be clever and use wildcards etc but this was easier to hardcode two rules :(
# rule get_prebuild_bowtie_mouse:
#     output:
#         ref=directory(f"{idxs}/mouse/GRCm39/GCA_000001635.27/bowtie/"),
#         index_1bt2 = f"{idxs}/mouse/GRCm39/GCA_000001635.27/bowtie/GRCm39.1.bt2",
#     params:
#         url="https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip"
#     shell:"""
#     mkdir -p {output.ref}
#     cd {output.ref}
#     wget {params.url}
#     # the -j junks the paths, so stuff ends up in the current working dir
#     unzip -j  $(basename {params.url})
#     find .
#     """

# rule get_prebuild_bowtie_human:
#     output:
#         ref=directory(f"{idxs}/human/CHM13/v1.0/bowtie/"),
#         index_1bt2 = f"{idxs}/human/CHM13/v1.0/bowtie/chm13.draft_v1.0_plusY.1.bt2",
#     params:
#         url="https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip"
#     shell:"""
#     mkdir -p {output.ref}
#     cd {output.ref}
#     wget {params.url}
#     # the -j junks the paths, so stuff ends up in the current working dir
#     unzip -j  $(basename {params.url})
#     find .
#     """


rule snap:
    container: "docker://ghcr.io/vdblab/snap-aligner:2.0.1"
    input:
        ref=f"{refs}/{{org}}/{{name}}/{{version}}/{{thisref}}",
    output:
        ref=f"{idxs}/{{org}}/{{name}}/{{version}}/snap/{{thisref}}/Genome",
    threads: 16
    resources:
        mem_mb=64 * 1024
    params:
        dbdir=lambda wildcards, output: os.path.dirname(output.ref),
    shell:"""
    snap-aligner index {input.ref} {params.dbdir}/ -t{threads}
    """

rule bowtie_index:
    container: "docker://staphb/bowtie2:2.4.4"
    input:
        ref=f"{refs}/{{org}}/{{name}}/{{version}}/{{thisref}}",
    output:
        ref=f"{idxs}/{{org}}/{{name}}/{{version}}/bowtie2/{{thisref}}.1.bt2",
    threads: 32
    resources:
        mem_mb=64 * 1024
    params:
        dbpre=lambda wildcards, output: output.ref.replace(".1.bt2", "")
    shell:"""
    bowtie2-build --threads {threads} {input.ref} {params.dbpre}
    """
