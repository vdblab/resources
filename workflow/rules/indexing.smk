import os
import sys
from pathlib import Path


refs = Path("references") / "sequences"
indx = Path("references") / "indexes"



rule minimap:
    """
    we use Minimap to check for host contamination.
    """
    input:
        ref=f"{refs}/{{org}}/{{name}}/{{thisref}}.gz",
    output:
        ref=f"{indexes}/{{org}}/{{name}}/minimap2/{{thisref}}.gz.mmi",
    container:
        "docker://evolbioinfo/minimap2:v2.24"
    resources:
        mem_mb=64 * 1024,
    threads: 16
    shell:"""
    minimap2  -x sr --sr -u both --secondary=no -N 30 -c -d {output.ref} {input.ref}
    """
