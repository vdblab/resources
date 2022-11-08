


# these aren't provided compressed :(
paths = {"CMMCv1": "https://ezmeta.unige.ch/CMMG/Kraken2db/cmmg/",
         "UHGGv1": "https://ezmeta.unige.ch/CMMG/Kraken2db/uhgg/"}
all_kraken_untarred_dbs =  expand(
    f"{dbs}/{{name}}/{{version}}/{{name}}.db_ready", zip,
    org=["mouse", "human"], name=["kraken", "kraken"],
    version=["CMMCv1", "UHGGv1"]
)


rule get_ftp_dir:
    container: default_container
    output:
        f"{dbs}/{{name}}/{{version}}/{{name}}.db_ready",
    params:
        url=lambda wildcards: paths[wildcards.version],
        cuts = lambda wildcards: len(paths[wildcards.version].replace("https://ezmeta.unige.ch/", "").split("/"))
    shell:"""
    cd $(dirname {output})
    wget --continue --no-parent  -nH --cut-dirs={params.cuts}  --read-timeout=1 --recursive {params.url}
    touch $(basename {output})
    """


rule get_phanta_db:
    """TODO: download both the normal and the prokaryote masked dbs from
    somewhere other than dropbox.  As of 2022-10-17, only the normal db is available
    (on dropbox), and the masked db has not been released.
    """
    container: default_container
    output:
        f"{dbs}/phanta/v1/phanta.db_ready",
    shell:"""
    wget https://github.com/bhattlab/phanta#advanced-usage
    cd $(dirname {output})
    wget --continue --read-timeout=1 -r wget https://www.dropbox.com/sh/3ktsdqlcph6x95r/AACGSj0sxYV6IeUQuGAFPtk8a/database_V1.tar.gz

    tar xvzf database_V1.tar.gz
    touch $(basename {output})
    """
