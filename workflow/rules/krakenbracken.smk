
# these aren't provided compressed :(
paths = {"CMMCv1": "https://ezmeta.unige.ch/CMMG/Kraken2db/cmmg/",
         "UHGGv1": "https://ezmeta.unige.ch/CMMG/Kraken2db/uhgg/"}
all_kraken_untarred_dbs =  expand(
    f"{dbs}/{{name}}/{{version}}/{{name}}.db_ready", zip,
    org=["mouse", "human"], name=["kraken", "kraken"],
    version=["CMMCv1", "UHGGv1"]
)

bracken_readlens=[50,75,100,150,200,250,300]

all_bracken_dbs = expand(
    f"{dbs}/kraken/{{mock}}/database{{readlens}}mers.kmer_distrib",
    readlens=bracken_readlens,
    mock=["zymogutv0", "zymogutv1_plusspikes"]

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


# rule get_phanta_db:
#     """TODO: download both the normal and the prokaryote masked dbs from
#     somewhere other than dropbox.  As of 2022-10-17, only the normal db is available
#     (on dropbox), and the masked db has not been released.
#     """
#     container: default_container
#     output:
#         f"{dbs}/phanta/v1/phanta.db_ready",
#     shell:"""
#     wget https://github.com/bhattlab/phanta#advanced-usage
#     cd $(dirname {output})
#     wget --continue --read-timeout=1 -r wget https://www.dropbox.com/sh/3ktsdqlcph6x95r/AACGSj0sxYV6IeUQuGAFPtk8a/database_V1.tar.gz

#     tar xvzf database_V1.tar.gz
#     touch $(basename {output})
#     """

# rule get_nt:
#     container: "docker://ghcr.io/vdblab/kraken2-patch:2.1.2",
#     output:
#         f"{dbs}/kraken/zymogutv0/library/nt/library.fna"


# rule map_and_mask:

rule make_zymo_mock_profiling_db:
    # this doesn't use the normal kraken container image because of ftp issues
    # https://github.com/DerrickWood/kraken2/pull/637
    container: "docker://ghcr.io/vdblab/kraken2-patch:2.1.2",
    output:
        f"{dbs}/kraken/zymogutv0/hash.k2d"
    threads: 8
    resources:
        mem_mb=8*1024,
        runtime=24*60,
    shell: """
    DBNAME=$(dirname {output[0]})
    # we do this rm step because the snakemake outputs dont track every file in the db, just the trigger file
    rm -rf $DBNAME/library/ D6331.refseq/
#    kraken2-build --threads {threads} --download-library nt --db $DBNAME --use-ftp
#    kraken2-build --threads {threads} --download-library nt --db $DBNAME --use-ftp
#    kraken2-build --threads {threads} --download-library human --db $DBNAME --use-ftp
#    kraken2-build --threads {threads} --download-library viral --db $DBNAME --use-ftp
#    kraken2-build --threads {threads} --download-library archaea --db $DBNAME --use-ftp

    kraken2-build --download-taxonomy --db $DBNAME

    wget -N https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip
    # so if a build breaks, we want a clean dir to start over
    kraken2-build --download-library UniVec_Core --db $DBNAME

    unzip D6331.refseq.zip
    mkdir D6331.refseq/tmp/
    cat D6331.refseq/genomes/Akkermansia_muciniphila.fasta | sed "s,>,>kraken:taxid|239935,g" > D6331.refseq/tmp/Akkermansia_muciniphila.fasta
    cat D6331.refseq/genomes/Bacteroides_fragilis.fasta | sed "s,>,>kraken:taxid|817,g" > D6331.refseq/tmp/Bacteroides_fragilis.fasta
    cat D6331.refseq/genomes/Bifidobacterium_adolescentis.fasta | sed "s,>,>kraken:taxid|1680,g" > D6331.refseq/tmp/Bifidobacterium_adolescentis.fasta
    cat D6331.refseq/genomes/Candida_albican.fasta | sed "s,>,>kraken:taxid|5476y,g" > D6331.refseq/tmp/Candida_albican.fasta
    cat D6331.refseq/genomes/Clostridioides_difficile.fasta | sed "s,>,>kraken:taxid|1496,g" > D6331.refseq/tmp/Clostridioides_difficile.fasta
    cat D6331.refseq/genomes/Clostridium_perfringens.fasta | sed "s,>,>kraken:taxid|1502,g" > D6331.refseq/tmp/Clostridium_perfringens.fasta
    cat D6331.refseq/genomes/Enterococcus_faecalis.fasta | sed "s,>,>kraken:taxid|1351,g" > D6331.refseq/tmp/Enterococcus_faecalis.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B1109.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B1109.fasta
    cat D6331.refseq/genomes/Escherichia_coli_b2207.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_b2207.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B3008.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B3008.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B766.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B766.fasta
    cat D6331.refseq/genomes/Escherichia_coli_JM109.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_JM109.fasta
    cat D6331.refseq/genomes/Faecalibacterium_prausnitzii.fasta | sed "s,>,>kraken:taxid|853,g" > D6331.refseq/tmp/Faecalibacterium_prausnitzii.fasta
    cat D6331.refseq/genomes/Fusobacterium_nucleatum.fasta | sed "s,>,>kraken:taxid|851,g" > D6331.refseq/tmp/Fusobacterium_nucleatum.fasta
    cat D6331.refseq/genomes/Lactobacillus_fermentum.fasta | sed "s,>,>kraken:taxid|1613,g" > D6331.refseq/tmp/Lactobacillus_fermentum.fasta
    cat D6331.refseq/genomes/Methanobrevibacter_smithii.fasta | sed "s,>,>kraken:taxid|2173,g" > D6331.refseq/tmp/Methanobrevibacter_smithii.fasta
    cat D6331.refseq/genomes/Prevotella_corporis.fasta | sed "s,>,>kraken:taxid|28128,g" > D6331.refseq/tmp/Prevotella_corporis.fasta
    cat D6331.refseq/genomes/Roseburia_hominis.fasta | sed "s,>,>kraken:taxid|301301,g" > D6331.refseq/tmp/Roseburia_hominis.fasta
    cat D6331.refseq/genomes/Saccharomyces_cerevisiae.fasta | sed "s,>,>kraken:taxid|4932,g" > D6331.refseq/tmp/Saccharomyces_cerevisiae.fasta
    cat D6331.refseq/genomes/Salmonella_enterica.fasta | sed "s,>,>kraken:taxid|28901,g" > D6331.refseq/tmp/Salmonella_enterica.fasta
    cat D6331.refseq/genomes/Veillonella_rogosae.fasta | sed "s,>,>kraken:taxid|423477,g" > D6331.refseq/tmp/Veillonella_rogosae.fasta
    echo "Building DB"
    for file in D6331.refseq/tmp/*.fasta
    do
        echo "   - adding $file"
        kraken2-build --threads {threads} --add-to-library $file --db $DBNAME
    done
    rm -r D6331.refseq*
    kraken2-build --build --threads {threads} --db $DBNAME
    touch {output}
    """

rule make_zymo_mock_spike_profiling_db:
    # this doesn't use the normal kraken container image because of ftp issues
    # https://github.com/DerrickWood/kraken2/pull/637
    container: "docker://ghcr.io/vdblab/kraken2-patch:2.1.2",
    output:
        f"{dbs}/kraken/zymogutv1_plusspikes/hash.k2d"
    threads: 8
    resources:
        mem_mb=lambda wc, attempt: 8*1024 * attempt,
        runtime=lambda wc, attempt: 24*60* attempt,
    shell: """
    DBNAME=$(dirname {output[0]})
    # we do this rm step because the snakemake outputs dont track every file in the db, just the trigger file
    # so if a build breaks, we want a clean dir to start over
#    rm -rf $DBNAME/library/ D6331.refseq/
    kraken2-build --download-taxonomy --db $DBNAME

    wget -N https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip
    kraken2-build --download-library UniVec_Core --db $DBNAME
    kraken2-build --download-library bacteria --db $DBNAME
    kraken2-build --download-library archaea --db $DBNAME
    kraken2-build --download-library fungi --db $DBNAME
    kraken2-build --threads {threads} --download-library human --db $DBNAME

    unzip D6331.refseq.zip
    mkdir D6331.refseq/tmp/
    cat D6331.refseq/genomes/Akkermansia_muciniphila.fasta | sed "s,>,>kraken:taxid|239935,g" > D6331.refseq/tmp/Akkermansia_muciniphila.fasta
    cat D6331.refseq/genomes/Bacteroides_fragilis.fasta | sed "s,>,>kraken:taxid|817,g" > D6331.refseq/tmp/Bacteroides_fragilis.fasta
    cat D6331.refseq/genomes/Bifidobacterium_adolescentis.fasta | sed "s,>,>kraken:taxid|1680,g" > D6331.refseq/tmp/Bifidobacterium_adolescentis.fasta
    cat D6331.refseq/genomes/Candida_albican.fasta | sed "s,>,>kraken:taxid|5476y,g" > D6331.refseq/tmp/Candida_albican.fasta
    cat D6331.refseq/genomes/Clostridioides_difficile.fasta | sed "s,>,>kraken:taxid|1496,g" > D6331.refseq/tmp/Clostridioides_difficile.fasta
    cat D6331.refseq/genomes/Clostridium_perfringens.fasta | sed "s,>,>kraken:taxid|1502,g" > D6331.refseq/tmp/Clostridium_perfringens.fasta
    cat D6331.refseq/genomes/Enterococcus_faecalis.fasta | sed "s,>,>kraken:taxid|1351,g" > D6331.refseq/tmp/Enterococcus_faecalis.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B1109.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B1109.fasta
    cat D6331.refseq/genomes/Escherichia_coli_b2207.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_b2207.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B3008.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B3008.fasta
    cat D6331.refseq/genomes/Escherichia_coli_B766.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_B766.fasta
    cat D6331.refseq/genomes/Escherichia_coli_JM109.fasta | sed "s,>,>kraken:taxid|562,g" > D6331.refseq/tmp/Escherichia_coli_JM109.fasta
    cat D6331.refseq/genomes/Faecalibacterium_prausnitzii.fasta | sed "s,>,>kraken:taxid|853,g" > D6331.refseq/tmp/Faecalibacterium_prausnitzii.fasta
    cat D6331.refseq/genomes/Fusobacterium_nucleatum.fasta | sed "s,>,>kraken:taxid|851,g" > D6331.refseq/tmp/Fusobacterium_nucleatum.fasta
    cat D6331.refseq/genomes/Lactobacillus_fermentum.fasta | sed "s,>,>kraken:taxid|1613,g" > D6331.refseq/tmp/Lactobacillus_fermentum.fasta
    cat D6331.refseq/genomes/Methanobrevibacter_smithii.fasta | sed "s,>,>kraken:taxid|2173,g" > D6331.refseq/tmp/Methanobrevibacter_smithii.fasta
    cat D6331.refseq/genomes/Prevotella_corporis.fasta | sed "s,>,>kraken:taxid|28128,g" > D6331.refseq/tmp/Prevotella_corporis.fasta
    cat D6331.refseq/genomes/Roseburia_hominis.fasta | sed "s,>,>kraken:taxid|301301,g" > D6331.refseq/tmp/Roseburia_hominis.fasta
    cat D6331.refseq/genomes/Saccharomyces_cerevisiae.fasta | sed "s,>,>kraken:taxid|4932,g" > D6331.refseq/tmp/Saccharomyces_cerevisiae.fasta
    cat D6331.refseq/genomes/Salmonella_enterica.fasta | sed "s,>,>kraken:taxid|28901,g" > D6331.refseq/tmp/Salmonella_enterica.fasta
    cat D6331.refseq/genomes/Veillonella_rogosae.fasta | sed "s,>,>kraken:taxid|423477,g" > D6331.refseq/tmp/Veillonella_rogosae.fasta

    # Salinibacter ruber DSM 13855,
    # make kraken-friendly headers.
    mkdir -p tmp
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.fna.gz | gunzip | sed "s,>,>kraken:taxid|309807,g" > D6331.refseq/tmp/Salinibacter_ruber.fasta
    # getting Trichoderma reesei  QM6a instead of  Trichoderma reesei ATCC 13631 .  Might be fine license-wise but not sure
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip | sed "s,>,>kraken:taxid|51453,g" > D6331.refseq/tmp/Trichoderma_reesei.fasta
    # getting CBA1121 instead of Haloarcula hispanica ATCC 33960  for the same reason.
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/729/095/GCF_008729095.1_ASM872909v1/GCF_008729095.1_ASM872909v1_genomic.fna.gz | gunzip | sed "s,>,>kraken:taxid|634497,g" > D6331.refseq/tmp/Haloarcula_hispanica.fasta


    echo "Building DB"
    for file in D6331.refseq/tmp/*.fasta
    do
        echo "   - adding $file"
        kraken2-build --threads {threads} --add-to-library $file --db $DBNAME
    done
    rm -r D6331.refseq*
    rm -r tmp/
    kraken2-build --build --threads {threads} --db $DBNAME
    touch {output}
    """

rule make_bracken_dbs_for_zymo_mock:
    container: "docker://ghcr.io/vdblab/kraken2:2.1.2",
    input:
        f"{dbs}/kraken/{{mock}}/hash.k2d"
    output:
        f"{dbs}/kraken/{{mock}}/database{{readlen}}mers.kmer_distrib",
    resources:
        mem_mb=lambda wc, attempt: 150 *1024 * attempt,
        runtime=lambda wc, attempt: 24*60* attempt,
    threads: 32
    shell: """
    DB=$(dirname {output})
    # see
    # bracken-build -d $DB -t {threads} -k 35 -l {wildcards.readlen}
    # 1a
    kraken2 --db=$DB --threads={threads} <( find -L $DB/library \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) -exec cat {{}} + ) > database.kraken
    # 1b
    kmer2read_distr --seqid2taxid $DB/seqid2taxid.map --taxonomy $DB/taxonomy --kraken database.kraken --output database{wildcards.readlen}mers.kraken -k 35 -l {wildcards.readlen} -t {threads}
    # 1c
    generate_kmer_distribution.py -i database{wildcards.readlen}mers.kraken -o $DB/database{wildcards.readlen}mers.kmer_distrib
    """

rule complete_bracken_db:
    """this rule exists to ensure that all the kmer distribution files are
    present before "oking" the db
    """
    container: "docker://ghcr.io/vdblab/kraken2:2.1.2",
    input:
        expand(f"{dbs}/kraken/{{mock}}/database{{readlen}}mers.kmer_distrib",
               readlen=bracken_readlens,
               mock=["zymogutv0", "zymogutv1_plusspikes"]
               ),

    output:
        f"{dbs}/kraken/zymogutv1/mock.db_ready",
    threads: 1
    shell: """
    rm database.kraken
    rm database*mers.kraken
    touch {output}
    """
