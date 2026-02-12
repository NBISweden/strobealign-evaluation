# Snakefile for generating the input datasets.
# It downloads reference FASTAs and simulates reads from them.
#
# This creates the downloads/, datasets/ and genomes/ directories.
# (downloads/ can be deleted.)

# Additional 'genomes' that can be added to this list:
# - ecoli50 (fifty E. coli genomes)
# - chrY (chromosome Y of CHM13)

GENOMES = ("fruitfly", "maize", "CHM13", "rye", "chrY")
N_READS = {
    50: 1_000_000,
    75: 1_000_000,
    100: 1_000_000,
    150: 1_000_000,
    200: 1_000_000,
    300: 1_000_000,
    500: 1_000_000,
    1000: 500_000,
    5000: 100_000,
    10000: 50_000,
}
LONG_READ_LENGTHS = tuple(n for n in N_READS if n >= 1000)  # single-end only
READ_LENGTHS = tuple(n for n in N_READS if n < 1000)
MODELS = {"clr": "data/pbsim3/QSHMM-RSII.model", "ont": "data/pbsim3/QSHMM-ONT-HQ.model", "hifi": "data/pbsim3/QSHMM-RSII.model"}
DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)
LONG_DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=LONG_READ_LENGTHS)
ENDS = ("pe", "se")

# sim3/4/5 were used in earlier studies, but to get a wider spread of
# error/variation rates, we switched to sim0, sim4, sim6.
VARIATION_SETTINGS = {
    "sim1": "",
#    "sim3": "--snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50",
    "sim4": "--snp-rate 0.005 --small-indel-rate 0.0005 --max-small-indel-size 50",
#    "sim5": "--snp-rate 0.005 --small-indel-rate 0.001 --max-small-indel-size 100",
    "sim6": "--snp-rate 0.05 --small-indel-rate 0.002 --max-small-indel-size 100",
}
SIM = ["sim0"] + list(VARIATION_SETTINGS)
LONG_SIM = ["ont", "hifi", "clr"]


wildcard_constraints:
    read_length=r"\d{2,3}",
    long_read_length=r"\d{4,5}",


localrules:
    download_fruitfly, download_maize, download_chm13, download_rye, download_ecoli50, filter_ecoli50, filter_fruitfly, clone_seqan, samtools_faidx


rule:
    input:
        expand("datasets/{sim}/{ds}/{r}.fastq.gz", sim=SIM, ds=DATASETS, r=(1, 2)),
        expand("datasets/{sim}/{ds}/truth.bam", sim=SIM, ds=DATASETS + LONG_DATASETS),
        expand("datasets/{sim}/{ds}/1.fastq.gz", sim=SIM + LONG_SIM, ds=LONG_DATASETS),
        expand("datasets/{sim}/{ds}/truth.maf.gz", sim=LONG_SIM, ds=LONG_DATASETS)

# Download genomes

rule download_fruitfly:
    output: "downloads/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz"
    shell:
        "curl ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz > {output}"

rule download_maize:
    output: "downloads/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
    shell:
        "curl https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz > {output}"

rule download_chm13:
    output: "downloads/chm13v2.0.fa.gz"
    shell:
        "curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz > {output}"

rule download_rye:
    output: "downloads/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz"
    shell:
        "curl https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Secale_cereale/latest_assembly_versions/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz > {output}"

rule download_ecoli50:
    output:
        "downloads/ecoli50/done.txt"
    input: "ecoli-accessions.txt"
    shell:
        "mkdir -p downloads/ecoli50; "
        "head -n 50 ecoli-accessions.txt > downloads/ecoli50/accessions.txt; "
        "ncbi-genome-download -o downloads/ecoli50 -A downloads/ecoli50/accessions.txt --formats fasta --flat-output bacteria; "
        "touch downloads/ecoli50/done.txt"

rule filter_ecoli50:
    output: "genomes/ecoli50.fa"
    input: "downloads/ecoli50/done.txt"
    shell:
        "python noplasmids.py downloads/ecoli50/*.fna.gz > {output}"


# Uncompress and filter downloaded genomes

rule filter_fruitfly:
    output: "genomes/fruitfly.fa"
    input: rules.download_fruitfly.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 16 kbp
        awk '$2>=16000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
        samtools faidx -r {output}.tmp.regions.txt {output}.tmp.fa > {output}.tmp2.fa
        mv {output}.tmp2.fa {output}

        rm {output}.tmp.fa {output}.tmp.fa.fai {output}.tmp.regions.txt
        """

rule uncompress_maize:
    output: "genomes/maize.fa"
    input: rules.download_maize.output
    shell:
        "zcat {input} > {output}"

rule uncompress_chm13:
    output: "genomes/CHM13.fa"
    input: rules.download_chm13.output
    shell:
        "zcat {input} > {output}"

rule filter_rye:
    output: temp("genomes/rye_raw.fa")
    input: rules.download_rye.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 50000
        awk '$2>=50000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
        samtools faidx -r {output}.tmp.regions.txt {output}.tmp.fa > {output}.tmp2.fa
        mv {output}.tmp2.fa {output}

        rm {output}.tmp.fa {output}.tmp.fa.fai {output}.tmp.regions.txt
        """

rule chunk_rye:
    input:
        fasta="genomes/rye_raw.fa"
    output:
        fasta="genomes/rye.fa"
    shell:
        "python chunker.py {input.fasta} > {output.fasta}"


rule extract_chry:
    output: "genomes/chrY.fa"
    input: fasta="genomes/CHM13.fa", fai="genomes/CHM13.fa.fai"
    shell:
        "samtools faidx {input.fasta} chrY > {output}"


# Generate simulated reads and BAM files with expected alignments (truth)

rule mason_variator:
    output:
        vcf="variants/{sim,sim[1-9]}-{genome}.vcf"
    input:
        fasta="genomes/{genome}.fa",
        fai="genomes/{genome}.fa.fai",
        mason_variator="bin/mason_variator"
    params:
        variation_settings=lambda wildcards: VARIATION_SETTINGS[wildcards.sim]
    shell:
        """
        if [ "{wildcards.sim}" = "sim1" ]; then
            touch {output.vcf}
        else
            {input.mason_variator} -ir {input.fasta} {params.variation_settings} -ov {output.vcf}.tmp.vcf
            mv -v {output.vcf}.tmp.vcf {output.vcf}
        fi
        """


def mason_simulator_parameters(wildcards):
    read_length = int(wildcards.read_length)
    result = f"--illumina-read-length {read_length}"
    if read_length >= 250:
        result += " --fragment-mean-size 700"
    return result


rule mason_simulator:
    output:
        r1_fastq="datasets/{sim,sim[1-9]}/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{sim,sim[1-9]}/{genome}-{read_length}/2.fastq.gz",
        bam="datasets/{sim,sim[1-9]}/{genome}-{read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf="variants/{sim}-{genome}.vcf",
        mason_simulator="bin/mason_simulator"
    params:
        extra=mason_simulator_parameters,
        n_reads=lambda wildcards: N_READS[int(wildcards.read_length)],
        vcf_arg=lambda wildcards: "" if wildcards.sim == "sim1" else f"-iv variants/{wildcards.sim}-{wildcards.genome}.vcf"
    log: "logs/mason_simulator/{sim}-{genome}-{read_length}.log"
    shell:
        "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
        "\n{input.mason_simulator}"
        " --num-threads 1"  # Output depends on number of threads, leave at 1 for reproducibility
        " -ir {input.fasta}"
        " -n {params.n_reads}"
        " {params.vcf_arg}"
        " {params.extra}"
        " -o {output.r1_fastq}.tmp.fastq.gz"
        " -or {output.r2_fastq}.tmp.fastq.gz"
        " -oa {output.bam}.tmp.bam"
        " 2>&1 | tee {log}"
        "\nmv -v {output.r1_fastq}.tmp.fastq.gz {output.r1_fastq}"
        "\nmv -v {output.r2_fastq}.tmp.fastq.gz {output.r2_fastq}"
        "\nmv -v {output.bam}.tmp.bam {output.bam}"


rule mason_simulator_long:
    output:
        fastq="datasets/{sim,sim[1-9]}/{genome}-{long_read_length}/1.fastq.gz",
        bam="datasets/{sim,sim[1-9]}/{genome}-{long_read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf="variants/{sim}-{genome}.vcf",
        mason_simulator="bin/mason_simulator"
    params:
        n_reads=lambda wildcards: N_READS[int(wildcards.long_read_length)],
        fragment_length=lambda wildcards: int(int(wildcards.long_read_length) * 1.5),
        vcf_arg=lambda wildcards: "" if wildcards.sim == "sim1" else f"-iv variants/{wildcards.sim}-{wildcards.genome}.vcf"
    log: "logs/mason_simulator/{sim}-{genome}-{long_read_length}.log"
    shell:
        "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
        "\n{input.mason_simulator}"
        " --num-threads 1"  # Output depends on number of threads, leave at 1 for reproducibility
        " --illumina-read-length {wildcards.long_read_length}"
        " --fragment-mean-size {params.fragment_length}"
        " -ir {input.fasta}"
        " -n {params.n_reads}"
        " {params.vcf_arg}"
        " -o {output.fastq}.tmp.fastq.gz"
        " -oa {output.bam}.tmp.bam"
        " 2>&1 | tee {log}"
        "\nmv -v {output.fastq}.tmp.fastq.gz {output.fastq}"
        "\nmv -v {output.bam}.tmp.bam {output.bam}"


def readsimulator_parameters(wildcards):
    read_length = int(wildcards.read_length)
    if read_length >= 250:
        return " --mean-insert-size 700"
    return ""


rule sim0:
    output:
        r1_fastq="datasets/sim0/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/sim0/{genome}-{read_length}/2.fastq.gz",
        bam="datasets/sim0/{genome}-{read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        extra=readsimulator_parameters,
        n_reads=lambda wildcards: N_READS[int(wildcards.read_length)],
    shell:
        "python readsimulator.py{params.extra} -n {params.n_reads} --read-length {wildcards.read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -1 {output.r1_fastq} -2 {output.r2_fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


rule sim0_long:
    output:
        fastq="datasets/sim0/{genome}-{long_read_length}/1.fastq.gz",
        bam="datasets/sim0/{genome}-{long_read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        n_reads=lambda wildcards: N_READS[int(wildcards.long_read_length)],
    shell:
        "python readsimulator.py --se -n {params.n_reads} --read-length {wildcards.long_read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -0 {output.fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


def pbsim_parameters(wildcards):
    mean_read_length = int(wildcards.long_read_length)
    result = f"--length-mean {mean_read_length}"

    fai_path = "genomes/" + wildcards.genome + ".fa.fai"
    ref_len = 0
    with open(fai_path) as f:
        for line in f:
            fields = line.split()
            ref_len += int(fields[1])
    num_reads = N_READS[mean_read_length]
    depth = (num_reads * mean_read_length) / ref_len
    result += f" --depth {depth}"

    if wildcards.sim == "hifi":
        result += " --pass-num 10"
    return result


rule pbsim:
    output:
        maf="datasets/{sim,clr|ont}/{genome}-{long_read_length}/truth.maf.gz",
        fastq="datasets/{sim,clr|ont}/{genome}-{long_read_length}/1.fastq.gz"
    input:
        fasta="genomes/{genome}.fa",
        fai="genomes/{genome}.fa.fai",
        model=lambda wildcards: MODELS[wildcards.sim]
    params:
        extra=pbsim_parameters,
        outprefix="datasets/{sim}/{genome}-{long_read_length}/tmp",
        outid="S"
    log: "logs/pbsim3/{sim}-{genome}-{long_read_length}.log"
    shell:
        "pbsim"
        " --strategy wgs"
        " --genome {input.fasta}"
        " --method qshmm"
        " --qshmm {input.model}"
        " --prefix {params.outprefix}"
        " --id-prefix {params.outid}"
        " {params.extra}"
        " --length-sd 0"
        "\ncat {params.outprefix}_*.fq.gz > {output.fastq}"
        "\ncat {params.outprefix}_*.maf.gz > {output.maf}"
        "\nrm {params.outprefix}_*.ref"
        "\nrm {params.outprefix}_*.maf.gz"
        "\nrm {params.outprefix}_*.fq.gz"


rule pbsim_hifi:
    output:
        maf="datasets/{sim,hifi}/{genome}-{long_read_length}/truth.maf.gz",
        bam=temp("datasets/{sim,hifi}/{genome}-{long_read_length}/1.bam")
    input:
        fasta="genomes/{genome}.fa",
        model=lambda wildcards: MODELS[wildcards.sim]
    params:
        extra=pbsim_parameters,
        outprefix="datasets/{sim}/{genome}-{long_read_length}/tmp",
        outid="S"
    log: "logs/pbsim3/{sim,hifi}-{genome}-{long_read_length}.log"
    shell:
        "pbsim"
        " --strategy wgs"
        " --genome {input.fasta}"
        " --method qshmm"
        " --qshmm {input.model}"
        " --prefix {params.outprefix}"
        " --id-prefix {params.outid}"
        " {params.extra}"
        " --length-sd 0"
        "\ncat {params.outprefix}_*.maf.gz > {output.maf}"
        "\nrm {params.outprefix}_*.ref"
        "\nrm {params.outprefix}_*.maf.gz"
        "\nsamtools merge -o {output.bam} {params.outprefix}_*.bam"
        "\nrm {params.outprefix}_*.bam"


rule ccs:
    output:
        fastq="datasets/hifi/{genome}-{long_read_length}/1.fastq.gz"
    input:
        bam="datasets/hifi/{genome}-{long_read_length}/1.bam"
    log:
        "datasets/hifi/{genome}-{long_read_length}/ccs.log"
    threads:
        32
    shell:
        """
        ulimit -n 16384
        unset TMPDIR; ccs --log-file {log} -j {threads} {input.bam} {output.fastq} 
        """

# Misc

rule samtools_faidx:
    output: "{genome}.fa.fai"
    input: "{genome}.fa"
    shell: "samtools faidx {input}"


# Build our own Mason binaries because the one from Conda crashes
rule clone_seqan:
    output: "seqan/cloned"
    shell:
        "git clone https://github.com/seqan/seqan.git"
        "; ( cd seqan && git checkout seqan-v2.5.0rc2 )"
        "; touch seqan/cloned"


rule build_mason:
    output: "bin/mason_variator", "bin/mason_simulator"
    input: "seqan/cloned"
    threads: 99
    shell:
        "cmake -DSEQAN_BUILD_SYSTEM=APP:mason2 -DSEQAN_ARCH_SSE4=1 -B build-seqan seqan; "
        "cmake --build build-seqan -j {threads}; "
        "mv build-seqan/bin/mason_simulator build-seqan/bin/mason_variator bin/"
        #"; rm -r seqan"


rule download_pbsim_models:
    output: "data/pbsim3/QSHMM-ONT-HQ.model", "data/pbsim3/QSHMM-RSII.model"
    threads: 99
    shell:
        "mkdir -p data/pbsim3"
        "; cd data/pbsim3"
        "; wget -O- https://github.com/yukiteruono/pbsim3/archive/refs/tags/v3.0.5.tar.gz"
        " | tar xz --strip-components=2 pbsim3-3.0.5/data"
