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
DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)
LONG_DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=LONG_READ_LENGTHS)
ENDS = ("pe", "se")

VARIATION_SETTINGS = {
    "sim1": "",
    "sim3": "--snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50",
    "sim4": "--snp-rate 0.005 --small-indel-rate 0.0005 --max-small-indel-size 50",
    "sim5": "--snp-rate 0.005 --small-indel-rate 0.001 --max-small-indel-size 100",
    "sim6": "--snp-rate 0.05 --small-indel-rate 0.002 --max-small-indel-size 100",
}
SIM = ["sim0", "sim1"] + list(VARIATION_SETTINGS)


wildcard_constraints:
    read_length=r"\d{2,3}",
    long_read_length=r"\d{4,5}",
    sim01=r"sim(0|0p1)"


localrules:
    download_fruitfly, download_maize, download_chm13, download_rye, download_ecoli50, filter_ecoli50, filter_fruitfly, clone_seqan, samtools_faidx


rule:
    input:
        expand("datasets/{sim}/{ds}/{r}.fastq.gz", sim=SIM, ds=DATASETS, r=(1, 2)),
        expand("datasets/{sim}/{ds}/truth.bam", sim=SIM, ds=DATASETS + LONG_DATASETS),
        expand("datasets/{sim}/{ds}/1.fastq.gz", sim=SIM, ds=LONG_DATASETS),

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
        vcf="variants/{sim}-{genome}.vcf"
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


rule sim01:
    output:
        r1_fastq="datasets/{sim01}/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{sim01}/{genome}-{read_length}/2.fastq.gz",
        bam="datasets/{sim01}/{genome}-{read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        extra=readsimulator_parameters,
        n_reads=lambda wildcards: N_READS[int(wildcards.read_length)],
        error_rate=lambda wildcards: {"sim0": 0.0, "sim0p1": 0.1}[wildcards.sim01]
    shell:
        "python readsimulator.py{params.extra} -e {params.error_rate} -n {params.n_reads} --read-length {wildcards.read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -1 {output.r1_fastq} -2 {output.r2_fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


rule sim01_long:
    output:
        fastq="datasets/{sim01}/{genome}-{long_read_length}/1.fastq.gz",
        bam="datasets/{sim01}/{genome}-{long_read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        n_reads=lambda wildcards: N_READS[int(wildcards.long_read_length)],
        error_rate=lambda wildcards: {"sim0": 0.0, "sim0p1": 0.1}[wildcards.sim01]
    shell:
        "python readsimulator.py --se -e {params.error_rate} -n {params.n_reads} --read-length {wildcards.long_read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -0 {output.fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


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
