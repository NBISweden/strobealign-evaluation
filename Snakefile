# Snakefile for generating the input datasets.
# It downloads reference FASTAs and simulates reads from them.
#
# This creates the downloads/, datasets/ and genomes/ directories.
# (downloads/ can be deleted.)

N_READS = 1_000_000

# Additional 'genomes' that can be added to this list:
# - ecoli50 (fifty E. coli genomes)
# - chrY (chromosome Y of CHM13)

GENOMES = ("drosophila", "maize", "CHM13", "rye")
READ_LENGTHS = (50, 75, 100, 150, 200, 300, 500)

DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)
ENDS = ("pe", "se")

VARIATION_SETTINGS = {
    "sim3": "--sv-indel-rate 0.000005 --snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50",
    "sim4": "--sv-indel-rate 0.00001 --snp-rate 0.005 --small-indel-rate 0.0005 --max-small-indel-size 50",
    "sim5": "--sv-indel-rate 0.00002 --snp-rate 0.005 --small-indel-rate 0.001 --max-small-indel-size 100",

    # SIM6 uses a different fruit fly genome (fruitfly.fa instead of drosophila.fa), which does not contain
    # contigs smaller than 10 kbp. This is to avoid crashes in mason_simulator.
    "sim6": "--sv-indel-rate 0.001 --snp-rate 0.05 --small-indel-rate 0.002 --max-small-indel-size 100",
}
SIM = ["sim0"] + list(VARIATION_SETTINGS)


localrules:
    download_drosophila, download_maize, download_chm13, download_rye, download_ecoli50, filter_ecoli50, filter_drosophila, clone_seqan

rule:
    input:
        expand("datasets/{sim}/{ds}/{r}.fastq.gz", sim=SIM, ds=DATASETS, r=(1, 2)),
        expand("datasets/{sim}/{ds}/truth.bam", sim=SIM, ds=DATASETS, ends=ENDS)

# Download genomes

rule download_drosophila:
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

rule filter_drosophila:
    output: "genomes/drosophila.fa"
    input: rules.download_drosophila.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 1000 bp
        awk '$2>=1000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
        samtools faidx -r {output}.tmp.regions.txt {output}.tmp.fa > {output}.tmp2.fa
        mv {output}.tmp2.fa {output}

        rm {output}.tmp.fa {output}.tmp.fa.fai {output}.tmp.regions.txt
        """

# A variant of drosophila where all contigs smaller than 10 kbp are removed
# (instead of 1 kbp). This is used for SIM6.
rule filter_fruitfly:
    output: "genomes/fruitfly.fa"
    input: rules.download_drosophila.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 1000 bp
        awk '$2>=10000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
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
    output: "genomes/rye.fa"
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
        "{input.mason_variator} -ir {input.fasta} {params.variation_settings} -ov {output.vcf}.tmp.vcf"
        "\n mv -v {output.vcf}.tmp.vcf {output.vcf}"


def mason_simulator_parameters(wildcards):
    read_length = int(wildcards.read_length)
    if read_length == 500 and wildcards.sim == "sim5" and wildcards.genome == "drosophila":
        # Workaround for crash with length 500
        result = "--illumina-read-length 460"
    else:
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
        fasta=lambda wildcards: "genomes/fruitfly.fa" if wildcards.genome == "drosophila" and wildcards.sim == "sim6" else "genomes/{genome}.fa".format(genome=wildcards.genome),
        vcf=lambda wildcards: "variants/sim6-fruitfly.vcf" if wildcards.genome == "drosophila" and wildcards.sim == "sim6" else "variants/{sim}-{genome}.vcf".format(sim=wildcards.sim, genome=wildcards.genome),
        mason_simulator="bin/mason_simulator"
    params:
        extra=mason_simulator_parameters
    log: "logs/mason_simulator/{sim}-{genome}-{read_length}.log"
    shell:
        "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
        "\n{input.mason_simulator}"
        " --num-threads 1"  # Output depends on number of threads, leave at 1 for reproducibility
        " -ir {input.fasta}"
        " -n {N_READS}"
        " -iv {input.vcf}"
        " {params.extra}"
        " -o {output.r1_fastq}.tmp.fastq.gz"
        " -or {output.r2_fastq}.tmp.fastq.gz"
        " -oa {output.bam}.tmp.bam"
        " 2>&1 | tee {log}"
        "\nmv -v {output.r1_fastq}.tmp.fastq.gz {output.r1_fastq}"
        "\nmv -v {output.r2_fastq}.tmp.fastq.gz {output.r2_fastq}"
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
        extra=readsimulator_parameters
    shell:
        "python readsimulator.py{params.extra} -n {N_READS} --read-length {wildcards.read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -1 {output.r1_fastq} -2 {output.r2_fastq} {output.bam}.tmp.bam"
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
