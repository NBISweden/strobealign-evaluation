# Snakefile for generating the input datasets.
# It downloads reference FASTAs and simulates reads from them.
#
# This creates the downloads/, datasets/ and genomes/ directories.
# (downloads/ can be deleted.)

N_READS = 1_000_000

GENOMES = ("drosophila", "maize", "CHM13", "rye", "ecoli50")
READ_LENGTHS = (50, 75, 100, 150, 200, 300, 500)

DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)
ENDS = ("pe", "se")

VARIATION_SETTINGS = {
    "sim3": "--sv-indel-rate 0.000005 --snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50",
    "sim5": "--sv-indel-rate 0.00002 --snp-rate 0.005 --small-indel-rate 0.001 --max-small-indel-size 100",
}
SIM = list(VARIATION_SETTINGS)


localrules:
    download_drosophila, download_maize, download_chm13, download_rye, download_ecoli50, filter_ecoli50, filter_drosophila, clone_seqan

rule:
    input:
        expand("datasets/{sim}/{ds}/{r}.fastq.gz", sim=SIM, ds=DATASETS, r=(1, 2)),
        expand("datasets/{sim}/{ds}/truth.{ends}.bam", sim=SIM, ds=DATASETS, ends=ENDS)

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


# Generate simulated reads and BAM files with expected alignments (truth)

rule mason_variator:
    output:
        vcf="variations/{sim}-{genome}.vcf"
    input:
        fasta="genomes/{genome}.fa",
        fai="genomes/{genome}.fa.fai",
        mason_variator="bin/mason_variator"
    params:
        variation_settings=lambda wildcards: VARIATION_SETTINGS[wildcards.sim]
    shell:
        "{input.mason_variator} -ir {input.fasta} {params.variation_settings} -ov {output.vcf}.tmp.vcf"
        "\n mv -v {output.vcf}.tmp.vcf {output.vcf}"

rule mason_simulator:
    output:
        r1_fastq="datasets/{sim}/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{sim}/{genome}-{read_length}/2.fastq.gz",
        sam="datasets/{sim}/{genome}-{read_length}/truth.pe.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf=rules.mason_variator.output.vcf,
        mason_simulator="bin/mason_simulator"
    run:
        extra = ""
        if wildcards.read_length in {"250", "300", "500"}:
            extra="--fragment-mean-size 700"
        shell(
            "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
            "\n {input.mason_simulator} -ir {input.fasta} -n {N_READS} -iv {input.vcf} --illumina-read-length {wildcards.read_length} -o {output.r1_fastq}.tmp.fastq.gz -or {output.r2_fastq}.tmp.fastq.gz -oa {output.sam}.tmp.bam {extra}"
            "\n mv -v {output.r1_fastq}.tmp.fastq.gz {output.r1_fastq}"
            "\n mv -v {output.r2_fastq}.tmp.fastq.gz {output.r2_fastq}"
            "\n mv -v {output.sam}.tmp.bam {output.sam}"
        )

rule single_end_truth:
    output:
        bam="datasets/{dataset}/truth.se.bam"
    input:
        bam="datasets/{dataset}/truth.pe.bam"
    shell:
        "samtools view -f 64 --remove-flags 235 -o {output.bam} {input.bam}"


# Misc

rule samtools_faidx:
    output: "{genome}.fa.fai"
    input: "{genome}.fa"
    shell: "samtools faidx {input}"


# Build our own Mason binaries because the one from Conda crashes
rule clone_seqan:
    output: "seqan/cloned"
    shell:
        "git clone --branch db5e0ce7e0b7946ff5d1ca22e652faa0b5b9603c https://github.com/seqan/seqan.git"
        "; touch seqan/cloned"

rule build_mason:
    output: "bin/mason_variator", "bin/mason_simulator"
    input: "seqan/cloned"
    threads: 99
    shell:
        "cmake -DSEQAN_BUILD_SYSTEM=APP:mason2 -B build-seqan seqan"
        "; make -s -C build-seqan -j {threads}"
        "; mv build-seqan/bin/mason_simulator build-seqan/bin/mason_variator bin/"
        #"; rm -r seqan"
