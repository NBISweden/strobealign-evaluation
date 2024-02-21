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


rule:
    input:
        expand("datasets/{ds}/{r}.fastq.gz", ds=DATASETS, r=(1, 2)),
        expand("datasets/{ds}/truth.{ends}.bam", ds=DATASETS, ends=("se", "pe"))

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
        vcf="variations/{genome}.vcf"
    input:
        fasta="genomes/{genome}.fa",
        fai="genomes/{genome}.fa.fai"
    shell:
        "mason_variator -ir {input.fasta} --sv-indel-rate 0.000005 --snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50 -ov {output.vcf}.tmp.vcf"
        "\n mv -v {output.vcf}.tmp.vcf {output.vcf}"

rule mason_simulator:
    output:
        r1_fastq="datasets/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{genome}-{read_length}/2.fastq.gz",
        sam="datasets/{genome}-{read_length}/truth.pe.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf=rules.mason_variator.output.vcf
    run:
        extra = ""
        if wildcards.read_length in {"250", "300", "500"}:
            extra="--fragment-mean-size 700"
        shell(
            "mason_simulator -ir {input.fasta} -n {N_READS} -iv {input.vcf} --illumina-read-length {wildcards.read_length} -o {output.r1_fastq}.tmp.fastq.gz -or {output.r2_fastq}.tmp.fastq.gz -oa {output.sam}.tmp.bam {extra}"
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
