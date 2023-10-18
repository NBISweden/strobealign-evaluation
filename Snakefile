N_READS = 1_000_000

# GENOMES = ("drosophila", "maize", "CHM13", "rye")
# READ_LENGTHS = (50, 75, 100, 150, 200, 300, 500)
GENOMES = ("drosophila", "maize", "CHM13")
READ_LENGTHS = (50, 100, 200)

DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)

COMMITS = {
    "min": "3223dc5946d9f38814e25a25149548dd146cc8d0",  # original
    "max": "6a837431f29fc3be3c3a74bea507538d9ea5abe7",
}

rule:
    input: "table.tex"


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
        "samtools view -f 64 --remove-flags 239 -o {output.bam} {input.bam}"


# Map reads with BWA-MEM

rule bwa_index:
    output:
        "{genome}.fa.amb",
        "{genome}.fa.ann",
        "{genome}.fa.bwt",
        "{genome}.fa.pac",
        "{genome}.fa.sa"
    input: "{genome}.fa"
    shell:
        "bwa index {input}"

rule map_bwa_mem_single_end:
    output:
        bam="bwamem/{genome}-{read_length}/se.bam"
    input:
        fasta="genomes/{genome}.fa",
        index="genomes/{genome}.fa.bwt",
        r1_fastq="datasets/{genome}-{read_length}/1.fastq.gz",
    threads: 20
    log:
        "bwamem/{genome}-{read_length}/se.bam.log"
    shell:
        "/usr/bin/time -v bwa mem -t {threads} {input.fasta} {input.r1_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"

rule map_bwa_mem_paired_end:
    output:
        bam="bwamem/{genome}-{read_length}/pe.bam"
    input:
        fasta="genomes/{genome}.fa",
        index="genomes/{genome}.fa.bwt",
        r1_fastq="datasets/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{genome}-{read_length}/2.fastq.gz",
    threads: 20
    log:
        "bwamem/{genome}-{read_length}/pe.bam.log"
    shell:
        "/usr/bin/time -v bwa mem -t {threads} {input.fasta} {input.r1_fastq} {input.r2_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"


# Compile and run the two versions of strobealign

rule compile_strobealign:
    output: "bin/strobealign-{name}"
    params:
        commit=lambda wildcards: COMMITS[wildcards.name]
    threads: 20
    shell:
        """
        wget https://github.com/ksahlin/strobealign/archive/{params.commit}.zip
        unzip {params.commit}.zip
        rm {params.commit}.zip
        cd strobealign-{params.commit}
        cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
        make -j {threads} -C build strobealign
        mv build/strobealign ../{output}
        cd ..
        rm -rf strobealign-{params.commit}
        """

rule run_strobealign_paired_end:
    output:
        bam="strobealign-{program,(min|max)}/{genome}-{read_length}/pe.bam"
    input:
        fasta="genomes/{genome}.fa",
        r1_fastq="datasets/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{genome}-{read_length}/2.fastq.gz",
    threads: 20
    log:
        "strobealign-{program}/{genome}-{read_length}/pe.bam.log"
    shell:
        "/usr/bin/time -v bin/strobealign-{wildcards.program} -t {threads} {input.fasta} {input.r1_fastq} {input.r2_fastq} 2> {log}.tmp"
        " | grep -v '^@PG'"
        " | samtools view --no-PG -o {output.bam}.tmp.bam -"
        "\n mv -v {output.bam}.tmp.bam {output.bam}"
        "\n mv -v {log}.tmp {log}"

rule combine_strobealign_min_max:
    output:
        bam="strobealign-combined/{dataset}/{ends}.bam"
    input:
        bam1="strobealign-min/{dataset}/{ends}.bam",
        bam2="strobealign-max/{dataset}/{ends}.bam"
    params:
        extra=lambda wildcards: "-s" if wildcards.ends == "se" else ""
    shell:
        "python combine.py {params.extra} --output {output.bam} {input.bam1} {input.bam2}"

# Accuracy for all programs

rule get_accuracy:
    output:
        txt="{program}/{genome}-{read_length}/accuracy.{ends}.txt"
    input:
        truth="datasets/{genome}-{read_length}/truth.{ends}.bam",
        bam="{program}/{genome}-{read_length}/{ends}.bam"
    shell:
        "python get_accuracy.py --truth {input.truth} --predicted_sam {input.bam} > {output.txt}.tmp"
        "\n mv -v {output.txt}.tmp {output.txt}"




def read_accuracy(path) -> float:
    with open(path) as f:
        line = next(iter(f))
    fields = line.split()
    return float(fields[1])


# TODO se!

rule accuracy_table:
    output:
        tex="table.tex"
    input:
        expand(
            "{program}/{dataset}/accuracy.{ends}.txt",
            program=("bwamem", "strobealign-min", "strobealign-max", "strobealign-combined"),
            dataset=DATASETS,
            ends=("pe", )
        )
    run:
        with open(output.tex, "w") as f:
            for ends in ("pe", ): #"se"):
                title = "Single-end" if ends == "se" else "Paired-end"
                print(f"\n# {title}", file=f)
                print(r"\begin{tabular}{lrrrr}", file=f)
                print(r"dataset &     min &     max & combined& combined minus min & BWA minus combined\\", file=f)
                for dataset in DATASETS:
                    accuracies = {
                        program: read_accuracy(f"{program}/{dataset}/accuracy.{ends}.txt")
                        for program in ("bwamem", "strobealign-min", "strobealign-max", "strobealign-combined")
                    }
                    delta = accuracies["strobealign-combined"] - accuracies["strobealign-min"]
                    bwa_delta = accuracies["bwamem"] - accuracies["strobealign-combined"]
                    print(
                        f"{dataset:>14s} & "
                        f"{accuracies['strobealign-min']:.4f} & "
                        f"{accuracies['strobealign-max']:.4f} & "
                        f"{accuracies['strobealign-combined']:.4f} & "
                        f"{delta:+.4f} & "
                        f"{bwa_delta:+.4f}\\\\",
                        file=f
                    )
                print("\\end{tabular}", file=f)


# Misc

rule samtools_faidx:
    output: "{genome}.fa.fai"
    input: "{genome}.fa"
    shell: "samtools faidx {input}"
