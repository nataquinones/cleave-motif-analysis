import pandas as pd
import re
import os

def data_dir(path):
    return os.path.join('data', path)
df = pd.read_csv(data_dir('sample_metadata.tsv'), sep='\t')
samples = df['sample']
read_dict = dict(zip(df['sample'], df['read_file']))
references = ['ecobl21']

# -----------------------------------------------------------------------------

rule all:
    input:
        expand(data_dir('aligned/to_{reference}/{sample}.red.bed'), sample=samples, reference=references),
        expand(data_dir('aligned/to_{reference}/{sample}.cov'), sample=samples, reference=references)


# -----------------------------------------------------------------------------
# remove adaptors with cutadapt
rule trim_reads_with_cutadapt:
    input:
        reads_R1 = lambda wc: data_dir(f'reads/{read_dict[wc.sample]}_R1_001.fastq'),
        reads_R2 = lambda wc: data_dir(f'reads/{read_dict[wc.sample]}_R2_001.fastq')
    output:
        cuta_trimmed_reads_R1 = data_dir('reads/trimmed/{sample}_R1.cutadapt.trim.fastq'),
        cuta_trimmed_reads_R2 = data_dir('reads/trimmed/{sample}_R2.cutadapt.trim.fastq')
    conda:
        'envs/mshg.yml'
    resources:
        mem_mb=6000
    shell:
        'cutadapt '\
            '-a CTGTCTCTTATACACATCT '\
            '-A CTGTCTCTTATACACATCT '\
            '--revcomp '\
            '--minimum-length 1 '\
            '--poly-a '\
            '-o {output.cuta_trimmed_reads_R1} '\
            '-p {output.cuta_trimmed_reads_R2} '\
            '{input.reads_R1} '\
            '{input.reads_R2}'


# trim reads for quality with fastp
rule trim_reads_with_fastp:
    input:
        cuta_trimmed_reads_R1 = data_dir('reads/trimmed/{sample}_R1.cutadapt.trim.fastq'),
        cuta_trimmed_reads_R2 = data_dir('reads/trimmed/{sample}_R2.cutadapt.trim.fastq')
    output:
        fastp_trimmed_reads_R1 = data_dir('reads/trimmed/{sample}_R1.fastp.trim.fastq'),
        fastp_trimmed_reads_R2 = data_dir('reads/trimmed/{sample}_R2.fastp.trim.fastq'),
        report = data_dir('reads/fastp/{sample}.html'),
        json = data_dir('reads/fastp/{sample}.json')
    conda:
        'envs/mshg.yml'
    resources:
        mem_mb=6000
    shell:
        'fastp '\
            '--low_complexity_filter '\
            '-i {input.cuta_trimmed_reads_R1} '\
            '-I {input.cuta_trimmed_reads_R2} '
            '-o {output.fastp_trimmed_reads_R1} '\
            '-O {output.fastp_trimmed_reads_R2} '\
            '-h {output.report} '\
            '-j {output.json}'


# align reads to reference with bwa mem
rule align_reads_to:
    input:
        trimmed_reads_R1 = data_dir('reads/trimmed/{sample}_R1.fastp.trim.fastq'),
        trimmed_reads_R2 = data_dir('reads/trimmed/{sample}_R2.fastp.trim.fastq'),
        reference = 'data/references/{reference}.fa.gz'
    output:
        sam = data_dir('aligned/to_{reference}/{sample}.sam')
    conda:
        'envs/mshg.yml'
    shell:
        'bwa mem '\
            '{input.reference} '\
            '{input.trimmed_reads_R1} '\
            '{input.trimmed_reads_R2} '\
            '> {output.sam}'


# filter reads:
# 1. Flag -F 3332 exludes: read unmapped, not primary alignment, 
#                          supplementary alignment, read is PCR or optical duplicate
# 2. bioawk -Hc sam '$cigar !~ /^[0-9]+[SH]/': removes all the reads that are 
#                                              soft or hard clipped at the start
# 3. Sort and convert to bam (for IGV and for betools) 
# 4. Index the bam (for IGV)

rule filter_aligned_reads:
    input:
        sam = data_dir('aligned/to_{reference}/{sample}.sam')
    output:
        dedup_bam = temp(data_dir('aligned/to_{reference}/{sample}.dedup.bam')),
        filtered_sam = data_dir('aligned/to_{reference}/{sample}.filtered.sam'),
        filtered_bam = data_dir('aligned/to_{reference}/{sample}.filtered.bam'),
    conda:
        'envs/mshg.yml'
    shell:
        'samtools collate -@ 4 -O -u {input.sam} | '\
        'samtools fixmate -@ 4 -m -u - - | '\
        'samtools sort -@ 4 -u - | '\
        'samtools markdup -@ 4 - {output.dedup_bam} &&'\
        'samtools view -h -F 3332 {output.dedup_bam} | '\
        'bioawk -Hc sam \'$cigar !~ /^[0-9]+[SH]/\' > {output.filtered_sam} && '\
        'samtools view -bS {output.filtered_sam} | '\
        'samtools sort - > {output.filtered_bam} && '\
        'samtools index {output.filtered_bam}'


# make bed file with only first in pair (-f 64), excluding the ones mapped to the reverse strand (-F 16)
# -----------------------------------------------------------------------------

rule make_bed:
    input:
        filtered_bam = data_dir('aligned/to_{reference}/{sample}.filtered.bam')
    output:
        red_bam = data_dir('aligned/to_{reference}/{sample}.filtered.red.bam'),
        red_bed = data_dir('aligned/to_{reference}/{sample}.red.bed')
    conda:
        'envs/mshg.yml'
    shell:
        'samtools view '\
            '-b '\
            '-f 64 '\
            '-F 16 '\
            '{input.filtered_bam} '\
            '> {output.red_bam} && '\
        'samtools index '\
            '{output.red_bam} && '\
        'bedtools bamtobed '\
            '-i {output.red_bam} '\
            '> {output.red_bed}'

# make a coverage file

rule coverage:
    input:
        sam = data_dir('aligned/to_{reference}/{sample}.sam')
    output:
        depth = data_dir('aligned/to_{reference}/{sample}.cov')
    conda:
        'envs/mshg.yml'
    shell:
        'samtools view -bS {input.sam} | samtools sort | samtools depth -a - > {output.depth}'
