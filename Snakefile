from pathlib import Path
import re
import subprocess
import os
import pandas as pd
import numpy as np
import pyfastx
import janitor
from Bio.Seq import Seq
from io import StringIO


def revcomp(dna_seq):
    return dna_seq[::-1].translate(str.maketrans("ATGC", "TACG"))

def find_orfs(dna):
    pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA|$|.$|..$))')
    orfs = pattern.findall(dna) + pattern.findall(revcomp(dna))
    longest = max(orfs, key=len, default="")
    
    if longest[-3:] in ["TAG", "TAA", "TGA"]:
        return longest[:-3]
    
    remove = len(longest) % 3
    if remove == 0:
        return longest
    elif remove != 0:
        return longest[:-remove]
    else:
        raise AssertionError("Unreachable")

def common_suffix(folder: str) -> str:
    samples = sorted(
        [file.name for file in Path(folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", file.name)]
    )

    suffix = ""
    test = samples[0]
    
    for i in range(1, len(test) + 1):
        index = -i
        if any(sample[index] != test[index] for sample in samples):
            break
        suffix += test[index]
        
    return suffix[::-1]

def paired_reads(folder: str) -> list[str]:
    def common_name(str1: str, str2: str) -> str:
        name = ""
        for a, b in zip(str1, str2):
            if a != b:
                break
            name += a
        return name
    
    samples = sorted(
        [x.stem for x in Path(folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", x.name)]
    )
    
    prefixes = []
    for i in range(0, len(samples), 2):
        read1, read2 = samples[i], samples[i + 1]
        common = common_name(read1, read2)
        prefixes.append(common)
        
    return prefixes

def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1]] for x in fastx]))

    df = (
        pd.DataFrame({"name": reads[0], "sequence": reads[1]})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )
    
    return df

def run_blastn(contigs: str, db: str, temp_file: str) -> pd.DataFrame:
    # In order for blastn to find the files
    os.environ["BLASTDB"] = config["BLASTDB_ENVIRON_VARIABLE"]
    df = pd.read_csv(contigs)
    if df.shape[0] == 0:
        return df
    
    matches = []
    
    for contig in df.itertuples():
        with open(temp_file, "w+") as f:
            print(f">{contig.name}", file=f)
            print(contig.sequence, file=f)
        command = [
            "blastn", "-query", temp_file, "-db", db, "-max_target_seqs", "1",
            "-outfmt", "6 stitle sacc pident slen sstart send sstrand"
        ]
        
        match = subprocess.check_output(command, universal_newlines=True).strip()
        matches.append(match)
        
    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.shape[0] == 0:
        return df
    
    df[["match_name", "accession", "percent_identity", "subject_len", "align_start", "align_end", "subject_strand"]] = (
       df.matches.str.split("\t", expand=True).loc[:, 0:6]
    )
    df = ( 
        df
        .assign(subject_strand=lambda x: [y[0] for y in x.subject_strand.str.split("\n")])
        .assign(sequence=lambda x: [y.sequence if y.subject_strand == "plus" else revcomp(y.sequence) for y in x.itertuples()])
        .assign(align_start_temp=lambda x: x.align_start)
        .assign(align_start=lambda x: [y.align_start if y.subject_strand == "plus" else y.align_end for y in x.itertuples()])
        .assign(align_end=lambda x: [y.align_end if y.subject_strand == "plus" else y.align_start_temp for y in x.itertuples()])
        .drop(columns="align_start_temp")
        .sort_values("align_start")
    )
    return df

    
# --- CONFIG AND SETUP
configfile: "config.yaml"

SAMPLES = paired_reads(config["SAMPLES"])
SUFFIX = common_suffix(config["SAMPLES"])

SAMPLES_FOLDER = config["SAMPLES"]
RESULT_FOLDER = f"{config['RESULTS_FOLDER']}/{Path(SAMPLES_FOLDER).parts[-1]}" 

# --- RUN EVERYTHING

rule all:
    input:
        #blast
        expand(f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.pilon.csv", sample=SAMPLES),
        expand(f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.megahit.csv", sample=SAMPLES),
        #idxstats
        expand(f"{RESULT_FOLDER}/{{sample}}/idxstats/{{sample}}.idxstats", sample=SAMPLES),

        
        
# --- RULES
rule fastp:
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/fastp.log"
    input:
        r1=f"{SAMPLES_FOLDER}/{{sample}}1{SUFFIX}",
        r2=f"{SAMPLES_FOLDER}/{{sample}}2{SUFFIX}",
    output:
        r1=temp(f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r1_trimmed.fq"),
        r2=temp(f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r2_trimmed.fq"),
        html_report=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.html",
        json_report=temp(f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.json"),
    threads: 
        config["THREADS"]
    shell:
        """
        fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1} \
        --out2 {output.r2} \
        --report_title {wildcards.sample} \
        --thread {threads} \
        --html {output.html_report} \
        --json {output.json_report} \
        > {log} 2>&1
        """

# Align to virus db
rule bwa_virusdb:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    params:
        index=config["VIRUS_INDEX_REFSEQ"]
    threads: 
        config["THREADS"]
    output:
        mapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_viraldb.bam",
        mapped_sorted=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_viraldb_sorted.bam",
        filtered_mapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_viraldb_filtered.bam",
        r1=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_mapped_to_virusdb_1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_mapped_to_virusdb_2.fastq",
        stats=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_viraldb.stats",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_viraldb.flagstat",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_virusdb.log"
    run:
        shell("bwa mem -t {threads} {params.index} {input.r1} {input.r2} > {output.mapped} 2> {log}")
        # stats
        shell("samtools flagstat {output.mapped} > {output.flagstat}")
        shell("samtools stats {output.mapped} > {output.stats}")

        # filter out mapped reads
        shell("samtools view -O bam -h -f 3 {output.mapped} > {output.filtered_mapped}")

        # Make into fastq
        shell("samtools fastq -1 {output.r1} -2 {output.r2} {output.filtered_mapped}")
        # if the fastq file is empty, make dummy fastq
        if Path(output.r1).read_text() == "":
            with open(output.r1, "w+") as f:
                print("@M00568:503:000000000-DK3D6:1:1101:15676:1346", file=f)
                print("GGACCCGGCCGTCCCCGGGTCCCCGAGGCCAGAGAGGTGCCGCGCCGGCGCATGTTGGAAAAGGCAGAGCTGGGTCTGGAGTCGGTGATGGGGGAAGGCGGTGGAGAGGCGTCCACGTCACTGGCCTCCTCGTCCGTCCGCCACTGGGCCG", file=f)
                print("+", file=f)
                print("AA@1>A>1ADDD00AFCCEEG/F//EGG?CHAE/0E/A/EGECE/?EG/E/EE/GFDG1<1EFEHFCEHGDBG?GHGF0<FF0<//?<111//<<CFGHGGCG?CGGCGG-:?-EB9CEAA0;F?FB.E/9EECGG?F?@-;@-B/:A?E@", file=f)

            with open(output.r2, "w+") as f:
                print("@M00568:503:000000000-DK3D6:1:1101:15676:1346", file=f)
                print("CTGGAGCCGGTTGTGTTTGGAGCCAAGGCCATCCTGGCCCGCACGACGGCCCAGTGCCGGACTGACGAGGAGGCCAGTGACGTGGACGCCTCTCCACCGCCTTCCCCCATCACCGACTCCAGACCCAGCTCTGCCTTTTCCAACATGCGCC", file=f)
                print("+", file=f)
                print("1>>AACF11>1>FAE0BF0C10000B0AA0FFFH11GFFFECCE?EEECG?@EHC1BF??////11E//>EEE/F/?1FGFG??//F/>CG/1FGBHG?C@/?FHFGGGGBHEHC-<<ACC/0..CC.C0/;:FFHBFF0;00/CG/CE??", file=f)

        shell("samtools sort {output.filtered_mapped} > {output.mapped_sorted}")
        shell("samtools index {output.mapped_sorted}")

# INPUT FOR COMING RULES
FASTQ_1 = rules.bwa_virusdb.output.r1
FASTQ_2 = rules.bwa_virusdb.output.r2
        
rule megahit:
    input:
        r1=FASTQ_1,
        r2=FASTQ_2,
    output:
        # make temp
        contigs=f"{RESULT_FOLDER}/{{sample}}/MEGAHIT/final.contigs.fa",
        csv=f"{RESULT_FOLDER}/{{sample}}/MEGAHIT/{{sample}}.contigs.csv",
    params:
        f"{RESULT_FOLDER}/{{sample}}/MEGAHIT",
        min_len=config["CONTIG_LENGTH"],
    threads: 
        config["THREADS"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/megahit.log"
    run:
        # megahit cannot continue if the dir exists
        # NOTE: Added --presets meta-sensitive --no-mercy
        shell("rm -rf {params}")
        shell("megahit -1 {input.r1} -2 {input.r2} --presets meta-sensitive --no-mercy -o {params} --out-prefix {wildcards.sample} 2> {log}")
        # remove everything except contigs file
        #shell("ls -d -1 {params}/* | grep -v .fa | xargs rm -rf")
        
        # if contigs.fa is empty
        from pathlib import Path
        if Path(output.contigs).read_text() == "":
            with open(output.contigs, "w+") as f:
                print(">DUMMY_CONTIG", file=f)
                print("TTAACCTTGG" * 20, file=f)

        # inline wrangle
        df = fastx_file_to_df(output.contigs)
        df = df.assign(sample_id=wildcards.sample)
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.csv, index=False)

# Assembly with spades
rule spades:
    input:
        r1=FASTQ_1,
        r2=FASTQ_2,
    output:
        contigs_fasta=f"{RESULT_FOLDER}/{{sample}}/SPADES/{{sample}}.contigs.fa",
        contigs_fasta_trimmed=f"{RESULT_FOLDER}/{{sample}}/SPADES/{{sample}}.contigs_trimmed.fa",
        contigs_csv=f"{RESULT_FOLDER}/{{sample}}/SPADES/{{sample}}.contigs.csv",
    params:
        out_dir=f"{RESULT_FOLDER}/{{sample}}/SPADES",
        min_len=config["CONTIG_LENGTH"],
        max_ram=100,
    threads: 
        config["THREADS"]
    run:
        shell("rm -rf {params}")
        shell("rnaviralspades.py -t {threads} -m {params.max_ram} -1 {input.r1} -2 {input.r2} -o {params.out_dir}")
        # remove everything except contigs file
        shell("ls -d -1 {params.out_dir}/* | grep -v contigs.fasta | xargs rm -rf")
        #rename contig file
        shell("mv {params.out_dir}/contigs.fasta {output.contigs_fasta}")
        
        # if contigs.fa is empty
        if Path(output.contigs_fasta).read_text() == "":
            with open(output.contigs_fasta, "w+") as f:
                print(">DUMMY_CONTIG", file=f)
                print("TTAACCTTGG" * 20, file=f)
        
        # inline the wrangle rule
        df = fastx_file_to_df(output.contigs_fasta)
        df = df.assign(sample_id=wildcards.sample)
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.contigs_csv, index=False)

        with open(output.contigs_fasta_trimmed, "a+") as f:
            if df.shape[0] == 0:
                print(">DUMMY_CONTIG", file=f)
                print("TTAACCTTGG" * 20, file=f)
            else:
                for x in df.head(20).itertuples():
                    print(f">{x.name}", file=f)
                    print(f"{x.sequence}", file=f)

rule pilon:
    input:
        contigs=rules.spades.output.contigs_fasta_trimmed,
        r1=FASTQ_1,
        r2=FASTQ_2,
    output:
        contigs_bam=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_contigs.bam",
        improved_contigs=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_improved_contigs.fasta",
        csv=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_improved_contigs.csv"
    params:
        index_folder=f"{RESULT_FOLDER}/{{sample}}/PILON/bwa",
        pilon_folder=f"{RESULT_FOLDER}/{{sample}}/PILON",
        min_len=config["CONTIG_LENGTH"],
    threads: 
        config["THREADS"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_contigs.log"
    run:
        # bwa index
        # need to create dir
        shell("rm -rf {params.index_folder}")
        shell("mkdir {params.index_folder}")
        index_folder = f"{params.index_folder}/{wildcards.sample}"
        shell("bwa index -p {index_folder} {input.contigs}")

        # map reads to index of contigs
        shell("bwa mem -t {threads} {index_folder} {input.r1} {input.r2} 2> {log} | samtools view -h -O bam | samtools sort -o {output.contigs_bam}")
        shell("samtools index {output.contigs_bam}")
        # pilon
        shell("pilon -Xmx50G --threads {threads} --genome {input.contigs} --frags {output.contigs_bam} --outdir {params.pilon_folder} --output {wildcards.sample}_improved_contigs")

        # remove the index
        shell("rm -rf {params.index_folder}")

        #inline wrangle
        df = fastx_file_to_df(output.improved_contigs)
        df = df.assign(sample_id=wildcards.sample)
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.csv, index=False)



# --- Classify contigs
rule blastn:
    input:
        pilon_contigs=rules.pilon.output.csv,
        spades_contigs=rules.spades.output.contigs_csv,
        megahit_contigs=rules.megahit.output.csv,
    output:
        blast_pilon=f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.pilon.csv",
        blast_spades=f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.spades.csv",
        blast_megahit=f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.megahit.csv",
    params:
        temp_file=f"{RESULT_FOLDER}/{{sample}}/BLASTN/temp.blastn",
        db=config["BLASTN_DB"],
    run:
        df_spades = (
            run_blastn(contigs=input.spades_contigs, db=params.db, temp_file=params.temp_file)
            .sort_values("read_len", ascending=False)
            .assign(file_name=f"{wildcards.sample}")
            .assign(assembler="spades")
            .assign(longest_orf_nt=lambda x: [find_orfs(y.sequence) for y in x.itertuples()])
            .assign(longest_orf_AA=lambda x: ["".join(Seq(y.longest_orf_nt).translate()) for y in x.itertuples()])
        )
        
        if df_spades.shape[0] == 0:
            df_spades.to_csv(output.blast_spades, index=False)

        
        else:
            (
                df_spades
                .assign(orf_nt_len=lambda x: x.longest_orf_nt.str.len())
                .assign(orf_AA_len=lambda x: x.longest_orf_AA.str.len())
                .to_csv(output.blast_spades, index=False)
            )

        df_pilon = (
            run_blastn(contigs=input.pilon_contigs, db=params.db, temp_file=params.temp_file)
            .sort_values("read_len", ascending=False)
            .assign(file_name=f"{wildcards.sample}")
            .assign(assembler="pilon_spades")
            .assign(longest_orf_nt=lambda x: [find_orfs(y.sequence) for y in x.itertuples()])
            .assign(longest_orf_AA=lambda x: ["".join(Seq(y.longest_orf_nt).translate()) for y in x.itertuples()])
        )
        
        if df_pilon.shape[0] == 0:
            df_pilon.to_csv(output.blast_pilon, index=False)


        else:
            (
                df_pilon
                .assign(orf_nt_len=lambda x: x.longest_orf_nt.str.len())
                .assign(orf_AA_len=lambda x: x.longest_orf_AA.str.len())
                .to_csv(output.blast_pilon, index=False)
            )

        
        
       

        df_megahit = (
            run_blastn(contigs=input.megahit_contigs, db=params.db, temp_file=params.temp_file)
            .sort_values("read_len", ascending=False)
            .assign(file_name=f"{wildcards.sample}")
            .assign(assembler="megahit")
            .assign(longest_orf_nt=lambda x: [find_orfs(y.sequence) for y in x.itertuples()])
            .assign(longest_orf_AA=lambda x: ["".join(Seq(y.longest_orf_nt).translate()) for y in x.itertuples()])
        )
        
        if df_megahit.shape[0] == 0:
            df_megahit.to_csv(output.blast_megahit, index=False)
            
        else:
            (
                df_megahit
                .assign(orf_nt_len=lambda x: x.longest_orf_nt.str.len())
                .assign(orf_AA_len=lambda x: x.longest_orf_AA.str.len())
                .to_csv(output.blast_megahit, index=False)
            )
        
        if Path(params.temp_file).exists():
            os.remove(params.temp_file)

            
            
rule samtools_idxstats:
    input:
        bam=rules.bwa_virusdb.output.mapped_sorted,
        stats=rules.bwa_virusdb.output.stats,
    output:
        idxstats=f"{RESULT_FOLDER}/{{sample}}/idxstats/{{sample}}.idxstats",
        csv=f"{RESULT_FOLDER}/{{sample}}/idxstats/{{sample}}.idxstats.csv"
    run:
        shell("samtools idxstats {input.bam} > {output.idxstats}")

        total_reads = (
            pd.read_csv(
            f"{input.stats}",
            sep="\t",
            skiprows=7,
            header=None,
            names=["x", "y", "z"],
            on_bad_lines="skip"
            )
            .loc[lambda x: x.x == "SN"]
            .loc[lambda x: x.y == "sequences:"]
            ["z"]
            .squeeze()
        )

        idxstats = (
            pd.read_csv(
            f"{output.idxstats}",
            sep="\t",
            header=None,
            names=["genome", "genome_len", "mapped", "unmapped"]
            )
            .sort_values("mapped", ascending=False)
            .assign(file=f"{wildcards.sample}")
            .assign(total_reads=total_reads)
        )

        idxstats.to_csv(f"{output.csv}", index=False)
        