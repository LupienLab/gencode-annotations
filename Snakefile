# ==============================================================================
# Environment
# ==============================================================================
import os.path as path

vGENCODE = 33

# ==============================================================================
# Meta Rules
# ==============================================================================
rule all:
    input:
        expand(
            "gencode.v{ver}.{ext}.bed",
            ver=vGENCODE,
            ext=["genes.all", "promoters.all"]
        )

# ==============================================================================
# Rules
# ==============================================================================
rule download_gencode:
    output:
        "gencode.v{vGENCODE}.annotation.gff3.gz",
    shell:
        "curl -O ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{vGENCODE}/gencode.v{vGENCODE}.annotation.gff3.gz"

rule all_genes:
    input:
        "gencode.v{vGENCODE}.annotation.gff3",
    output:
        "gencode.v{vGENCODE}.genes.all.bed",
    shell:
        # use either tab or ";" as field separators
        # only keep genes (protein-coding + others)
        "awk '{{FS=\"(\\t|;)\"; OFS=\"\\t\"}}{{if (NR > 5 && $3 == \"gene\"){{gsub(/gene_id=/, \"\", $10); gsub(/gene_name=/, \"\", $12); gsub(/\"/, \"\", $11); print $1, $4, $5, $10, \".\", $7, $12}} }}' {input} | sort -k1,1 -V -k2,2n > {output}"

rule promoters:
    input:
        "gencode.v{vGENCODE}.genes.all.bed",
    output:
        "gencode.v{vGENCODE}.promoters.all.bed",
    params:
        dnstream = 500,
        upstream = 1500
    shell:
        "awk '{{FS=OFS=\"\\t\"}}{{if ($6 == \"+\") {{ print $1, $2 - {params.upstream}, $2 + {params.dnstream}, $4, $5, $6, $7 }} else {{ print $1, $3 - {params.dnstream}, $3 + {params.upstream}, $4, $5, $6, $7 }} }}' {input} > {output}"


# ==============================================================================
# Tools
# ==============================================================================
rule gunzip:
    input:
        "{file}.gz",
    output:
        "{file}"
    wildcard_constraints:
        file = "((?!\\.gz).)*"
    shell:
        "gunzip {input}"
