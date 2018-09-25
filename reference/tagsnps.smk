
rule all:
    input:
        expand("{dir}/{prefix}{chr}{suffix}.{kb}kb.r2{r2}.tags",
            dir=config['dir]',
            prefix=config['prefix'],
            suffix=config['suffix'],
            kb=config['kb'],
            r2=config[r2],
            chr=range(1,23))


rule extractSNPs:
    input:
        "{dir}/{prefix}{chr}{suffix}.bim"
    output:
        "{dir}/{prefix}{chr}{suffix}_chunk{chunk}.txt"
    shell:
        """
        awk '{print $2}' {input} | sed  -n ''$start','$end'p;'$end'q' > {output}
        """

rule computeTags:
    input:
        chunk="{dir}/{prefix}{chr}{suffix}.txt",
        bim="{dir}/{prefix}{chr}{suffix}.bim",
        fam="{dir}/{prefix}{chr}{suffix}.fam",
        bed="{dir}/{prefix}{chr}{suffix}.bed"
    output:
         {dir}/{prefix}{chr}.{kb}.{r2}.tags
    shell:
    """
    plink --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --tag-r2 {wildcards.r2} \
          --tag-kb {wildcards.kb} \
          --show-tags {input.chunk} \
          --list-all \
          --out {wildcards.dir}/{wildcards.prefix}{wildcards.chr}{wildcards.suffix}.{wildcards.kb}kb.r2{wildcards.r2}.tags
    """
