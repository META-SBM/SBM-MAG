import glob

#samples = [s.split('/')[-1].replace('_contigs.sam', '') for s in glob.glob('/mnt/disk1/PROJECTS/DAVID_WGS/bwa_mem_2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/*/*_contigs.sam')]
#samples = ["V350166995_57","V350166995_58" , "V350166995_59"  ,"V350166995_60" , "V350166995_61" , "V350166995_62",  "V350166995_63" , "V350166995_64",  "V350166995_65", 
# "V350166995_66" , "V350166995_67" , "V350166995_68" , "V350166995_69" , "V350166995_70" , "V350166995_71" , "V350166995_72" , "V350166995_81" , "V350166995_82" ,
# "V350166995_83" , "V350166995_84" , "V350166995_85" , "V350166995_86" , "V350166995_87" , "V350166995_88" , "V350166995_89" , "V350166995_90" , "V350166995_91" ,
# "V350166995_92" , "V350166995_93" , "V350166995_94" , "V350166995_95" , "V350166995_96" , "V350174473_1"  , "V350174473_10" , "V350174473_11" , "V350174473_12" ,
# "V350174473_13" , "V350174473_14" , "V350174473_15" , "V350174473_16" , "V350174473_17" , "V350174473_18" , "V350174473_19" , "V350174473_2"  , "V350174473_21" ,
# "V350174473_22" , "V350174473_23" , "V350174473_24" , "V350174473_25" , "V350174473_26" , "V350174473_27" , "V350174473_28" , "V350174473_29" , "V350174473_3"  ,
# "V350174473_30" , "V350174473_31" , "V350174473_32" , "V350174473_33" , "V350174473_34" , "V350174473_35" , "V350174473_36" , "V350174473_4"  , "V350174473_5"  ,
# "V350174473_6"  , "V350174473_7"  , "V350174473_8"  , "V350174473_9"  , "V350181201_29" , "V350181201_30" , "V350181201_31" , "V350181201_32" , "V350181201_33" ,
# "V350181201_34" , "V350181201_35" , "V350181201_36" , "V350181201_37" , "V350181201_38" , "V350181201_39" , "V350181201_40" , "V350181201_41" , "V350181201_42" ,
# "V350181201_43" , "V350181201_44" , "V350181201_45" , "V350181262_100", "V350181262_101", "V350181262_102", "V350181262_103", "V350181262_104", "V350181262_105",
# "V350181262_106", "V350181262_107", "V350181262_108", "V350181262_109", "V350181262_110", "V350181262_111", "V350181262_112", "V350181262_113", "V350181262_114",
# "V350181262_115", "V350181262_116", "V350181262_117", "V350181262_118", "V350181262_119", "V350181262_121", "V350181262_122", "V350181262_124",
# "V350181262_126", "V350181262_127", "V350181262_128", "V350181262_97" , "V350181262_98" , "V350181262_99" , "V350181321_100", "V350181321_101", "V350181321_102",
# "V350181321_103", "V350181321_104", "V350181321_105", "V350181321_106", "V350181321_107", "V350181321_108", "V350181321_109", "V350181321_110", "V350181321_111",
# "V350181321_112", "V350181321_113", "V350181321_114", "V350181321_115", "V350181321_116", "V350181321_118", "V350181321_119", "V350181321_120", "V350181321_121",
# "V350181321_122", "V350181321_123", "V350181321_124", "V350181321_125", "V350181321_126", "V350181321_97" , "V350181321_98" , "V350181321_99" , "V350181323_10" ,
# "V350181323_11" , "V350181323_12" , "V350181323_13" , "V350181323_14" , "V350181323_15" , "V350181323_16" , "V350181323_161", "V350181323_162", "V350181323_163",
# "V350181323_164", "V350181323_165", "V350181323_166", "V350181323_167", "V350181323_17" , "V350181323_18" , "V350181323_19" , "V350181323_20" , "V350181323_22" ,
# "V350181323_23" , "V350181323_24" , "V350181323_25" , "V350181323_26" , "V350181323_27" , "V350181323_28" , "V350181323_29" , "V350181323_30" , "V350181323_31" ,
# "V350181323_32" , "V350181323_9" ]
samples = ["V350166995_161", "V350166995_162" ,"V350166995_163", "V350166995_164", "V350166995_165", "V350166995_166", "V350166995_167", "V350166995_168", "V350166995_169", "V350166995_170",
 "V350166995_171", "V350166995_172", "V350166995_173", "V350166995_174", "V350166995_175",
 "V350166995_176", "V350166995_177", "V350166995_178", "V350166995_179", "V350166995_180",
 "V350166995_181", "V350166995_182", "V350166995_183", "V350166995_184", "V350166995_185",
 "V350166995_186", "V350166995_187", "V350166995_188", "V350166995_189", "V350166995_190",
 "V350166995_191", "V350166995_192", "V350174473_37" , "V350174473_38" , "V350174473_39" ,
 "V350174473_40" , "V350174473_41" , "V350174473_42" , "V350174473_44" , "V350174473_45" ,
 "V350174473_46" , "V350174473_47" , "V350174473_48" , "V350174473_49" , "V350174473_50" ,
 "V350174473_51" , "V350174473_52" , "V350174473_53" , "V350174473_54" , "V350174473_55" ,
 "V350174473_56" , "V350174473_57" , "V350174473_58" , "V350174473_59" , "V350174473_60" ,
 "V350174473_61" , "V350174473_62" , "V350174473_63" , "V350174473_64" , "V350174473_65" ,
 "V350174473_66" , "V350174473_67" , "V350174473_68" , "V350174473_69" , "V350174473_70" ,
 "V350174473_71" , "V350174473_72" , "V350181201_10" , "V350181201_11" , "V350181201_12" ,
 "V350181201_13" , "V350181201_14" , "V350181201_15" , "V350181201_16" , "V350181201_17" ,
 "V350181201_18" , "V350181201_19" , "V350181201_20" , "V350181201_21" , "V350181201_22" ,
 "V350181201_23" , "V350181201_24" , "V350181201_25" , "V350181201_26" , "V350181201_27" ,
 "V350181201_28" , "V350181201_9"  , "V350181262_129", "V350181262_130", "V350181262_131",
 "V350181262_132", "V350181262_133", "V350181262_134", "V350181262_135", "V350181262_136",
 "V350181262_137", "V350181262_138", "V350181262_139", "V350181262_140", "V350181262_141",
 "V350181262_142", "V350181262_143", "V350181262_144", "V350181262_145", "V350181262_146",
 "V350181262_147", "V350181262_148", "V350181262_149", "V350181262_150", "V350181262_151",
 "V350181262_152", "V350181262_153", "V350181262_154", "V350181262_156", "V350181262_158",
 "V350181262_159", "V350181262_160", "V350181321_129", "V350181321_130", "V350181321_131",
 "V350181321_132", "V350181321_133", "V350181321_134", "V350181321_135", "V350181321_136",
 "V350181321_137", "V350181321_138", "V350181321_139", "V350181321_140", "V350181321_141",
 "V350181321_142", "V350181321_143", "V350181321_144", "V350181321_145", "V350181321_147",
 "V350181321_148", "V350181321_149", "V350181321_152", "V350181321_153", "V350181321_154",
 "V350181321_155", "V350181321_156", "V350181321_157", "V350181321_158", "V350181321_159",
 "V350181321_160", "V350181323_33" , "V350181323_34" , "V350181323_35" , "V350181323_36" ,
 "V350181323_37" , "V350181323_38" , "V350181323_39" , "V350181323_40" , "V350181323_41" ,
 "V350181323_42" , "V350181323_43" , "V350181323_44" , "V350181323_45" , "V350181323_46" ,
 "V350181323_47" , "V350181323_48" , "V350181323_49" , "V350181323_50" , "V350181323_51" ,
 "V350181323_53" , "V350181323_54" , "V350181323_55" , "V350181323_57" , "V350181323_58" ,
 "V350181323_59" , "V350181323_60" , "V350181323_61" , "V350181323_62" , "V350181323_63" ,
 "V350181323_64"] 
# 

rule all:
    input:
#       expand("DAVID_WGS/binning/metabat2/megahit/raw__cutadapt_kyhMgi__Trim_L25T25S4_15M50__biobloom_no_match__megahit__bwa_mem2/V350174473_UDB_68_10/V350174473_UDB_68.binning.done")
        expand('/mnt/SSD7TB/HOT_CACHE_MIGHT_BE_DELETED_ANYTIME/DAVID_WGS/taxa/inStrain/{sample}',sample = samples)
        #expand('/mnt/disk1/PROJECTS/DAVID_WGS/binning/semibin2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/drep',sample = samples)
        #expand('/mnt/SSD7TB/HOT_CACHE_MIGHT_BE_DELETED_ANYTIME/DAVID_WGS/binning/semibin2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/CheckM/checkm.txt',sample = samples)
        #expand('/mnt/SSD7TB/HOT_CACHE_MIGHT_BE_DELETED_ANYTIME/DAVID_WGS/assembly/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/contigs.fasta',sample = samples)
# ==================================== MEGAHIT =====================================
#rule megahit:
#    input:
#        r1 = '{fs_prefix}/reads/{preproc}/{sample}_R1.fastq.gz',
#        r2 = '{fs_prefix}/reads/{preproc}/{sample}_R2.fastq.gz'
#    output:
#        out_dir = directory('{fs_prefix}/assembly/megahit/{preproc}/{sample}'),
#        fasta = "{fs_prefix}/assembly/megahit/{preproc}/{sample}/final.contigs.fa"
#    threads:16
#    shell:
#            """
#            megahit -1 {input.r1} -2 {input.r2} --force --out-prefix {output.out_dir} -t {threads}
#  	    """
#
#rule quast:
#    input: 
#        contig = '{fs_prefix}/co_assemble/megahit/{preproc}/{sample}/final.contigs.fa'
#    output:
#        results_assemble = '{fs_prefix}/co_assemble/megahit/{preproc}/{sample}/report.html',
#        out_dir = '{fs_prefix}/co_assemble/megahit/{preproc}/{sample}'
#    shell:
#        "~/quast-5.2.0/metaquast.py {input.contig} -o {output.out_dir} "
#
# ====================================  SPADES =====================================

#rule spades:
#    input: 
#        r1 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R1.fastq.gz',
#        r2 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R2.fastq.gz'
#    output:
#        out = '{fs_prefix}/assembly/spades/{preproc}/{sample}/contigs.fasta'
#    params:
#        out_dir = '{fs_prefix}/assembly/spades/{preproc}/{sample}',
#        correct = '{fs_prefix}/assembly/spades/{preproc}/{sample}/corrected'
#    threads: 10
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/spades.yaml'
#    shell: """/home/kuzmichenko_pa/miniconda3/envs/spades/bin/spades.py --meta --threads {threads} -1 {input.r1} -2 {input.r2} -o {params.out_dir};
#            rm -r {params.correct}"""

# ==================================== ANVIO =====================================
#rule anvio_simplify_names:
#    input: contigs = '{fs_prefix}/assembly/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/contigs.fasta'
#    output: 
#        db         = '{fs_prefix}/anvio/spades/{sample}/anvio_len{length}/{sample}_contigs.db',
#        formatted  = '{fs_prefix}/anvio/spades/{sample}/anvio_len{length}/{sample}_contigs.fasta'
#    threads: 10
#    wildcard_constraints:
#        length= 1500 
#    conda: 'anvio'
#    params:
#        prefix = '{sample}'
#    shell: """
#            anvi-script-reformat-fasta {input.contigs}  -o {output.formatted} -l {wildcards.length}  --simplify-names --prefix {params.prefix}; 
#            anvi-gen-contigs-database -T {threads} -f {output.formatted} -o {output.db};   anvi-run-hmms -T {threads} -c {output.db}; anvi-run-ncbi-cogs -T {threads} -c {output.db};
#            anvi-run-kegg-kofams -c {output.db} -T {threads}; 
#            anvi-run-scg-taxonomy -c {output.db} -T {threads};
#        """
# ==================================== BWA-MEM2 =====================================
#rule bwa_mem2_index:
#    input: 
#        contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta'
#    output:
#        amb = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/contigs.amb'
#    params:
#        index_folder = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/contigs'
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/bwa.yaml'
#    shell: """
#    bwa-mem2 index -p {params.index_folder} {input.contigs};
#    """
#
#
#rule bwa_mem2:
#    input:
#        r1 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R1.fastq.gz',
#        r2 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R2.fastq.gz',
#        contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,
#        amb = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/contigs.amb'
#    output:
#        bam = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bam',
#        bai = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bai'
#    params:
#        index_folder = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/contigs',
#        #sam = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.sam'
#    threads:10
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/bwa.yaml'
#    shell: """
#    bwa-mem2 mem -t {threads} {params.index_folder} {input.r1}  {input.r2} | samtools sort -@ {threads} -o {output.bam}; 
#    samtools index {output.bam} -@ {threads} -o {output.bai};
#    """
#rule samtools_sort_index:
#     input:
#         sam = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.sam'
#     output:
#         bai = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bai',
#         stat = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.idxstats'
#     params: 
#         bam = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bam',
#     shell: """
#     samtools sort {input.sam} -@ 10 -o {params.bam}; samtools index {params.bam} -@ 10 -o {output.bai};
#     samtools idxstats {params.bam} > {output.stat}"""

## ==================================== ANVIO PROFILE =====================================
#rule anvio_profile:
#    input: 
#        db  = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.db',
#        bam = '{fs_prefix}/bwa_mem_2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/{sample}_contigs.bam',
#    output: 
#        profile = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/profile_db/PROFILE.db',
#    params:
#        out_pref = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/profile_db/'
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/anvio.yaml'
#    threads : 10
#    shell: """ 
#            anvi-profile -i {input.bam} -c {input.db}  --output-dir {params.out_pref} --sample-name {wildcards.sample} -W -T {threads}; 
#        """
##rule anvio_merge:
##    input:
##        profile = '{fs_prefix}/bwa_mem_2/megahit/{preproc}/{sample}/anvio_len{length}/{df_sample}/PROFILE.db',
##    output:
##
##    conda:'/mnt/disk1/PROJECTS/DAVID_WGS/envs/anvio.yaml'
##    shell:
##        "anvi-merge */PROFILE.db -o SAMPLES-MERGED -c contigs.db;"
#
#rule anvio_description:
#    input:
#        profile  = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/profile_db/PROFILE.db',
#        contig   = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.db',
#        bin_desc = '{fs_prefix}/binning/semibin2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/{sample}_bins.txt'
#    output:
#        tax  = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_taxonomy.txt',
#        comp = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_completeness.txt'
#    params:
#        collection = 'Arhangelsk'
#    shell:""" conda activate anvio-8;
#        anvi-import-collection  {input.bin_desc} --contigs-mode -p {input.profile} -c {input.contig} -C {params.collection};
#        anvi-estimate-scg-taxonomy -p {input.profile} -c {input.contig} -C {params.collection} -o {output.tax};
#        anvi-estimate-genome-completeness -p {input.profile} -c {input.contig} -C {params.collection} -o {output.comp}
#    """ 
## ==================================== METABAT2 =====================================
#
#rule jgi_summarize_bam_contig_depths:
#    input:
#        contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,
#        bam  = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bam',
#    output:
#        depth = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}.contigs.depth.txt'
#
#    shell: " jgi_summarize_bam_contig_depths --outputDepth {output.depth}  {input.bam}  "

#
#
#rule metabat:
#    input:
#        contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,
#        depth = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}.contigs.depth.txt'
#    output:
#        done = '{fs_prefix}/binning/metabat2/spades/{preproc}/{sample}/{sample}.binning.done'
#    params:
#        out_pref = '{fs_prefix}/binning/metabat2/spades/{preproc}/{sample}/{sample}'
#    threads: 10
#    shell: "metabat2 --numThreads {threads} --minCV 2 --minClsSize 200000 --minContig 1500 -i {input.contigs} -a {input.depth} -o {params.out_pref} ; touch {output.done};"
#
##rule bat:
##    input:
##        contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,
##        folder  = '{fs_prefix}/binning/metabat2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}'
##    output:
##        cat = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/CAT/out.CAT.',       
##        bat = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/BAT/out.BAT.',
##        txt = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/BAT/out.BAT.bin2classification.txt'
##    threads: 10
##    params:
##        db_d = '/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_CAT_database',
##        db_t = '/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_taxonomy',
##        prefix = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/BAT/out.BAT.',
##        faa = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/CAT/out.CAT.concatenated.predicted_proteins.faa',
##        alg = '{fs_prefix}/binning/metabat2/spades/{sample}/anvio_len1500/cat_bat/{sample}/CAT/out.CAT.alignment.diamond',
##    shell:"""
##        CAT contigs -c {input.contigs} -d {params.db_d} -t {params.db_t} -n {threads}  --force --out_prefix {output.cat};
##        CAT bins --bin_folder {input.folder} --out_prefix {output.bat} -d {params.db_d} -t {params.db_t} -n {threads} -s .fa  -p {params.faa} -a {params.alg} --force;
##    """
#
#rule SemiBin2:
#    input:
#        bam     = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.bam',
#        contigs  = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,     
#    output:
#        out_pref = directory('{fs_prefix}/binning/semibin2/spades/{preproc}/{sample}')
#    threads: 10
#    wildcard_constraints:
#        length= 1500 
#    conda:'/mnt/disk1/PROJECTS/DAVID_WGS/envs/SemiBin.yaml'
#    shell: """
#    SemiBin2 \
#    single_easy_bin\
#    --input-fasta {input.contigs} \
#    --input-bam {input.bam} \
#    --output {output.out_pref}  \
#    --environment human_gut"""
#
#rule dastool:
#    input:
#        table1 = '{fs_prefix}/binning/DASTool/TABLES/{sample}_metabat2.tsv',
#        table2 = '{fs_prefix}/binning/DASTool/TABLES/{sample}_semibin2.tsv',
#        contigs  = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' 
#    output:
#        bins = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_DASTool_bins'),
#        seql = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}.seqlength',
#        protein = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_proteins.faa',
#        log = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_DASTool.log',
#        contig2bin = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_DASTool_contig2bin.tsv',
#        summary = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_DASTool_summary.tsv',
#        b6 = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_proteins.faa.all.b6',
#        archaea = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_proteins.faa.archaea.scg',
#        bacteria = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_proteins.faa.bacteria.scg',
#        SCG = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_proteins.faa.findSCG.b6'
#    threads: 5
#    params:
#        basename='{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}'
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/das_tool.yaml'
#    shell:
#        """
#        /home/kuzmichenko_pa/DAS_Tool-1.1.7/DAS_Tool -i {input.table1},{input.table2} \
#         -l MetaBAT2,SemiBin2 \
#         -t {threads} --write_bins --search_engine diamond \
#         -c {input.contigs} \
#         -o {params.basename}
#        """

#rule drep:
#    input:
#        bins_best = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/{sample}_DASTool_bins',
#    output:
#        folder = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep'),
#        data = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/data'),
#        dereplicated_genomes = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/dereplicated_genomes'),
#    threads: 20
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/drep.yaml'
#    shell:""" 
#    dRep dereplicate {output.folder} -g {input.bins_best}/*.fa -p {threads} -comp 75 -con 10 --P_ani 0.9 --S_ani 0.98 ;
#    """ 
#    #dRep dereplicate {output.folder} -g /mnt/disk1/PROJECTS/DAVID_WGS/binning/semibin2/spades/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{wildcards.sample}/output_bins/*.fa -p {threads}
    #mkdir -p {output.folder}
    #rm -rf {output.folder}/*
    
#rule gtdbtk:
#    input:
#        mag_dir = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/dereplicated_genomes/'
#    output:
#        classify_dir = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/classify')
#    params:
#        identify_dir = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/identify'),
#        align_dir = directory('{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/align')
#    threads:10
#    shell:"""
#    gtdbtk identify --genome_dir {input.mag_dir} --out_dir {params.identify_dir} --extension fa --cpus {threads}
#    gtdbtk align --identify_dir {params.identify_dir} --out_dir {params.align_dir} --cpus {threads}
#    gtdbtk classify --genome_dir {input.mag_dir} --align_dir {params.align_dir} --out_dir {output.classify_dir} --skip_ani_screen --extension fa --cpus {threads}
#    """
#
#rule bwa_mem2_bins:
#    input:
#        r1 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R1.fastq.gz',
#        r2 = '/mnt/disk1/PROJECTS/DAVID_WGS/reads/{preproc}/{sample}_R2.fastq.gz',
#        #contigs = '{fs_prefix}/anvio/spades/{sample}/anvio_len1500/{sample}_contigs.fasta' ,
#        #amb = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/contigs.amb'
#        combined_fasta = '{fs_prefix}/binning/DASTool/dastool_out/combined_output.fasta' ,
#        amb = '{fs_prefix}/binning/DASTool/dastool_out/combined_output.fasta.amb'
#    output:
#        bam = '{fs_prefix}/bwa_mem_2/spades/mag/{preproc}/{sample}/{sample}_contigs.bam',
#        bai = '{fs_prefix}/bwa_mem_2/spades/mag/{preproc}/{sample}/{sample}_contigs.bai'
#    #params:
#        #index_folder = '{fs_prefix}/bwa_mem_2/spades/mag/{preproc}/{sample}/contigs',
#        #sam = '{fs_prefix}/bwa_mem_2/spades/{preproc}/{sample}/{sample}_contigs.sam'
#    threads:10
#    conda: '/mnt/disk1/PROJECTS/DAVID_WGS/envs/bwa.yaml'
#    shell: """
#    bwa-mem2 mem -t {threads} {input.combined_fasta} {input.r1}  {input.r2} | samtools sort -@ {threads} -o {output.bam}; 
#    samtools index {output.bam} -@ {threads} -o {output.bai};
#    """

#rule parse:
#    input:
#        dereplicated_genomes = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/dereplicated_genomes'
#    output:
#        stb_file = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/representitve_genomes.stb'
#    conda:'/mnt/disk1/PROJECTS/DAVID_WGS/envs/drep.yaml'
#    shell:"""
#    parse_stb.py --reverse -f {input.dereplicated_genomes}/*.fa  -o {output.stb_file}
#    """
#
rule instrain:
    input:
        bam = '{fs_prefix}/bwa_mem_2/spades/mag/raw__cutadapt_kyhMgi__tmtic_kyhMgi__biobloom_no_match/{sample}/{sample}_contigs.bam',
        reference_fasta = '{fs_prefix}/binning/DASTool/dastool_out/combined_output.fasta',
        stb_file = '{fs_prefix}/binning/DASTool/dastool_out/{sample}/drep/representitve_genomes.stb'
    output:
        mag_output = directory('{fs_prefix}/taxa/inStrain/{sample}'),
    threads: 5
    conda:'/mnt/disk1/PROJECTS/DAVID_WGS/envs/instrain.yaml'
    shell:"""
    inStrain profile {input.bam} {input.reference_fasta} -p {threads} -s {input.stb_file}  -o {output.mag_output} --database_mode 
    """

