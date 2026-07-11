#plinkGwasLPDcombined

import gwaslab as gl


#gl.options.set_option("data_directory", "/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/")

threads=8

input_file="../eur.assoc_single.chr@.txt"
output_file="genesis"
plot="genesis.png"
title="GP2 EUR downsampled GWAS, case/control and sex matched with LPD phase 1"

ss = gl.Sumstats(
    input_file,
    snpid="snpID",
    build="38",
    chrom="chr",
    pos="pos",
    ea="alt",
    eaf="freq",
    nea="ref",
    p="Score.pval",
    beta="Est",
    se="Est.SE",
    n="n.obs",
    sep="\t",
    na_values=["NA", "."],
    other=["MAC", "Score", "Score.SE","Score.Stat"],
    verbose=True,
)


#strict data clean

ss.fix_chr(remove=True)
ss.fix_pos(remove=True)
ss.fix_allele(remove=True)
ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True)
ss.normalize_allele(threads=threads)
ss.sort_coordinate()


#harmonize = align NEA to the reference allele in the genome and infer the strand of palindromic SNPs and indels

#align NEA with REF in the refernce fasta file 
#ss.check_ref(ref_seq="/home/duartej3/isilon/0.JF/Backups/LPD_DATA/mergeLPD_trinity_GP2underperfSNPs_pseudoGWAS_QCpipeline/processPhase1/Homo_sapiens_assembly38.fasta")


#QC after variant allignment 
#ss.flip_allele_stats() 
# ss.data.drop(labels=["SNPID"],axis=1,inplace=True) #drop ID and construct a fresh one 
# ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True) #construct new ID
# ss.normalize_allele(threads=threads) 




#infer strand of palindromcs and indels 
# New optimized method: bulk lookup table (faster for large datasets)
#remove_snp 78 means remove palindrmic snps wih mAF unabled to inver (7), and remove palindromic SNPs with no iformation in reference vcf file
#remove_indel 8 means remove indels with no information in reference vcf file
#mode pi means infer strand for palindromic SNPs and indels 

# ss.infer_strand(ref_infer="/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/AMR.ALL.split_norm_af.1kg_30x.hg38.vcf.gz",
#                         ref_alt_freq="AF",
#                         remove_snp="78",
#                         remove_indel="8",
#                         maf_threshold=0.40,
#                         daf_tolerance=0.20,
#                         threads=threads,)


#QC after variant allignment 
# ss.flip_allele_stats() 
# ss.data.drop(labels=["SNPID"],axis=1,inplace=True) #drop ID and construct a fresh one 
# ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True) #construct new ID
# ss.normalize_allele(threads=threads) 


#QC before filling columns and output 
#ss.remove_dup(mode="c",keep='last',keep_col="P", keep_ascend=True)
#ss.check_sanity(eaf=(0, 1), p=(0, 1))
#ss.check_data_consistency()

#fill needed dataa and drop unnecessary columns if wanted 
#ss.fill_data(to_fill=["OR","OR_95L","OR_95U"], overwrite=True)
#ss.data.drop(labels=["BETA","SE"],axis=1,inplace=True)

ss.plot_mqq(
    mode='mqq',
    #cut=14,
    #skip=5,
    title=title, 
    sig_level=1e-6,
    anno_sig_level=1e-6,
    anno="GENENAME",
    build="38",
    sig_line=True,
    additional_line=[5e-8],
    additional_line_color=['black'],
    font_family="DejaVu Sans",
    fontsize=11,
    anno_fontsize=12,
    colors=["#000000","#ABABAB"],
    save=plot, save_kwargs={"dpi":300, "facecolor":"white"}
)

#Tidy output 
ss.fill_data(to_fill=["OR","OR_95L","OR_95U"], overwrite=True)
ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=output_file, fmt="gwaslab", cols=["MAC", "Score", "Score.SE","Score.Stat"])
lead_variants = ss.get_lead(
    windowsizekb=500,
    sig_level=5e-8,
    anno=True,
    build="38",
    source="ensembl")
lead_variants.to_csv("lead_vars.txt", index=False, sep='\t')


