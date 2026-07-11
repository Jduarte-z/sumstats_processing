#designed in version 4.1.6
import os
import gwaslab as gl


gl.options.set_option("data_directory", "/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/")

threads=4

input_file="../all_saige_lpdCombined_maf0.01ExceptSncaLrrk2Gba1"
output_file="lpd_combined_maf0.01ExceptGba1SncaLrrk2"


ss = gl.Sumstats(input_file, build="38",    
    rsid="MarkerID",
    chrom="CHR",
    pos="POS",
    ea="Allele2",
    eaf="AF_Allele2",
    nea="Allele1",
    n=7344,
    ncase="N_case",
    ncontrol="N_ctrl",
    beta="BETA",
    se="SE",
    p="p.value",
    sep="\t",
    other=["AC_Allele2"],
    na_values=["NA",".","nan",""],
    verbose=True
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
ss.check_ref(ref_seq="/home/duartej3/isilon/0.JF/Backups/LPD_DATA/mergeLPD_trinity_GP2underperfSNPs_pseudoGWAS_QCpipeline/processPhase1/Homo_sapiens_assembly38.fasta")


#QC after variant allignment 
ss.flip_allele_stats() 
ss.data.drop(labels=["SNPID"],axis=1,inplace=True) #drop ID and construct a fresh one 
ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True) #construct new ID
ss.normalize_allele(threads=threads) 




#infer strand of palindromcs and indels 
# New optimized method: bulk lookup table (faster for large datasets)
#remove_snp 78 means remove palindrmic snps wih mAF unabled to inver (7), and remove palindromic SNPs with no iformation in reference vcf file
#remove_indel 8 means remove indels with no information in reference vcf file
#mode pi means infer strand for palindromic SNPs and indels 

ss.infer_strand(ref_infer="/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/AMR.ALL.split_norm_af.1kg_30x.hg38.vcf.gz",
                        ref_alt_freq="AF",
                        remove_snp="78",
                        remove_indel="8",
                        maf_threshold=0.40,
                        daf_tolerance=0.20,
                        threads=threads,)


#QC after variant allignment 
ss.flip_allele_stats() 
ss.data.drop(labels=["SNPID"],axis=1,inplace=True) #drop ID and construct a fresh one 
ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True) #construct new ID
ss.normalize_allele(threads=threads) 


#QC before filling columns and output 
ss.remove_dup(mode="c",keep='last',keep_col="P", keep_ascend=True)
ss.check_sanity(eaf=(0, 1), p=(0, 1))
ss.check_data_consistency()

#fill needed dataa and drop unnecessary columns if wanted 
ss.fill_data(to_fill=["OR","OR_95L","OR_95U"], overwrite=True)
#ss.data.drop(labels=["BETA","SE"],axis=1,inplace=True)


#Tidy output 
ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=output_file, fmt="gwaslab")
