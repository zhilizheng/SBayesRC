
#' @title Tidy GWAS summary data
#' @usage tidy(ma_path, LD_folder, output)
#' @param ma_file string, summary data path, COJO format
#' @param ld_folder string,  path to LD folder
#' @param output string, output path
#' @param freq_thresh number, max difference in allele frequency
#' @param N_sd_range, number, filter the per SNP sample size in mean +- N_sd_range * sd 
#' @return none, the output is in the output file
#' @export
tidy = function(ma_file, ld_folder, output, freq_thresh=0.2, N_sd_range=3, rate2pq=0.5){

    message("All messages are redirected to ", output, ".log")
    zz = file(paste0(output, ".log"), "wt")
    sink(zz)
    sink(zz, type="message")
    message("Summary data: ", ma_file)
    message("LD path: ", ld_folder)
    message("Output: ", output)
    message("Freq_thresh: ", freq_thresh)
    message("N_range: ", N_sd_range)
    message("rate2pq: ", rate2pq)

    snpinfo = fread(paste0(ld_folder, "/snp.info"))
    message(nrow(snpinfo), " SNPs in LD information")

    ma = fread(ma_file)

    ma_col_names = colnames(ma)
    val_col_names = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
    if(!all(val_col_names %in% ma_col_names)){
        stop("The summary data is not a valid COJO format")
    }

    message(nrow(ma), " SNPs in summary data")

    ma_val = ma[is.finite(N) & is.finite(b) & is.finite(se) & is.finite(freq) & p>=0 & p<=1]
    message(nrow(ma_val), " valid SNPs in summary data")

    com_snp = intersect(snpinfo$ID, ma_val$SNP)
    message(length(com_snp), " SNPs in common with LD information")

    idx1 = match(com_snp, snpinfo$ID)
    snpinfo_val = snpinfo[idx1]

    idx2 = match(com_snp, ma_val$SNP)
    ma_val2 = ma_val[idx2]

    if(all(ma_val2$SNP==snpinfo_val$ID)){
        message("Consistent SNP IDs")
    }else{
        stop("Duplicated SNP name in summary data")
    }

    bA1A1 = (ma_val2$A1 == snpinfo_val$A1) & (ma_val2$A2 == snpinfo_val$A2)
    bA1A2 = (ma_val2$A1 == snpinfo_val$A2) & (ma_val2$A2 == snpinfo_val$A1)

    bAll = bA1A1 | bA1A2

    message(sum(bAll), " SNPs have consistent alleles (A1, A2) between the summary data and LD")

    ma_val2[, FRQ_ref:=snpinfo_val$A1Freq]
    ma_val2[bA1A2, FRQ_ref:=1-FRQ_ref]

    ma_val3 = ma_val2[bAll]
    ma_val3[, frq_diff:=FRQ_ref - freq]
    ma_val4 = ma_val3[abs(frq_diff) <= freq_thresh]
    message(nrow(ma_val4), " SNPs passed the allele frequency checking with threshold ", freq_thresh)

    fwrite(ma_val4[, ..val_col_names], file=paste0(output, ".full"), quote=F, sep="\t", na="NA")

    mean_N = mean(ma_val4$N)
    sd_N =  sd(ma_val4$N)

    message("Mean N = ", mean_N, ", SD = ", sd_N)

    up_N = mean_N + N_sd_range * sd_N
    low_N = mean_N - N_sd_range * sd_N

    ma_val5 = ma_val4[N>=low_N & N<=up_N]
    message(nrow(ma_val5), " SNPs have sample size within mean +- ", N_sd_range, "SD")
    message("After QC, Median sample size: ", median(ma_val5$N))

    ma_val5[, D:=2 * freq * (1 - freq) * N]
    ma_val5[, varps:=D*(N * se^2 + b^2)/N]
    vary = median(ma_val5$varps)
    ma_val5[, indic:=sqrt( 2 * freq * (1 - freq) / vary * (N*se^2 + b^2))]
    message(" Quantile of rate2pq:")
    print(quantile(ma_val5$indic, na.rm=TRUE))
    message("Var_y: ", vary)
    ma_val6 = ma_val5[indic > (1 - rate2pq) & indic < (1 + rate2pq)] 
    message(nrow(ma_val6), " SNPs remained after QC by rate2pq ", rate2pq)

    fwrite(ma_val6[, ..val_col_names], file=output, quote=F, sep="\t", na="NA")
    message("Done")
}

