
#' @title Impute summary data
#' @usage impute(mafile, LD_folder, output)
#' @param mafile string, the path of summary data in COJO format
#' @param LD_folder string,  path of LD folder
#' @param output string, path to output
#' @param thresh number, eigen cutoff threshold for imputation, default 0.995
#' @return none, the output is in the output file
#' @export
impute = function(mafile, LD_folder, output, thresh=0.995){

    if(file.exists(output)){
        message("don't need to run: ", output, "exists")
        stop()
    }

    message("messages are redirected to ", output, ".log")
    zz = file(paste0(output, ".log"), "wt")
    sink(zz)
    sink(zz, type="message")

    ma = fread(mafile)
    ma[, D:=2*freq * (1-freq) * N]
    ma[, varps:=D*(N * se^2 + b^2)/N]

    vp_median = median(ma$varps)
    N_median = median(ma$N)


    snpinfo = fread(paste0(LD_folder, "/snp.info"), head=TRUE)
    setnames(snpinfo, c("CHR", "SNP", "GD", "BP", "A1", "A2", "freq", "N", "blk"))
    snpfinal = snpinfo[, .(SNP, A1, A2, freq, blk)]

    comSNP = intersect(ma$SNP, snpfinal$SNP)

    ma_run = ma[SNP %in% comSNP]

    idx = match(ma_run$SNP, snpfinal$SNP)
    message(length(idx), " SNPs in common between GWAS summary and LD")

    bA1A1 = (ma_run$A1 == snpfinal[idx]$A1) & (ma_run$A2 == snpfinal[idx]$A2)
    message(sum(bA1A1), " SNPs set from summary data")

    bA1A2 = (ma_run$A1 == snpfinal[idx]$A2) & (ma_run$A2 == snpfinal[idx]$A1)
    message(sum(bA1A2), " SNPs flipped alleles")

    message(sum(bA1A1 | bA1A2), " SNPs are typed SNPs")

    snpfinal[, b:=NA_real_]
    snpfinal[, se:=NA_real_]
    snpfinal[, p:=NA_real_]
    snpfinal[, N:=NA_real_]
    snpfinal[idx[bA1A1], freq:=ma_run[bA1A1]$freq]
    snpfinal[idx[bA1A1], b:=ma_run[bA1A1]$b]
    snpfinal[idx[bA1A1], se:=ma_run[bA1A1]$se]
    snpfinal[idx[bA1A1], p:=ma_run[bA1A1]$p]
    snpfinal[idx[bA1A1], N:=ma_run[bA1A1]$N]

    snpfinal[idx[bA1A2], freq:=(1.0 - ma_run[bA1A2]$freq)]
    snpfinal[idx[bA1A2], b:=(-ma_run[bA1A2]$b)]
    snpfinal[idx[bA1A2], se:=ma_run[bA1A2]$se]
    snpfinal[idx[bA1A2], p:=ma_run[bA1A2]$p]
    snpfinal[idx[bA1A2], N:=ma_run[bA1A2]$N]

    snpfinal[is.na(N), N:=N_median]


    message("Start imputation...")
    all_ma = list()

    idxBlocks = unique(snpfinal$blk)

    ### check the LDs
    curVars = c(0.9995, 0.999, 0.995)
    usedVar = 0
    for(curVar in curVars){
        ldfile = paste0(LD_folder, "/eig_block1", "_var", curVar, ".bin3")
        if(file.exists(ldfile)){
            usedVar = curVar
            break
        }
            
    }
    if(usedVar == 0){
        stop("The LD path is not valid")
    }

    for(idxBlk in idxBlocks){
        message("==========", idxBlk, "=========")
        ldfile = paste0(LD_folder, "/eig_block", idxBlk, "_var", usedVar, ".bin3")

        ma_Block = snpfinal[blk==idxBlk]
        ma_Block[, r2:=1]

        m = nrow(ma_Block)

        idxtt = ma_Block[is.finite(b), which=TRUE]

        if(length(idxtt) != 0){
            ma_exist = ma_Block[idxtt]
            ma_exist[, z:=b/se]

            message(" Imputing... ")
            Zi = impGa(ldfile, ma_exist$z, idxtt, m, thresh)

            ma_block_n = ma_Block[-idxtt]
            ma_block_n[, z2:=Zi]
            ma_block_n[, base1:=sqrt(2*freq*(1-freq)*(N+z2*z2))]
            ma_block_n[, b2:=z2*sqrt(vp_median)/base1]
            ma_block_n[, se2:=sqrt(vp_median)/base1]
            ma_block_n[, p2:=pchisq(z2*z2, df=1, lower.tail=FALSE)]
            ma_block_n[, r2:=0] # indicate imputed

            ma_Block[is.na(p), b:=ma_block_n$b2]
            ma_Block[is.na(p), se:=ma_block_n$se2]
            ma_Block[is.na(p), r2:=ma_block_n$r2]
            ma_Block[is.na(p), p:=ma_block_n$p2]
        }else{
            ma_Block[is.na(p), b:=0]
            ma_Block[is.na(p), se:=1]
            ma_Block[is.na(p), r2:=-1]
            ma_Block[is.na(p), p:=1]
        }


        all_ma[[idxBlk]] = ma_Block
    }

    ma_out = rbindlist(all_ma)
    ma_out[, blk:=NULL]
    fwrite(ma_out, file=output, sep="\t", quote=F, na="NA")

    message("done")
}
