# Some helper functions
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

logger.begin <- function(outPrefix, log2file){
    bLoggerTOKEN <<- log2file
    if(bLoggerTOKEN){
        zzLOGGER <<- file(paste0(outPrefix, ".log"), "wt")
        message("The messages are redirected to ", outPrefix, ".log")
        message("  No output here...")
        sink(zzLOGGER)
        sink(zzLOGGER, type="message")
    }
}

logger.end <- function(){
    if(bLoggerTOKEN){
        sink()
        sink(type="message")
        close(zzLOGGER) 
        message("  Output reactivated.")
    }
}

