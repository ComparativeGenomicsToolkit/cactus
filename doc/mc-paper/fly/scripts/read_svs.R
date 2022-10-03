library(methods)
library(sveval)

args = commandArgs(TRUE)

multisamps = FALSE
if(length(args)>2){
  multisamps = as.logical(args[3])
}

if(multisamps){
  svs = readSVvcf.multisamps(args[1], keep.ins.seq=TRUE, keep.ref.seq=TRUE, check.inv=TRUE, keep.ids=TRUE, min.sv.size=30)
} else {
  svs = readSVvcf(args[1], keep.ins.seq=TRUE, keep.ref.seq=TRUE, check.inv=TRUE, keep.ids=TRUE, min.sv.size=30)
}

saveRDS(svs, args[2])
