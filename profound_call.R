library(celestial, quietly = TRUE)
library(foreach, quietly = TRUE)
library(magicaxis, quietly = TRUE)
library(stringr, quietly = TRUE)
library('data.table', quietly = TRUE)
library('Rfits', quietly = TRUE)
library('Rwcs', quietly = TRUE)
library(ProFound, quietly = TRUE)
library(ProPane, quietly = TRUE)

#Call ProFound Function, required to allow sourcing of profound_call.R externally.
#Don't change the input arguments foo and mask, its not necessary and will break things.
call_profound = function(mask, img, kv,mz, ref){

 
  
  # msk = dqs
  # 
  # for(i in 1:length(dqs)){
  #   msk[[i]]$imDat[] = 0L
  #   msk[[i]]$imDat[bitwAnd(dqs[[i]]$imDat, 4096) > 0]= 1L
  #   msk[[i]]$imDat = propaneDilate(msk[[i]]$imDat, size = 3L)
  #   
  #   msk[[i]]$imDat[bitwAnd(dqs[[i]]$imDat, 4) > 0] = 1L
  #   msk[[i]]$imDat[bitwAnd(dqs[[i]]$imDat, 16) > 0] = 1L
  #   msk[[i]]$imDat[bitwAnd(dqs[[i]]$imDat, 128) > 0] = 1L
  #   msk[[i]]$imDat[bitwAnd(dqs[[i]]$imDat, 512) > 0] = 1L
  # }
  
  
  mask = matrix(0, nrow = dim(ref)[1], ncol = dim(ref)[2])
  mask[bitwAnd(ref, 4096) > 0]= 1L
  mask = propaneDilate(mask, size = 3L)
  mask[bitwAnd(ref, 4) > 0]= 1L
  mask[bitwAnd(ref, 16) > 0]= 1L
  mask[bitwAnd(ref, 128) > 0]= 1L
  mask[bitwAnd(ref, 512) > 0]= 1L
  
  #these parameters can be changed.
 #mask = profoundDilate(mask, size = 3, expand = seq(4096,8192))
  
 mask[which(mask > 1, arr.ind = TRUE)] = 1
  
  if(any(is.null(mask))){
    cat(paste0("DQ mask contains NULLS, breaking out!"), "\n")
    break
  }
  
  #image, mask, and keyvalues inputs cannot be changed, Otherwise, 
  #every other profoundProFound setting is fair game.
  #cat(paste0("Input Magzero   ^^^^^^^^^^^^^^  ", mz, "\n"))
  #your image is actually a pointer, lets treat it as such

  profoundProFound(image = img, mask = mask, keepim = TRUE, box = c(floor(dim(img)[1]/3),floor(dim(img)[2]/3)),
                   keyvalues = kv, magzero = mz)
  
  
  
}

