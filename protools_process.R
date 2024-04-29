library(Rcpp, quietly = TRUE)
library(RANN, quietly = TRUE)
library(NISTunits, quietly = TRUE)
library(pracma, quietly = TRUE)

library(celestial, quietly = TRUE)
library(foreach, quietly = TRUE)
library(magicaxis, quietly = TRUE)
library(stringr, quietly = TRUE)
library('data.table', quietly = TRUE)
library('Rfits', quietly = TRUE)
library('Rwcs', quietly = TRUE)
library(ProFound, quietly = TRUE)
library(ProPane, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(Cairo, quietly = TRUE)
library(statip, quietly = TRUE)
start_time =  Sys.time()



inputargs = commandArgs(TRUE)
stub = as.character(inputargs[1]) #Path to your FITS files.
RAtarg = as.numeric(inputargs[2]) #RA targ for file organization stuff
DECtarg = as.numeric(inputargs[3]) #DEC targ for for file organization stuff
rad = as.numeric(inputargs[4])  #search radius for for file organization stuff
camera = as.character(inputargs[5])  #Camera for for file organization stuff
filter = as.character(inputargs[6])  #Filter for for file organization stuff
stage = as.character(inputargs[7])  #Product Stage for for file organization stuff
telescope = as.character(inputargs[8])  #Telescope for file organization stuff
Ncpu = as.numeric(inputargs[9])  #Number of CPUs to ProFound on. 2 by default if not provided.
pad = as.numeric(inputargs[10])  #Number of pixels to chop off all tile borders. MUST BE POSITIVE OR 0.
ext = as.numeric(inputargs[11]) #extension to look at for making footprint. Should be SCI (2) but left as an input in case your FITS files are weird.
ps = as.numeric(inputargs[12]) #Pixel Scale. If set to zero, Rfits will calculate it for us.
cores_warp = as.numeric(inputargs[13]) # cores_warp argument for ProPane.
cores = as.numeric(inputargs[14]) #number of cores for ProPane cores argument. Product of cores_warp and cores must be < # of images.
rem = as.logical(inputargs[15])  #Remove maps at the end? Useful if we only want the ProPane.
targ_wcs = as.character(inputargs[16]) #new feature. Target WCS. Must be a valid value or substring from the "WCSNAME" header keyword.

#From StackExchange...Custom Round Function Because R's base round(..) rounds 0.005 to 0. We can't
#let that happen!
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}



make_maps = function(stub, recursive = TRUE, RAtarg, DECtarg, rad, camera, filter, stage,
                     telescope, Ncpu = 2, pad = 0, targ_wcs, ext = 2){
  
  
  #before re-allocating stub, set wd files we will source(..)

  setwd(stub)
  #re-assign stub to the correct value given "TELESCOPE"
  stub = paste0(stub, telescope, "/")

  

  
  #these variables are all for file naming
  p1 = RAtarg
  p2 = DECtarg
  p3 = rad
  p1 = round2(RAtarg, 2) #RA targ directory name piece
  p2= round2(DECtarg, 2) #Dec targ directory name piece
  p3 = round2(rad, 2) #search radius directory name piece
  
  
  if(grepl("\\.", as.character(p1)) == FALSE){
    p1 = paste0(p1, ".0")
  }
  
  if(grepl("\\.", as.character(p2)) == FALSE){
    p2 = paste0(p2, ".0")
  }
  
  if(grepl("\\.", as.character(p3)) == FALSE){
    p3 = paste0(p3, ".0")
  }
  
  
  stsci_dir = paste0(p1, "_", p2, "_", p3)  #STSCI's API made this directory
  
  #Camera_Filter for organization. 
  top = (unlist(strsplit(camera, "/")))[1]  #grab "WFC3" instrument name from args.
  temp = (unlist(strsplit(camera, "/")))[2]
  p4 = paste0((unlist(strsplit(camera, "/")))[1],
              (unlist(strsplit(camera, "/")))[2], "_", filter, "_", stage)
  #reassign p4 if stringsplit fails.
  if(is.na(temp)){
    p4 = paste0((unlist(strsplit(camera, "/")))[1], "_", filter, "_", stage)
  }
  if(is.na(top)){
    top = paste0(camera)
  }
  
  dataloc = paste0(stub, stsci_dir, "/", stage, "/") 
  #make sure we only use input that we requested.
  drizzle = list.files(dataloc, pattern = "_info.csv", full.names = TRUE)
  drizzle = read.csv(drizzle)
  expected = length(drizzle[,1]) #length of first column
  
  #carefully point to where the files are

  maploc = paste0(stub, stsci_dir, "/profound_products/",p4,"_maps/") #the code will make this directory
  if(dir.exists(maploc) == FALSE){
    cat(paste0("Making Directory For Outputs At ", maploc, "\n"))
    #Location in which output will be directed.
    
    #Program actually begins here
    (system(paste0("mkdir -p ", maploc))) #Makes directory to write data to, will exit if unsuccessful.
  }
  #Rfits_make_list for all "SCI" and "DQ" extensions
  #goes without saying the number of DQ arrays and SCI arrays must be the same
  inlist = list.files(dataloc, pattern = (".fits"), full.names = TRUE, recursive = TRUE)
  
  Nimages = length(inlist)
  
  
  gdr3 = Rfits_key_scan(inlist, keylist = "WCSNAME", extlist = ext)
  good = which(str_detect(gdr3$WCSNAME, targ_wcs) == TRUE)
  
  
  #default behavior.
  #find most common WCS name if "N" (ignore option) is given and use that.
  if(targ_wcs == "N"){
    targ_wcs = mfv(gdr3$WCSNAME)[1] #ensures we choose "fairly" if there is a tie? MFV can't handle multiple
    #modes.
    good = which(str_detect(gdr3$WCSNAME, targ_wcs) == TRUE)
    inlist = inlist[good]
    
    cat(paste0("User elected to use most frequent WCS system in query. Using ", targ_wcs, " WCSNAME!!"), "\n")
  }
  
  if(length(inlist) == 0){
    cat(paste0("No Images Exist With Your Target WCS, Breaking Out!"), "\n")
   
  }
  
  if(length(inlist) != Nimages){
    cat(paste0("Using ", length(inlist), " of ", length(Nimages), " images based on target WCS."), "\n")
    
  }
 
  inlist = inlist[good]
  cat(paste0("Using ", Nimages, " of ", length(inlist), " images based on target WCS."), "\n")
  if(length(inlist) != 0){
  final_inlist = grep(inlist, pattern ="profound", invert = TRUE) #make sure we do not read in existing maps.
  
  if(length(final_inlist) != 0){
    inlist = inlist[final_inlist]
  }
  final_inlist = grep(inlist, pattern =filter, invert = FALSE, ignore.case = TRUE) #Filter by filter.
  if(length(final_inlist) != 0){
    inlist = inlist[final_inlist]
  }
  final_inlist = grep(inlist, pattern =(unlist(strsplit(camera, "/")))[1], invert = FALSE, ignore.case = TRUE) #Filter by camera
  if(length(final_inlist) != 0){
    inlist = inlist[final_inlist]
  }
  if(length(inlist) != expected){
   # cat(paste0(" Expected ", expected, " Input Frames, Got ", length(inlist), " Frames, Flagging and Skipping Query ", camera,"_",filter, "_",stsci_dir, "\n"))
   # stopifnot(length(inlist) == expected)
  }

  if(top == "WFC3"){
    zap=c('LOOKUP','DP[1-2]')
  }
  if(top == "ACS"){  
    #Zap a big list of depreciated keywords.
    zap=c('LOOKUP','DP[1-2]', 'Lookup', 'APERA1', 'APERA2', 'APERLKA1',
          'APERLKA2', 'AXISCORR', 'BOPO','APERA1','APERA2','APERLKA1','APERLKA2','AXISCORR','BOPOFFA1','BOPOFFA2','MAXCHCNT',
          'OCD1_1','OCD1_2','OCD2_1','OCD2_2','OCRPIX1','OCRPIX2','OCRVAL1','OCRVAL2','OCRTYPE1','OCRTYPE2','ONAXIS1',
          'ONAXIS2','OORIENTA','TARGA1',
          'TARGA2','TDDALPHA','TDDBETA','TDD_CTA','TDD_CYA','TDD_CYB','TDD_CXA','WCSDATE')

  }
  hdr = Rfits_read_all(inlist[1], zap=zap)
  #The reference file from which everything will be based. 
  #all files in "dataloc" must come from the same instrument, filter, etc... , for header scan reasons, and should be.

  scis =  which(names(hdr) == "SCI")  #Locate SCI extension positions
  dqs = which(names(hdr) == "DQ")     #Locate DQ extension positions 
  images = c()
  masks = c()

 # cat(paste0("inlist", inlist, "\n"))
 # cat(paste0("DQS  ", dqs, "\n"))
  for(j in scis){
    #Hopefully a foolproof way to get all extensions matching SCI/DQ. These lists must be the same length.
    images = append(images, Rfits_make_list(filelist = inlist, extlist = j, zap = zap))
    class(images) = "Rfits_list"
  }
  
  for(k in dqs){
    masks = append(masks, Rfits_make_list(filelist = inlist, extlist = k, zap = zap))
    class(masks) = "Rfits_list"
  }
  #If Number of Cores Requested is Unnecessary, reduce.
  if(Ncpu > length(images)){
    cat(paste0("User Requested ", Ncpu, " Cores For ", length(images), " images. Reverting to ", length(images), " cores! ", "\n"))
    Ncpu = length(images)
    
    
  }
  

  
  
  registerDoParallel(cl <- makeCluster(Ncpu, outfile = ""))
  foreach(i = 1:length(images), .verbose = FALSE,
          .packages = c("RANN", "NISTunits", "pracma", "Rcpp", "Rfits", "stringr",
                        "magicaxis", "data.table",
                        "Rwcs", "MASS", "ProFound", "celestial")) %dopar% {
             
                          #    cat(paste(images[[i]]$filename, "\n"))
                          if(top == "WFC3"){               
                            foo = Rfits_read_all(images[[i]]$filename, zap=zap)
                            img = Rfits_point(images[[i]]$filename, ext= images[[i]]$ext, zap=zap)
                            img$keyvalues$EQUINOX == foo[[1]]$keyvalues$EQUINOX 
                        #    msk = Rfits_point(masks[[i]]$filename, ext= images[[i]]$ext, zap=c('LOOKUP','DP[1-2]'))
                          }
                          
                          if(top == "ACS"){  
                            #Zap a big list of depreciated keywords.
                            zap=c('LOOKUP','DP[1-2]', 'Lookup', 'APERA1', 'APERA2', 'APERLKA1',
                                  'APERLKA2', 'AXISCORR', 'BOPO','APERA1','APERA2','APERLKA1','APERLKA2','AXISCORR','BOPOFFA1','BOPOFFA2','MAXCHCNT',
                                  'OCD1_1','OCD1_2','OCD2_1','OCD2_2','OCRPIX1','OCRPIX2','OCRVAL1','OCRVAL2','OCRTYPE1','OCRTYPE2','ONAXIS1',
                                  'ONAXIS2','OORIENTA','TARGA1',
                                  'TARGA2','TDDALPHA','TDDBETA','TDD_CTA','TDD_CYA','TDD_CYB','TDD_CXA','WCSDATE')
                            foo = Rfits_read_all(images[[i]]$filename, zap = zap)
                            img = Rfits_point(images[[i]]$filename, ext= images[[i]]$ext, zap=zap)
                            img$keyvalues$EQUINOX == foo[[1]]$keyvalues$EQUINOX
                         #   msk = Rfits_point(masks[[i]]$filename, ext= images[[i]]$ext, zap=zap)
                            
                            
                          }
                          step = which(inlist == i)
                          
                          ccd = images[[i]]$keyvalues$CCDCHIP
                          if(is.null(ccd)){
                            ccd = 1
                          }
                          
                          #Error handling based on the location of the photometry keywords....
                          PHOTPLAM =  images[[i]]$keyvalues$PHOTPLAM
                          PHOTFLAM = images[[i]]$keyvalues$PHOTFLAM
                          if(!is.numeric(PHOTPLAM) | !is.numeric(PHOTFLAM)){
                            
                            PHOTPLAM =  foo[[1]]$keyvalues$PHOTPLAM
                            PHOTFLAM = foo[[1]]$keyvalues$PHOTFLAM
                          }
                          
                          if(is.numeric(foo[[1]]$keyvalues$EXPTIME) == FALSE){
                            cat(paste0("Your File Has No EXPTIME Flag, Something is very wrong."))
                            a = paste0(step)
                          }
                          
                          
                          if(foo[[2]]$keyvalues$BUNIT == "ELECTRONS/S"){
                            mz = -2.5*log10(PHOTFLAM)-5*log10(PHOTPLAM) - 2.408 
                         #   mz = -2.5*log10(PHOTFLAM)-5*log10(PHOTPLAM) - 2.408 + 2.5*log10(foo[[1]]$keyvalues$EXPTIME)
                          }
                          if(foo[[2]]$keyvalues$BUNIT == "ELECTRONS"){
                            mz = -2.5*log10(PHOTFLAM)-5*log10(PHOTPLAM) - 2.408 + 2.5*log10(foo[[1]]$keyvalues$EXPTIME)
                          }
                           
                        #  mz = -2.5*log10(PHOTFLAM)-5*log10(PHOTPLAM) - 2.408 #+ 2.5*log10(foo[[1]]$keyvalues$EXPTIME)
                          
                          
                          a =  foo[[1]]$keyvalues$ROOTNAME #Assumes your global header has ROOTNAME
                          
                          
                          
                          #check for valid ROOTNAME. Must be a character. File will be named (loop step)_profound if missing.
                          if(is.character(a) == FALSE){
                            cat(paste0("Your FITS File Has No ROOTNAME! Output File Will Be ", step, "_profound.fits"), "\n")
                            a = paste0(step)
                          }
                          
                          pf_filename = paste0(p4, "_", a, "_chip_", ccd, "_profound.fits")
                          
                          output = paste0(maploc, pf_filename)
               
                          if(file.exists(output) == FALSE){
                            
                          #  cat(paste0("Making Map ", (output), "\n"))
                            
                 
                            #some testing to remove SIP CTYPE flags

                            img$keyvalues$EQUINOX = foo[[1]]$keyvalues$EQUINOX
                            
                            
                            
                            
                            source("profound_call.R")  #Source profound_call. 
                              
                            
                            if(pad <= 0){
                              pf = call_profound(mask = masks[[i]][,]$imDat, img = img, kv = images[[i]]$keyvalues,
                                                 mz = mz, ref = masks[[i]][,]$imDat)
                            }
                            
                            if(pad > 0){
                              xrun = images[[i]]$keyvalues$NAXIS1 - pad
                              yrun = images[[i]]$keyvalues$NAXIS2 - pad
                              
                              #update NAXIS for what was chopped off.
                              images[[i]]$keyvalues$NAXIS1 = xrun
                              images[[i]]$keyvalues$NAXIS2 = yrun
                              
                              pf = call_profound(mask = masks[[i]][,]$imDat[pad:xrun,pad:yrun], img = images[[i]][,]$imDat[pad:xrun,pad:yrun], kv = images[[i]]$keyvalues,
                                                 mz = mz, ref = masks[[i]][,]$imDat)
                              
                              
                            }
                            pf$value  #store the output the way the rest of the script expects it to be formatted as.
                            #name the profound map file with enough info to make it unique...hopefully.
                            
                            #Sky Subtracted Image
                            
                            Rfits_write_image((pf$image - pf$sky), filename = paste0(maploc, pf_filename), ext = 1,
                                              create_ext = TRUE, create_file = TRUE, compress = FALSE,
                                              keyvalues = img$keyvalues, keynames = names(img$keyvalues), overwrite_file = TRUE)
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 1,
                                            keyname = "EXTNAME", keyvalue = "SKY_SUB")
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 1,
                                            keyname = "MAGZERO", keyvalue = mz)
                            
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 1,
                                            keyname = "EXPTIME", keyvalue = foo[[1]]$keyvalues$EXPTIME)
                            
                            
                            #Mask Image 
                            Rfits_write_image(pf$mask, filename = paste0(maploc, pf_filename),
                                              keyvalues = img$keyvalues, integer = '16', ext = 2, 
                                              create_ext = TRUE, overwrite_file = FALSE, create_file = FALSE)
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 2,
                                            keyname = "EXTNAME", keyvalue = "MASK")
                            #ProFound Sky
                            Rfits_write_image(pf$sky, filename = paste0(maploc, pf_filename),
                                              keyvalues = img$keyvalues, ext = 3, create_ext = TRUE, overwrite_file = FALSE,
                                              compress = FALSE, create_file = FALSE)
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 3,
                                            keyname = "EXTNAME", keyvalue = "SKY")
                            #Inverse Variance Map
                            Rfits_write_image(pf$skyRMS^-2, filename = paste0(maploc, pf_filename),
                                              keyvalues = img$keyvalues, ext = 4, create_ext = TRUE, create_file = FALSE,
                                              compress = FALSE)
                            #
                            
                            Rfits_write_key(filename = paste0(maploc, pf_filename), ext = 4,
                                            keyname = "EXTNAME", keyvalue = "INVAR")
                            
                            
                            
                            
                            cat(paste0("Finished Step ", i, " of ", length(images), "\n")) 
                            
                          }
   
                        }
  
  invisible(stopCluster(cl))
  #end of function  
  }
}




#makes the footprint. Functional with no issues.
make_grid = function(stub, ext = 2, pad = 0, recursive = TRUE, ps = 0,
                     RAtarg, DECtarg, rad, camera, filter, stage,
                     telescope){
  
  #re-assign stub to the correct value given "TELESCOPE"
  stub = paste0(stub, telescope, "/")
  
  
  path = stub
  pad = 0  #force set edge pad to zero. Depreciated argument which originally padded borders to avoid chopping data out.
  #From StackExchange...Custom Round Function Because R's base round(..) rounds 0.005 to 0. We can't
  #let that happen!

  p1 = RAtarg
  p2 = DECtarg
  p3 = rad
  p1 = round2(RAtarg, 2) #RA targ directory name piece
  p2= round2(DECtarg, 2) #Dec targ directory name piece
  p3 = round2(rad, 2) #search radius directory name piece
  
  
  if(grepl("\\.", as.character(p1)) == FALSE){
    p1 = paste0(p1, ".0")
  }
  
  if(grepl("\\.", as.character(p2)) == FALSE){
    p2 = paste0(p2, ".0")
  }
  
  if(grepl("\\.", as.character(p3)) == FALSE){
    p3 = paste0(p3, ".0")
  }
  temp = (unlist(strsplit(camera, "/")))[2]
  p4 = paste0((unlist(strsplit(camera, "/")))[1],
              (unlist(strsplit(camera, "/")))[2], "_", filter, "_", stage)
  #reassign p4 if stringsplit fails.
  if(is.na(temp)){
    p4 = paste0((unlist(strsplit(camera, "/")))[1], "_", filter, "_", stage)
  }
  #make_maps.R should have made these directories.
  stsci_dir = paste0(p1, "_", p2, "_", p3)
  maploc = paste0(stub, stsci_dir, "/profound_products/",p4,"_maps/") 
  
  l1 = list.files(maploc, pattern = "profound", full.names = TRUE, recursive = TRUE)  #path to where your FITS files are located

  final_inlist = grep(l1, pattern =filter, invert = FALSE, ignore.case = TRUE) #Filter by filter.
  
  cat(paste0("looking for files final_inlist  ", maploc), "\n")
  
  if(length(final_inlist) != 0){
    l1 = l1[final_inlist]
  }
  
  sci = Rfits_make_list(l1, extlist = ext) #make list of frames to draw footprint for
  

  
  
  #User provided pixel scale (in arcseconds per pixel)
  ps = ps
  #if user does NOT provide pixel scale, it will be calculated from the files.
  if(ps == 0){
    ps = sapply(sci[1:length(l1)], function(x){pixscale(x$keyvalues)}) #get pixscale from each image
    ps = mean(ps) #calculate mean pixel scale
  }
  #define a list of RA's/DEC's to populate via loop. sapply / lapply won't work with Rwcs::corners.
  ras = c()
  decs = c()
  for(i in 1:length(l1)){  
    
    korners = corners(sci[[i]])   #Run Rwcs::corners on each input frame.
    ras = append(ras, korners[,1])
    decs = append(decs, korners[,2])
  }
  
  x1 = max(ras) #Define Boundaries
  x2 = min(ras)
  y1 = max(decs)
  y2 = min(decs)
  
  test_matrix = matrix(nrow = 4, ncol = 2, data =0)
  test_matrix[1,] = c(x1,y1) #top right
  test_matrix[2,] = c(x1,y2) #bottom right
  test_matrix[3,] = c(x2,y2) #top left
  test_matrix[4,] = c(x2,y1) #bottom left
  
  #Padded matrix for diagnostic visual purposes.
  padded_matrix = test_matrix
  padded_matrix[1:2,1] = test_matrix[1:2,1] + (pad*ps) / 3600 #pad right side
  padded_matrix[3:4,1] = test_matrix[3:4,1] - (pad*ps) / 3600 #pad left side
  padded_matrix[c(1,4),2] = test_matrix[c(1,4),2] + (pad*ps) / 3600 #pad top
  padded_matrix[c(2,3),2] = test_matrix[c(2,3),2] - (pad*ps) / 3600 #pad bottom
  
  hdr = Rfits_read_header(l1[1], ext = 2)
  RA_list = sapply(sci, function(y){y$keyvalues$CRVAL1})
  DEC_list = sapply(sci, function(y){y$keyvalues$CRVAL2})
  
  
  xrun = ceil(abs(x1 - x2)*(3600/ps)*cosd(mean(DEC_list, na.rm = TRUE))) + pad   #define the width of the grid your mosaic will be built onto
  yrun = ceil(abs(y1 - y2)*(3600/ps)) + pad   #define the height of the grid
  
  #Save a diagnostic 
  mosaics = unlist(strsplit(maploc, "maps/"))[1]
  
  #save PNG with a sensible size, preserve aspect ratio.

  
  
  frame_info = Rfits_key_scan(filelist = l1, keylist = c("CRVAL1", "CRVAL2", "NAXIS1", "NAXIS2", "CD1_1",
                                                         "CD2_1","CD1_2","CD2_2", "BITPIX", "EXPTIME"),
                              extlist = ext)
  
  
  
  target_WCS_test = propaneGenWCS(filelist = l1)

  CairoPNG(filename = paste0(mosaics, p1, "_", p2, "_",p3,"_WCS_footprint.png"), width = as.integer(xrun*(2000/xrun)),
           height = as.integer(yrun*(2000/xrun)), units = "px")
  
  for(i in 1:length(l1)){
    if(i == 1){
      Rwcs_overlap( 
        keyvalues_ref = target_WCS_test, keyvalues_test = sci[[i]]$keyvalues, plot = T
      )
      points(target_WCS_test$CRPIX1, target_WCS_test$CRPIX1, pch = "+", col = "black", cex = 8)
      
    }
    else
      
      Rwcs_overlap(
        keyvalues_ref = target_WCS_test, keyvalues_test = sci[[i]]$keyvalues, plot = T, add = T
      )
  }
  
  dev.off()
  
  CairoPNG(filename = paste0(mosaics, filter, "_", p1, "_", p2, "_",p3, "_footprint.png"), width = as.integer(xrun*(2000/xrun)),
           height = as.integer(yrun*(2000/xrun)), units = "px")
  
  #dummy plot with white points on which the lot is built.
  magplot(corners(sci[[1]])[,1], corners(sci[[1]])[,2], col = "white",
          xlim = c(padded_matrix[3,1],padded_matrix[1,1]), 
          ylim = c(padded_matrix[2,2],padded_matrix[1,2]), cex = 2, side = c(1,2,3,4),
          pch = 1, cex.axis = 1, cex.lab = 1) 
  
  for(i in 1:length(l1)){
    
    polygon(corners(sci[[i]])[,1], corners(sci[[i]])[,2], col = rgb(1,0,0,0.15), border = "black",
            lwd = 1)
    
  }
  
  points(target_WCS_test$CRPIX1, target_WCS_test$CRPIX1, pch = "+", col = "black", cex = 8)
  polygon(padded_matrix[,1], padded_matrix[,2], border = "black", col = rgb(0,0,1,0.15))
  
  #est = (prod(xrun,yrun)*8) / 1e6
  
  #new estimate of peak memory usage (64 bit float*Npix)
  est = prod(xrun+pad,yrun+pad)*64 / 1e6
  est = signif(est/8, 4)
  legend("topleft", bty = "n", border = "n",
         legend = c(paste0("Width ", xrun,  " pix"),
                    paste0("Height ", yrun, " pix")),
         #paste0("Est. Peak Memory Usage ", est, " Mb")),
         y.intersp = 0.9, cex = 4)
  
  
  dev.off()
  
  
}




#Make propane stack, has issues.
make_propane = function(stub, ps = 0, pad, ext = 2, RAtarg, DECtarg, rad, camera, filter, stage,
                        telescope, visit = 0, cores_warp = 1, cores = 1, rem = TRUE){
  

  setwd(stub)
  #re-assign stub to the correct value given "TELESCOPE"
  stub = paste0(stub, telescope, "/")
  
  #From StackExchange...Custom Round Function Because R's base round(..) rounds 0.005 to 0. We can't
  #let that happen!

  p1 = RAtarg
  p2 = DECtarg
  p3 = rad
  p1 = round2(RAtarg, 2) #RA targ directory name piece
  p2= round2(DECtarg, 2) #Dec targ directory name piece
  p3 = round2(rad, 2) #search radius directory name piece
  
  
  if(grepl("\\.", as.character(p1)) == FALSE){
    p1 = paste0(p1, ".0")
  }
  
  if(grepl("\\.", as.character(p2)) == FALSE){
    p2 = paste0(p2, ".0")
  }
  
  if(grepl("\\.", as.character(p3)) == FALSE){
    p3 = paste0(p3, ".0")
  }
  
  temp = (unlist(strsplit(camera, "/")))[2]
  p4 = paste0((unlist(strsplit(camera, "/")))[1],
              (unlist(strsplit(camera, "/")))[2], "_", filter, "_", stage)
  #reassign p4 if stringsplit fails.
  if(is.na(temp)){
    p4 = paste0((unlist(strsplit(camera, "/")))[1], "_", filter, "_", stage)
  }
  
  stsci_dir = paste0(p1, "_", p2, "_", p3)
  #the code should make these directories.
  maploc = paste0(stub, stsci_dir, "/profound_products/",p4,"_maps/") 
  
  l1 = list.files(maploc, pattern = "_profound.fits", full.names = TRUE)
  #temporary
 # l1 = l1[40:96]
  
  #read in the necessary info from the profound maps
  sci = Rfits_make_list(l1, extlist = 1) #sky-subtracted image
  masks = Rfits_make_list(l1, extlist = 2)  #dilated mask
  sky = Rfits_make_list(l1, extlist = 3)  #profound-sky from image
  invar = Rfits_make_list(l1, extlist = 4) #profound invar map from image


  
  
  top = (unlist(strsplit(camera, "/")))[1]  #grab "WFC3" instrument name from args.
  if(is.na(top)){
    top = paste0(camera)
  }
  
  filelist = matrix(data = 0, ncol = 2, nrow = length(sci))
  names = c("Image", "Path")
  filelist = data.frame(filelist)
  colnames(filelist) = names
  
  for(i in 1:length(sci)){
  filelist$Image[i] = i
  filelist$Path[i] = sci[[i]]$filename
    
  }
  write.table(filelist, file = paste0(stub, stsci_dir, "/", p4, "_", stsci_dir, "_image_list.csv"),
              row.names = FALSE, sep = ",")
  
  
  #keyscan the WCS info out of input headers, from the skysub image
  
  frame_info = Rfits_key_scan(filelist = l1, keylist = c("BITPIX", "MAGZERO", "EXPTIME"),
                              extlist = 1)
  
 # cat(paste0("Input Magzero   ^^^^^^^^^^^^^^  ", frame_info$MAGZERO, "\n"))
  
  ps = ps
  #if user does NOT provide pixel scale, it will be calculated from the files.
  if(ps == 0){
    ps = sapply(sci[1:length(l1)], function(x){pixscale(x$keyvalues)}) #get pixscale from each image
    ps = max(ps) #choose maximum pixel scale.
  }
  
  
  #identify what will become the target WCS of your mosaic
  
  target_WCS_test = propaneGenWCS(filelist = l1)

  if(prod(cores_warp, cores) > length(sci)){
    cat(paste0("User Requested ", prod(cores_warp, cores), " cores for ", length(sci), " images", "\n"))
    cat(paste0("Reverting to 1 cores and ", length(sci), " cores_warp.", "\n"))
    cores_warp = length(sci)
    cores = 1
    
  }
  
 

  
  drizzle = propaneStackWarpInVar(image_list = sci, 
                                  mask_list = masks,
                                  inVar_list = invar,
                                  magzero_in = frame_info$MAGZERO, #pass in vector of magzero_in
                                  cores = cores, cores_warp = cores_warp,
                                  keyvalues_out = target_WCS_test,
                                  magzero_out = 23.9, exp_list = frame_info$EXPTIME)
  
  
  
  #update NAXIS values 
  drizzle$image$keyvalues$NAXIS1# = xmax - xmin
  drizzle$image$keyvalues$NAXIS2 #= ymax - ymin
  drizzle$image$keyvalues$BITPIX = -32 #write.image whines if it doesn't have BITPIX, so we provide it.
  #  drizzle$weight$imDat$keyvalues$BITPIX = 16
  # drizzle$inVar$imDat$keyvalues$BITPIX = -32

  #write out the propane product as multi-extension FITS, name as RAtarg_DECtarg_propane.fits
  #since there is no other unique identifier.

  
  
  #write the image
  mosaics = unlist(strsplit(maploc, "maps/"))[1]
  out_filename = paste0(mosaics, p1, "_", p2, "_",p3,"_propane.fits")
  Rfits_write_image(drizzle$image$imDat, filename = out_filename,
                    keyvalues = drizzle$image$keyvalues,
                    create_ext = TRUE, create_file = TRUE, ext = 1,
                    keynames = names(drizzle$image$keyvalues),
                    compress = FALSE)
  Rfits_write_key(filename = out_filename, ext = 1,
                  keyname = "EXTNAME", keyvalue = "SCI")
  
  Rfits_write_image(drizzle$weight$imDat, filename = out_filename,
                    keyvalues = drizzle$image$keyvalues, ext = 2, 
                    create_ext = TRUE, overwrite_file = FALSE, create_file = FALSE,
                    integer = '16',
                    compress = FALSE)
  Rfits_write_key(filename = out_filename, ext = 2,
                  keyname = "EXTNAME", keyvalue = "WHT")
  
  Rfits_write_image(drizzle$inVar$imDat, filename = out_filename,
                    keyvalues = drizzle$image$keyvalues, 
                    create_ext = TRUE, overwrite_file = FALSE, create_file = FALSE, ext =3,
                    compress = FALSE)
  Rfits_write_key(filename = out_filename, ext = 3,
                  keyname = "EXTNAME", keyvalue = "INVAR")
  
  
  Rfits_write_image(drizzle$exp$imDat, filename = out_filename,
                    keyvalues = drizzle$image$keyvalues, 
                    create_ext = TRUE, overwrite_file = FALSE, create_file = FALSE, ext =3,
                    compress = FALSE, integer = 16)
  Rfits_write_key(filename = out_filename, ext = 4,
                  keyname = "EXTNAME", keyvalue = "EXPMAP")
  
  
  
  

 # invisible(gc(reset = TRUE))
  #Renaming ProPane not needed without SKYSURF lists.
  # if(visit != 0){
  #   setwd(paste0(stub, stsci_dir, "/profound_products/"))
  #   file2rename = list.files(path = paste0(stub, stsci_dir, "/profound_products"), pattern = paste0( filter, "_", stage, "_", stsci_dir,"_propane.fits"), full.names = FALSE,
  #                            recursive = FALSE)
  #   
  #   instrument = str_replace(camera, "/", "")
  #   cat(paste0("Renaming   ", file2rename, "  to   ", instrument, "_", filter,  "_Visit_", visit, "_propane.fits", "\n"))
  #   
  #   system(paste0("mv ", file2rename, " ", paste0(instrument, "_", filter,  "_Visit_", visit, "_propane.fits")))  
  #   
  # }
  
  if(rem == TRUE){
   l2 = str_replace_all(l1, "//", "/")

    cat(paste0("Removing Individual Maps", "\n"))
    for(i in 1:length(l2)){system(paste0("rm ", l2[i]))}
    
    
    
  }
  

  gc()

}



#Function calls
make_maps(stub, recursive = TRUE, RAtarg, DECtarg, rad, camera, filter, stage, telescope,
          Ncpu, pad, targ_wcs)  

make_grid(stub, ext = 2, pad, recursive = TRUE, ps, RAtarg, DECtarg, rad, camera, filter, stage,
          telescope)


make_propane(stub, ps, pad, ext, RAtarg, DECtarg, rad, camera, filter, stage, telescope, visit, cores_warp, cores, rem)

