##Query data from MAST and download a .sh script
## Execute shell script to download data

from astroquery.mast import Observations
import os
import numpy as np
import pandas as pd
import argparse
from astropy.io import ascii
import glob
import sys
import shutil
import astropy.units as u
import astropy.table
from astropy.coordinates import SkyCoord
import time
from astropy.table import Table
import subprocess
from astropy.table import vstack
parser = argparse.ArgumentParser(description='Download from MAST')

## Set up input args to run from command line
parser.add_argument('--RA', type=float, nargs="?",
                    help='Central right ascension of catalogue (deg)', default=110.8375)
parser.add_argument('--DEC', type=float, nargs="?", help="Central declination of catalogue (deg)", default=-73.4391)
parser.add_argument('--RAD', type=float, nargs="?", help="Search radius (deg)", default=1.0)
parser.add_argument('--CAMERA', type=str, nargs="*", help="Camera to download", default = ["ALL"])
parser.add_argument('--TELESCOPE', type=str, nargs="*", help="Telescope to download")

parser.add_argument('--NCPU', type=int, nargs="*",
                    help='Number of parallel ProFound map-making instances to run. Less intensive than ProPane', default = 1)

parser.add_argument('--FILTER', type=str, nargs="*", help="Filter to download", default = ["ALL"])

parser.add_argument('--STAGE', type=str, nargs="*",
                    help='Product level (e.g., UN/CAL for un/cal files)')

parser.add_argument('--PAD', type=int, nargs="*",
                    help='integer # of pixels to pad border of final mosaic with. Default of 0 but this can be changed or made negative if needed, zero is fine as of the latest ProTools updates.', default = 0)

parser.add_argument('--EXT', type=int, nargs="*",
                    help='Integer. Extension to reference for making footprints.', default = 2)

parser.add_argument('--PS', type=float, nargs="*",
                    help='Desired Pixel Scale In Arcseconds. Set to zero if you want ProTools to find it for you.', default = 0)

parser.add_argument('--CORES_WARP', type=int, nargs="*",
                    help='Number of cores to use for image warping. Use with caution, default is one.', default = 1)



parser.add_argument('--CORES', type=int, nargs="*",
                    help='Number of cores to use for image stacking. Use with caution, default is one.', default = 1)

parser.add_argument('--REM', type=str, nargs="*",
                    help='Remove Individual Maps After ProPane is made? Default  "True"', default = 'True')

parser.add_argument('--ASK', type=str, nargs="*", help='Ask User For Confirmation To Download Files?', default = ["True"])

args = parser.parse_args()

def query(ref_dir, ra, dec, rad, telescope, camera, filter, stage,  ASK, dl_products=False):
    print("Querying observation")
    ## Set up the search cone
    

    
   # coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")

    ## convenience stub name for saving the astropy tables
    ## replace slashes with underscores so OS doesn't think WFC3/IR is a directory path.

    
    stub_name = str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" + str(round(rad, 2)) + \
                "_" + telescope + "_" + camera.replace("/", "_") + "_" + filter + "_" + stage

    ## where we will put the downloaded frames
    dl_dir = os.path.join(ref_dir, telescope,
                          str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                          str(round(rad, 2)),
                          stage) + str("/")

    ## Make sure that we get all HST frames that are part of HAP or HLA products
    telescopes = [telescope]
    if "HST" in telescope:
        telescopes.append("HAP")
        telescopes.append("HLA")
    
    
    s_ra = [ra - rad, ra + rad]
    s_dec = [dec - rad, dec + rad]
    ## Now the meat of the code. Query MAST.
    if camera == filter == "ALL":
        obs_table = Observations.query_criteria(s_ra = s_ra,
                                                s_dec = s_dec,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"])
    elif filter == "ALL" and camera != "ALL":
        obs_table = Observations.query_criteria(s_ra = s_ra,
                                        s_dec = s_dec,
                                        project=telescopes,
                                        dataproduct_type=["IMAGE", "image"],
                                        instrument_name=[camera, camera + "/IMAGE"])
    elif camera == "ALL" and filter != "ALL":
        obs_table = Observations.query_criteria(s_ra = s_ra,
                                                s_dec = s_dec,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"],
                                                filters = filter)
    else:
        obs_table = Observations.query_criteria(s_ra = s_ra,
                                                s_dec = s_dec,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"],
                                                instrument_name=[camera, camera + "/IMAGE"],
                                                filters = filter)

    
    
    if len(obs_table) == 0:
        pass
    else:
        wd = os.getcwd()
        wd = wd + '/'
        
        # Observations.enable_cloud_dataset(provider='AWS')
        #data_products = Observations.get_product_list(obs_table)
        product_list = [Observations.get_product_list(obs) for obs in obs_table]
        data_products = vstack(product_list) ## replace data_products with this
        products_stage2 = Observations.filter_products(data_products,
                                                       extension="fits",
                                                       productSubGroupDescription=stage
                                                       )

        
        if "HAP-SVM" in set(products_stage2["project"]):
            products_stage2 = products_stage2[products_stage2["project"] == "HAP-SVM"]

        if len(products_stage2) == 0:
            pass
        else:
            
            df = products_stage2.to_pandas()
            ind = np.unique(products_stage2['obsID'], return_index = True)[1]
            df = df.iloc[ind]

            products_stage2 = Table.from_pandas(df)
            
            
            
            if ASK == "True":
                ind = np.unique(products_stage2['obsID'], return_index = True)[1]
              #  print(ind)
                TBD = np.round(np.sum(products_stage2['size'][ind]) / 10**9, decimals=2)
                Nfiles = len(products_stage2['size'][ind])
                
                print(f"You are about to download {Nfiles} images composing of {TBD} GB of Data.")
            while True:
                confirmation = input("Do you want to proceed with the download? (y/n): ").lower()
                if confirmation == "y":
                    break
                elif confirmation == "n":
                    print("Download cancelled by user.", file=sys.stderr)
                    exit(-1)

            os.makedirs(dl_dir, exist_ok=True)

            sh_files = glob.glob(dl_dir + "*.sh")
            for file in sh_files:
                os.remove(file)
            
            
     
           # Observations.download_products(products_stage2, download_dir=dl_dir, curl_flag=True)

            sh_files = glob.glob(dl_dir + "*.csv")
            for file in sh_files:
                os.remove(file)

            ascii.write(data_products, dl_dir + stub_name +
                        '_info_temp.csv', overwrite=True, format='csv')
            ascii.write(products_stage2, dl_dir + stub_name +
                        '_info.csv', overwrite=True, format='csv')

            if dl_products:
                for i in range(0,len(ind)):
                    #(print(df.iloc[i]['obsID']))

                  #  products_stage2 = Table.from_pandas(df.iloc[i])
                    Observations.download_products(products_stage2[i],
                                                   download_dir=dl_dir,
                                                   curl_flag=False,
                                                   cache=True)

                    
            return obs_table

def main(args_dict):

    JUMPROPE_MAST_TOKEN = os.getenv('JUMPROPE_MAST_TOKEN')
    JUMPROPE_DOWNLOAD_DIR = os.getenv('JUMPROPE_DOWNLOAD_DIR')

    if None in [JUMPROPE_MAST_TOKEN, JUMPROPE_DOWNLOAD_DIR]:
        print("Please set ENV variables. ")
        exit()

    my_session = Observations.login(token=str(JUMPROPE_MAST_TOKEN))
    ref_dir = str(JUMPROPE_DOWNLOAD_DIR)
    print(str(JUMPROPE_DOWNLOAD_DIR))
    
   # args_grid = [(a,b,c,d,e)for a in args_dict["TELESCOPE"]
    #                       for b in args_dict["CAMERA"]
     #                      for c in args_dict["FILTER"]
      #                     for d in args_dict["STAGE"]
       #                    for e in args_dict["ASK"]]

 #   for i in range(len(args_grid)):
     #   query(
        #    ref_dir, args_dict["RA"], args_dict["DEC"], args_dict["RAD"],
        #    args_grid[i][0], args_grid[i][1], args_grid[i][2], args_grid[i][3], args_grid[i][4],
        # #   dl_products=True
      #  )

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:

        args_dict = {
            "RA" : args.RA,
            "DEC" : args.DEC,
            "RAD" : args.RAD,
            "TELESCOPE" : args.TELESCOPE,
            "CAMERA" : args.CAMERA,
            "FILTER" : args.FILTER,
            "STAGE" : args.STAGE,
            "PAD" : args.PAD,
            "ASK" : args.ASK,
            "NCPU" : args.NCPU,
            "EXT" : args.EXT,
            "PS" : args.PS,
            "CORES_WARP" : args.CORES_WARP,
            "CORES" : args.CORES,
            "REM" : args.REM
        }
        

main(args_dict)



#make the footprint and the map after downloads and maps are done.
wd = os.getcwd()
ref_dir = wd
wd = wd + '/'


#query(ref_dir, ra, dec, rad, telescope, camera, filter, stage,  ASK, dl_products=False)
query(ref_dir, args.RA, args.DEC, args.RAD, args.TELESCOPE[0], args.CAMERA[0], args.FILTER[0], args.STAGE[0], args.ASK[0], dl_products=True)


print("doing something  at  location    " + str(args.RA) + ' ' + str(args.DEC))

cmd2 = ('Rscript ' + wd + 'protools_process.R '+ str(wd)  + ' ' + str(round(args.RA, 2)) + ' ' 
 + str(round(args.DEC, 2)) + ' ' + str(round(args.RAD, 2)) + ' ' + str(args.CAMERA) + ' ' + str(args.FILTER) + ' ' + str(args.STAGE)
 + ' ' + str(args.TELESCOPE) + ' ' + str(args.NCPU) + ' ' + str(args.PAD) + ' ' + str(args.EXT)+ ' ' + str(args.PS) 
       + ' ' + str(args.CORES_WARP) + ' ' + str(args.CORES) + ' ' + str(args.REM))      

cmd2 = cmd2.replace('[', '').replace(']', '')
print(cmd2)
os.system(cmd2)

#cmd2 = ('Rscript ' + wd + 'make_maps.R '+ wd + ' ' + str(2) + ' ' + str(ra) + ' ' 
# + str(dec) + ' ' + str(rad) + ' ' + camera + ' ' + filter + ' ' + stage
# + ' ' + telescope + ' ' + str(pad))

