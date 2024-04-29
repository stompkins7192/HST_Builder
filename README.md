# HST_Builder
A comprehensive tool to query the public HST database via MAST and return science ready mosiacs starting with FLC and FLT data products. Currently works for ACSWFC, WFC3UVIS, and WFC3IR data products. The code will automatically build and populate the necessary directory trees from where the script is called under the 'HST' directory.

R package dependencies located on Github, follow the installation instructions provided.

ProFound version 1.23.0 or newer, found at https://github.com/asgr/ProFound
ProPane version 1.7.2 or newer, found at https://github.com/asgr/ProPane
Celestial version 1.5.6 or newer, found at https://github.com/asgr/Celestial
magicaxis version 2.4.5 or newer, found at https://github.com/asgr/magicaxis
Rfits version 1.10.9 or newer, found at https://github.com/asgr/Rfits
Rwcs version 1.8.4 or newer, found at https://github.com/asgr/Rwcs

R package dependencies which can be installed directly from CRAN

stringr version 1.5.1 or newer
Cairo version 1.6-2 or newer
data.table version 1.15.0 or newer
doParallel version 1.0.17 or newer
pracma version 2.4.4 or newer
Cairo - Version 1.6-2 or newer
statip - Version 0.2.3 or newer


Python dependencies
Base Python 3.0 or newer 

The latest versions of these packages should be compatable with the code and should be updated prior to working. Running in Anaconda enviroments may require further development.

astroquery
numpy
pandas
astropy

Two example calls are provided below. 

 A small field using WFC3UVIS data and FLC data products on the most common WCS in the list.

    python3 protool.py --RA 54.206 --DEC -53.869 --RAD 0.05 --TELESCOPE 'HST' --CAMERA 'WFC3/UVIS' --FILTER 'F606W' --STAGE 'FLC' --PAD 0 --PS 0 --ASK 'True' --NCPU 1 --EXT 2 --PS 0 --CORES_WARP 1 --CORES 1 --REM 'TRUE' --WCS 'N'


 A small field using WFC3IR data products and FLT data products on the most common WCS in the list.

 python3 protool.py --RA 216.7736 --DEC 57.88066 --RAD 0.01 --TELESCOPE 'HST' --CAMERA 'WFC3/IR' --FILTER 'F140W' --STAGE 'FLT' --PAD 0 --PS 0 --ASK 'True' --NCPU 1 --EXT 2  --CORES_WARP 1 --CORES 1 --REM 'TRUE' --WCS 'N'



The same


 Arguments


   RA -- Target RA in degrees. Numeric.

   DEC -- Target DEC in degrees.  Numeric.

   RAD -- Box search radius in degrees. A value of 0.1 will query a region 0.2 degrees across.

   TELESCOPE -- MAST query argument.  Telescope, simple enough. String. Only supports 'HST' currently.

   CAMERA -- MAST query argument.  either "WFC3/UVIS"  ,  "WFC3/IR"   or    "ACS/WFC"  . Strict formatting requirements.   String.  

   FILTER -- MAST query argument. Filter, simple enough. String.

   STAGE -- MAST query argument. Stage of FITS file, eg 'FLC' or 'FLT' .  ACS and WFC3 UVIS data must be FLC and WFC3 IR data must be FLT.

   PAD -- Number of pixels to chop off each border before ProFound processing begins. A value of 50 will chop 50 pixels off the border of each CCD chip input image before processing begins. Must be either zero or positive.

   PS --  User-desired pixel scale in arcseconds. Set to zero if you want Rfits to calculate it. 

   ASK -- Python argument. Ask the user before MAST download? Useful if you are unsure of what you will get or are short on storage space or time.  'True'   or 'False'.   Case sensetive.

   NCPU -- Number of CPU's for doParallel when making the ProFound maps. There will be a map for each "SCI" extension in your input and ProFound is quite light on memory usage here.

   EXT --  Extension to reference for making footprints. Not terribly important and should be 2 (SCI).

   CORES_WARP -- The cores_warp arugment to ProPane. Integer.

   CORES    --  The cores argument to ProPane. Integer.

   REM  -- Remove indivudual skymaps after ProPane runs? Default true.   'TRUE' or 'FALSE'   . Must be in ALL CAPS because R.

   WCS -- Target WCS system. Default of 'N' will find see the most common WCS system in the headers and use that. Code will fail if an invalid system is chosen or no matches are found. May require research on older WCS systems if desired.


   The WCS argument must be a substring/string that will be found in the WCSNAME header keyword. The 'GAIAeDR3' string should appear in GAIA DR3 alligned WCS. Though one can subsest to "DR3" and it will look for the "DR3" substring in the list of WCSNAMES. 


