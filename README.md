# R_cml_rain
- R function for processing commercial microwave link data
- example of spatial reconstruction with virtual data

The repo is currently not ready to be directly compiled to R package, i.e. fun_CML.r file needs to be downloaded and sourced. Code of functions include comments which mostly provide basic info what function does and info about arguments and outputs of the function.

## CML functions
functions are in a fun_CML.r file. It is a set of functions used for CML processing including code used in fencl et al. 2017, Gauge-adjusted rainfall estimates from commercial microwave links, Hyd. Earth. Sys. Sci. or Goldhstein et al 2009, Rain Rate Estimation Using Measurements From Commercial Telecommunications Links, IEEE Transactions on Signal Processing

## Example file with data
Example file ex_CML_spatial.R  shows CML rainfall spatial reconstructoin according to  Goldshtein et al. (2009). Function used were coded by Martin Fencl and adjusted for workshop Urban Rainfall from Mobile Phone Networks which was held within the 14th International Conference on Urban Drainage 14th International Conference on Urban Drainage in Prague, on 10th Sep. 2017.

## Running the example file

1. Install R and R-Studio or any other editor.
2. Download the files from R_cml_rain folder and place them into one folder
3. Open the ex_cml_cpatial.R and run it
