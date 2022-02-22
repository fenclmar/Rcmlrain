# Rcmlrain
Functions for processing commercial microwave link (CML) data written in R language and an example of spatial reconstruction with virtual CML data.

## Instalation
Rcmlrain works with R language 3.0.0. and newer.  The repo CANNOT be directly compiled to R package, instead the file fun_CML.r needs to be sourced and loaded together with zoo package. 

```
source(url('https://raw.githubusercontent.com/fenclmar/Rcmlrain/master/Rcmlrain/fun_CML.r'))
install.packages('zoo')
library(zoo)
```

## CML functions
The file fun_CML.r contains set of functions used for CML processing, i.e. function for dry-wet classification, baseline identification, estimating wet antenna attenuation, and converting raindrop attenuation to rainfall intensity. The file contains also functions for CML adjustment and rainfall spatial reconstuction as suggested by Fencl et al. (2017) and Goldshtein et al. (2009). 

## Example file with data
The example file ex_CML_spatial.R  shows a CML rainfall spatial reconstructoin according to  Goldshtein et al. (2009). Functions used were coded for the workshop Urban Rainfall from Mobile Phone Networks which was held within the 14th International Conference on Urban Drainage in Prague, on 10th Sep. 2017.

## Running the example file

1. Install R and R-Studio or any other editor.
2. Download the files from R_cml_rain folder and place them into one folder
3. Open the ex_cml_spatial.R and run it

## References:
Fencl et al. 2019 Gauge-adjusted rainfall estimates from commercial microwave links, Hyd. Earth. Sys. Sci.
Goldhstein et al. 2009, Rain Rate Estimation Using Measurements From Commercial Telecommunications Links, IEEE Transactions on Signal Processing
