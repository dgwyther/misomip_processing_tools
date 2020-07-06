# misomip_processing_tools
This repository contains tools for processing output from the MISOMIP (ISOMIP+ and MISOMIP) experiments. Model output is converted into the format requested by Asay-Davis et al., 2016.

## Installation and requirements:
### Interactive scripts
The matlab scripts require the ROMS_Matlab toolbox (https://www.myroms.org/wiki/Matlab_Scripts) for the function `scoord.m`. No other requirements.
The python scripts require toolboxes `netCDF4` and `cython`.
#### Use:
1. Set history file, grid file and output file names at start of `do_make_isomip_output.m` script. This can than be run non-interactively in Matlab.
2. Run python2.7 script: `python ROMS_streamfunctions.py -r HISTORY_FILE.nc -i INTERMEDIATE_OUTPUT_FILE.nc -o FINAL_OUTPUT_FILE.nc`. For example, to make the Ocean2 file, this could be run: `python ROMS_streamfunctions.py -r Ocean2/ocean_his_ocean2.nc -i TEST_Ocean2_COM_ROMSUTAS.nc -o Ocean2_COM_ROMSUTAS.nc`
### Automated bash script
Reqs same as above. use only scripts in bash_script directory.
#### Use:
1. Set options in runOUTPUT.sh
2. execute with ./runOUTPUT.sh (might need to set permissions to u+x)

## To do:
- bash script still needs tidying up. Add path to matlab executable.
- remove naughty code duplication (same .m and .py files for both main directory and bash_script)

## Authors:
Matlab script: David Gwyther
Python script: Xylar Asay-Davis @xylar
