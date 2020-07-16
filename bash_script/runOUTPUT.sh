#!/bin/bash

# set output choices
hisName='ocean_his_ISOMIP_PLUS_Ocean0_5.00e-2.nc'
grdName='Ocean1/isomip_plus_ocean1.nc'
outName='Ocean0_COM_ROMSUTAS.nc'
logName='Log.txt'
errName='errLog.txt'

# run matlab
echo "starting matlab script"
( /usr/local/MATLAB/R2020a/bin/matlab -nodisplay -nodesktop -r "try; cd('/home/ubuntu/IceOceanVolume/ISOMIP_PLUS'); fun_do_make_isomip_output('${hisName}','${grdName}','TMP.${outName}'); catch; end; quit" ) > ${logName} 2> ${errName}

# run python
echo "starting python script"
(python ROMS_streamfunctions.py -r ${hisName} -i TMP.${outName} -o ${outName} ) > Py_${logName}

cat ${logName} Py_${logName} > ${logName}; rm Py_${logName}

echo "finished!"
