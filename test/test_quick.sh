#!/bin/bash
i=data/era5/cdf/S20160922_12
o=output/test_quick.nc
mkdir -pv output
python -m atmcirctools.bin.intp_lvl -i ${i} -o ${o} -l{310,320,330} || exit
ncview ${o} 2>/dev/null
