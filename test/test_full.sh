#!/bin/bash
i=data/era5/cdf/S20160922_12
o=output/test_full.nc
mkdir -pv output
python -m atmcirctools.bin.intp_lvl --direction both -i ${i} -o ${o} || exit
ncview ${o} 2>/dev/null
