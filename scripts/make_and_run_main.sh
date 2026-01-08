#!/bin/bash
export CUBACORES=0

cd /hpe-gr4/flore/Azimuthal-moments/
echo "dove sono"
pwd


if [ -f ../Azimuthal_moments ]; then rm ../Azimuthal_moments; fi
make main
if [ $? -ne 0 ]
then
    echo "Compilation failed!"
    exit 1
else
    time ./exe/Azimuthal_moments $(cat ./input/Azimuthal_moments.input)
fi
