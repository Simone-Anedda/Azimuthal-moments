#!/bin/bash
export CUBACORES=0

cd /hpe-gr4/flore/Azimuthal-moments/
echo "dove sono"
pwd


if [ -f ../Azimuthal_moments_predictions ]; then rm ../Azimuthal_moments_predictions; fi
make predictions
if [ $? -ne 0 ]
then
    echo "Compilation failed!"
    exit 1
else
    time ./exe/Azimuthal_moments_predictions $(cat ./input/Azimuthal_moments_predictions.input)
fi
