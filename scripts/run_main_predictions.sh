#!/bin/bash
#export CUBACORES=0

cd /hpe-gr4/flore/Azimuthal-moments/
echo "dove sono"
pwd

time ./exe/Azimuthal_moments_predictions $(cat ./input/Azimuthal_moments_predictions.input)
