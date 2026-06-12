#!/bin/bash
#export CUBACORES=0

cd /hpe-gr4/sanedda/Azimuthal-moments_yy/
echo "dove sono"
pwd

time ./exe/Azimuthal_moments_predictions_yy $(cat ./input/Azimuthal_moments_predictions.input)
