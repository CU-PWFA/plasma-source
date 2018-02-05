#!/bin/bash

clear

echo "Setting VSim directories"

export VSIM_SCRIPT_DIR=/home/robert/VSim-8.1
export VSIM_CONTENTS_DIR=$VSIM_SCRIPT_DIR/Contents

echo "Running VSimComposer.sh script"

source $VSIM_SCRIPT_DIR/VSimComposer.sh

echo "Setting Vorpal directories"

export VORPAL_BIN_DIR=${VSIM_CONTENTS_DIR}/engine/bin
export VORPAL_LIB_DIR=${VSIM_CONTENTS_DIR}/engine/lib
export VORPAL_SHARE_DIR=${VSIM_CONTENTS_DIR}/engine/share

echo "Running Vorpal"

spacing=(100 150 200)

for i in ${spacing[*]}; 
do
	cp AccelGradient.pre AccelGradient$i.pre
	/home/robert/VSim-8.1/Contents/engine/bin/mpiexec -np 4 /home/robert/VSim-8.1/Contents/engine/bin//vorpal -dt 5.7519025622e-15 -sd -d 40 -n 200 -i AccelGradient$i.pre -o AccelGradient$i -iargs WITNESS_DELAY=$i
done
