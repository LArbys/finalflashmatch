#!/bin/bash

homedir=$PWD
echo "Starting from $homedir"

source setup.sh
source configure.sh

echo "LARLITE: ${LARLITE_BASEDIR}"
echo "LARLITE: ${GEO2D_BASEDIR}"
echo "LARCV: ${LARCV_BASEDIR}"
echo "LAROPENCV: ${LAROPENCV_BASEDIR}"
echo "LARLITECV: ${LARLITECV_BASEDIR}"

cd $LARLITE_BASEDIR
make clean

cd $LARLITE_BASEDIR/UserDev/BasicTool
make clean

cd $LARLITE_BASEDIR/UserDev/RecoTool
make clean

cd $LARLITE_BASEDIR/UserDev/SelectionTool/OpT0Finder
make clean

cd $GEO2D_BASEDIR
make clean

cd $LAROPENCV_BASEDIR
make clean

cd $LARCV_BASEDIR
make clean

cd $LARLITECV_BASEDIR
make clean

cd ../

echo "DONE"
