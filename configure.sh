#!/bin/bash

if [ -z ${FINAL_FLASH_BASEDIR+x} ]; then
    export FINAL_FLASH_BASEDIR=$PWD
fi

# setup environment variables
source $FINAL_FLASH_BASEDIR/setup.sh

# setup larlite
source $FINAL_FLASH_BASEDIR/larlite/config/setup.sh

# setup laropencv
source $FINAL_FLASH_BASEDIR/LArOpenCV/setup_laropencv.sh

# setup Geo2D
source $FINAL_FLASH_BASEDIR/Geo2D/config/setup.sh

# setup LArCV
source $FINAL_FLASH_BASEDIR/LArCV/configure.sh

# setup larlitecv
source $FINAL_FLASH_BASEDIR/larlitecv/configure.sh

