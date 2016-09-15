#!/bin/bash

## script for setting up the directory in which the polarization framework will run
## At performs some VERY SIMPLE checks if copying and building are necessary but not at all if the files in polFit have actually changed!

basedir=$1
storagedir=$2
JobID=$3
datadir=$4

polDir=${basedir}/macros/polFit
resultDir=${storagedir}/${JobID}

DIRPRESENT=0
FILESPRESENT=0
EXECMADE=0

## check if the directory already exists and if the executable is built
if [ -d ${storagedir}/${JobID} ]; then
  echo "Result directory already built. Now checking if executable is present."
  DIRPRESENT=1

  if [ -x ${storagedir}/${JobID}/polGenRecFitPlot ]; then
    echo "'polGenRecFitPlot' is here and executable."
    EXECMADE=1
  else
    if [ $(ls -l ${storagedir}/${JobID}/*.cc | wc -l) -eq 3 ] && [ $(ls -l ${storagedir}/${JobID}/*.h | wc -l) -eq 5 ]; then
         echo "'polGenRecFitPlot' is not here but the source files are present."
         FILESPRESENT=1
    fi
  fi
fi

## do the setup
if [ ${DIRPRESENT} -ne 1 ]; then
   mkdir -p ${storagedir}/${JobID}
fi

if [ ${FILESPRESENT} -ne 1 ]; then
  cp ${polDir}/polGenRecFitPlot.cc ${resultDir}/polGenRecFitPlot.cc
  cp ${polDir}/polRapPtPlot.cc ${resultDir}/polRapPtPlot.cc
  cp ${polDir}/PlotFinalResults.cc ${resultDir}/PlotFinalResults.cc
  cp ${polDir}/Makefile ${resultDir}/Makefile
  cp ${polDir}/polGen.C ${resultDir}/polGen.C
  cp ${polDir}/polRec.C ${resultDir}/polRec.C
  cp ${polDir}/polFit.C ${resultDir}/polFit.C
  cp ${polDir}/polPlot.C ${resultDir}/polPlot.C

  cp ${basedir}/interface/rootIncludes.inc ${resultDir}/rootIncludes.inc
  cp ${basedir}${datadir}/commonVar.h ${resultDir}/commonVar.h
  cp ${basedir}/interface/ToyMC_chi.h ${resultDir}/ToyMC.h
  cp ${basedir}/interface/effsAndCuts_chi.h ${resultDir}/effsAndCuts.h
  cp ${basedir}/interface/clarg_parsing.h ${resultDir}/clarg_parsing.h
  cp ${basedir}/interface/optionStructs.h ${resultDir}/optionStructs.h
fi

if [ ${EXECMADE} -ne 1 ]; then
  thisdir=$(pwd)
  cd ${resultDir}
  touch polGenRecFitPlot.cc
  make
  if [ $? -ne 0 ]; then
    echo "ERROR WHILE BUILDING EXECUTABLES FOR POLFIT"
    exit 111
  fi

  cd ${thisdir}
fi
