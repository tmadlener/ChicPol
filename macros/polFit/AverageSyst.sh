#!/bin/bash
# source /afs/ihep.ac.cn/users/z/zhangll/fs/rootset.sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
# storagedir=${basedir}/Psi/Data
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/Systematics/test

########## INPUTS ##########

# for nState in 4 5;do
for nState in 6; do

  DataID=SetOfCuts11_chic_30June2016_rejCow_pt10Cut_chic1Binning

  cp ../../interface/rootIncludes.inc               rootIncludes.inc
  # cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
  # cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h
  cp ../../macros/DataFiles/${DataID}/commonVar.h   commonVar.h
  cp ../../interface/ToyMC_chi.h                    ToyMC.h

  SystID=FrameworkTest

  nSystematics=2

  # subtract JobID2 from JobID1 and store the result under SystID
  # Can be used to calculate the uncertainties associated to the background contamination when the lambdas for
  # fLSB=0 and fLSB=100 are present.
  # NOTE: also multiplies the result with a factor of sqrt(12), since this is an assumption that is needed by
  # AlterPPD and PlotResults
  subtractGraphs=true
  moveTo=Background # AverageSystematics stores the result in AverageSyst. this is where the subtracted Graph
                    # gets stored

  JobID1=fLSB0
  JobID2=fLSB100
  JobID3=FrameworkIII
  JobID4=
  JobID5=
  JobID6=
  JobID7=
  JobID8=
  JobID9=


  #SystID=TotalSyst
  #
  #nSystematics=7
  #
  #JobID1=BestSyst_Bkg
  #JobID2=BestSyst_FrameworkI
  #JobID3=BestSyst_Param
  #JobID4=BestSyst_Sig_NoUnpol
  #JobID5=BestSyst_TnP
  #JobID6=ConstSyst
  #JobID7=BestSyst_28p_SQRT12
  #JobID8=
  #JobID9=

  if [ $nState -eq 4 ]
  then
    ptBinMin=1
    ptBinMax=12
  fi
  if [ $nState -eq 5 ]
  then
    ptBinMin=1
    ptBinMax=5
  fi

  if [ $nState -eq 6 ]; then
    ptBinMin=1
    ptBinMax=5
  fi
  if [ $nState -eq 7 ]; then
    ptBinMin=1
    ptBinMax=4
  fi

  ########################################

  cd ${homedir}

  # touch AverageSystematics.cc
  make -B AverageSystematics

  mkdir Systematics
  mkdir Systematics/${SystID}

  SystDir=Systematics/${SystID}/AverageSyst

  mkdir ${SystDir}

  ./AverageSystematics JobID1=${JobID1} JobID2=${JobID2} JobID3=${JobID3} JobID4=${JobID4} JobID5=${JobID5} JobID6=${JobID6} JobID7=${JobID7} JobID8=${JobID8} JobID9=${JobID9} ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState ${nSystematics}nSystematics subtractGraphs=${subtractGraphs}

  if [ ${subtractGraphs} = "true" ]; then
    mkdir -p Systematics/${SystID}/${moveTo}
    mv Systematics/${SystID}/AverageSyst/TGraphResults*.root Systematics/${SystID}/${moveTo}/
  fi

  rm AverageSystematics
done
