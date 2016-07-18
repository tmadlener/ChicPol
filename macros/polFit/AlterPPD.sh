#!/bin/sh
# source /afs/ihep.ac.cn/users/z/zhangll/fs/rootset.sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/results

########## INPUTS ##########

for nState in 6;do

  cp ../../interface/rootIncludes.inc               rootIncludes.inc
  cp ../../interface/commonVar.h    commonVar.h
  cp ../../interface/ToyMC_chi.h        ToyMC.h

  #touch AlterPPD.cc
  #make

  JobID=AlterPPD_Fwk_TnP

  DefaultID=chic$[$nState-5]_30June2016_rejCow_cutDiMu10_chic1Binning_c1massL_0.05

  ShiftID=BiasCorrection_1S2S3SAug12_1Sig
  ShiftResults=0

  nSystematics=2

  #ProbDist 1...gauss, 2...linear
  # NOTE: tmadlener, 18.07.2016 until AlterPPD.cc is refactored the linear option assumes that the
  # Systematic uncertainty is already multiplied by sqrt(12) !!

  SystID1Base=FrameworkTest
  SystID1Specify=AverageSyst
  SystID1Title=FrameworkTotal
  SystID1ProbDist=1

  SystID2Base=FrameworkTest
  SystID2Specify=Efficiencies
  SystID2Title=TnP_Model
  SystID2ProbDist=1

  if [ $nState -eq 4 ]
  then
    rapBinMin=1
    rapBinMax=2
    ptBinMin=1
    ptBinMax=12
  fi
  if [ $nState -eq 5 ]
  then
    rapBinMin=1
    rapBinMax=3
    ptBinMin=1
    ptBinMax=5
  fi

  if [ $nState -eq 6 ]; then
    rapBinMin=1
    rapBinMax=1
    ptBinMin=1
    ptBinMax=5
  fi

  #SystID1Base=TheGreatRun_BKGmodel
  #SystID1Specify=BestSyst_28p_SQRT12 #in the cc file, this values are multiplied by sqrt12
  #SystID1Title=BKGmodel
  #SystID1ProbDist=2
  #
  #SystID2Base=TotalSyst
  #SystID2Specify=TotalSquaredSystMay21_ExceptBKGmodel
  #SystID2Title=TotalSyst_ExceptBKGmodel
  #SystID2ProbDist=1

  #SystID2Base=TheGreatRun_Rho
  #SystID2Specify=ConstSyst
  #SystID2Title=RhoFactor
  #SystID2ProbDist=1

  #SystID3Base=TheGreatRun_Param
  #SystID3Specify=BestSyst
  #SystID3Title=Parametrization
  #SystID3ProbDist=1
  #
  #SystID4Base=TheGreatRun_TnP
  #SystID4Specify=BestSyst
  #SystID4Title=TnP_model
  #SystID4ProbDist=1
  #
  #SystID5Base=TheGreatRun_FrameworkII
  #SystID5Specify=BestSyst_Bkg
  #SystID5Title=FrameworkIII
  #SystID5ProbDist=1
  #
  #SystID6Base=TheGreatRun_FrameworkII
  #SystID6Specify=BestSyst_Sig_NoUnpol
  #SystID6Title=FrameworkII
  #SystID6ProbDist=1
  #
  #SystID7Base=TheGreatRun_FrameworkI
  #SystID7Specify=BestSyst
  #SystID7Title=FrameworkI
  #SystID7ProbDist=1

  ########################################

  cd ${homedir}

  touch AlterPPD.cc
  make AlterPPD

  JobIDDir=${storagedir}/${DefaultID}
  JobIDDirAlterPPD=${storagedir}/${DefaultID}_${JobID}

  mkdir ${JobIDDirAlterPPD}

  rap_=${rapBinMin}
  while [ $rap_ -le ${rapBinMax} ]
  do
    pT_=${ptBinMin}
    while [ $pT_ -le ${ptBinMax} ]
    do

      rm ${JobIDDirAlterPPD}/results_chic$[$nState-5]_rap${rap_}_pT${pT_}.root
      cp ${JobIDDir}/results_chic$[$nState-5]_rap${rap_}_pT${pT_}.root ${JobIDDirAlterPPD}/results_chic$[$nState-5]_rap${rap_}_pT${pT_}.root

      ./AlterPPD ${DefaultID}=DefaultID ${ShiftID}=ShiftID ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title ${SystID1ProbDist}SystID1ProbDist ${SystID2ProbDist}SystID2ProbDist ${SystID3ProbDist}SystID3ProbDist ${SystID4ProbDist}SystID4ProbDist ${SystID5ProbDist}SystID5ProbDist ${SystID6ProbDist}SystID6ProbDist ${SystID7ProbDist}SystID7ProbDist ${SystID8ProbDist}SystID8ProbDist ${basedir}=basedir ${storagedir}=storagedir ${pT_}ptBinMin ${pT_}ptBinMax ${rap_}rapBinMin ${rap_}rapBinMax ${nSystematics}nSystematics ${nState}nState ShiftResults=${ShiftResults}


      pT_=$[pT_+1]
    done
    rap_=$[rap_+1]
  done

done
