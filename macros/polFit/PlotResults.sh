#!/bin/sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
#storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir
# storagedir=$basedir/Psi/Data
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/results

########## INPUTS ##########
NSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 6;do

  JobID=chic1_30June2016_rejCow_cutDiMu10_chic1Binning_c1massL_0.05
  additionalName=_chic$[nState-5]

  PlotMatt=0
  PlotCompare=0

  PlotAsymm=0
  PlotFinalData=0
  PlotSystematics=1
  PlotLegend=1
  PlotBG0plots=0
  DeltaTildeplots=0
  SBmSigPlots=0
  SteerIndividuals=0
  BGratioFits=0
  BGratioChi2Fits=0
  rapBinComb=0
  ExtendLegendInX=0
  ShiftInX=0
  PlotVsComp=0
  DrawLatexStuff=0
  DrawPreliminary=0

  DefaultID=${JobID}
  #DefaultID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
  # CompareID1=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
  # CompareID2=MCclosure_July27_Ups1S_MCtruthFineEta_0toP1Sigma
  # CompareID3=MCclosure_July27_Ups1S_MCtruthFineEta_P1toP3Sigma
  # CompareID4=MCclosure_Sept9_Ups3S_3DataSig_GEN_pT_Eff_RECOdata
  nComp=0

  LegendEntryDefID=with_RhoFactor
  # LegendEntryCompID1=no_RhoFactor
  # LegendEntryCompID2=1S_0toP1Sigma
  # LegendEntryCompID3=1S_P1toP3Sigma
  # LegendEntryCompID4=3Sig_RECOdata_GENeff


  nSystematics=4

  if [ $nState -eq 4 ]
  then
    ptBinMin=1
    ptBinMax=12
  fi

  if [ $nState -eq 5 ]
  then
    ptBinMin=1
    ptBinMax=6
  fi

  if [ $nState -eq 6 ]
  then
    ptBinMin=1
    ptBinMax=5
  fi

  if [ $nState -eq 7 ]
  then
    ptBinMin=1
    ptBinMax=4
  fi

  ### Background Polarization plots
  #DefaultID=BG0_Mar19_HighCtauSigCheck3p0
  #CompareID1=BG0_Mar19_LSB_LowCtauSigCheck3p0
  #CompareID2=BG0_Mar19_RSB_LowCtauSigCheck3p0
  #PlotBG0plots=1

  ### Background Polarization plots Split LSB
  #DefaultID=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB4
  #CompareID1=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB1
  #CompareID2=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB2
  #CompareID3=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB3
  #PlotBG0plots=1
  #rapBinComb=1


  SystID1Base=FrameworkTest
  SystID1Specify=FrameworkI
  SystID1Title=FrameworkI

  SystID2Base=FrameworkTest
  SystID2Specify=FrameworkII
  SystID2Title=FrameworkII

  SystID3Base=FrameworkTest
  SystID3Specify=FrameworkIII
  SystID3Title=FrameworkIII

  SystID4Base=FrameworkTest
  SystID4Specify=Efficiencies
  SystID4Title=Efficiencies

  ########################################

  cd ${homedir}

  cp ../../interface/ToyMC_chi.h ToyMC.h
  cp ../../interface/effsAndCuts_chi.h effsAndCuts.h
  cp ../../interface/rootIncludes.inc rootIncludes.inc
  cp ../../interface/commonVar.h commonVar.h

  touch PlotFinalResults.cc
  make

  mkdir -p FinalResults/${JobID}/Chic$[nState-5]

  JobIDDir=FinalResults/${JobID}

  mkdir -p ${JobIDDir}


  cp PlotFinalResults PlotFinalResults_chic$[nState-5]
  ./PlotFinalResults_chic$[nState-5] ${DefaultID}=DefaultID ${CompareID1}=CompareID1 ${CompareID2}=CompareID2 ${CompareID3}=CompareID3 ${CompareID4}=CompareID4 ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${basedir}=basedir ${storagedir}=storagedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nSystematics}nSystematics ${nComp}nComp ${nState}nState ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title PlotMatt=${PlotMatt} PlotAsymm=${PlotAsymm} PlotCompare=${PlotCompare} PlotFinalData=${PlotFinalData} PlotSystematics=${PlotSystematics} PlotLegend=${PlotLegend} PlotBG0plots=${PlotBG0plots} DeltaTildeplots=${DeltaTildeplots} SBmSigPlots=${SBmSigPlots} CompareSyst=${CompareSyst} SteerIndividuals=${SteerIndividuals} BGratioFits=${BGratioFits} BGratioChi2Fits=${BGratioChi2Fits} rapBinComb=${rapBinComb} SetCompStyle=${SetCompStyle} ${LegendEntryDefID}=LegendEntryDefID ${LegendEntryCompID1}=LegendEntryCompID1 ${LegendEntryCompID2}=LegendEntryCompID2 ${LegendEntryCompID3}=LegendEntryCompID3 ${LegendEntryCompID4}=LegendEntryCompID4 ExtendLegendInX=${ExtendLegendInX} ShiftInX=${ShiftInX} PlotVsComp=${PlotVsComp} DrawLatexStuff=${DrawLatexStuff} DrawPreliminary=${DrawPreliminary}
  rm PlotFinalResults_chic$[nState-5]
  rm PlotFinalResults

  cd ${homedir}/FinalResults/${JobID}/Chic$[nState-5]
  if [ ${PlotFinalData} -eq 1 ]
  then
    cp ${basedir}/latex/FinalDataResults_chic$[nState-5].tex ./FinalDataResults.tex
    pdflatex FinalDataResults.tex
    mv FinalDataResults.pdf FinalDataResults${additionalName}.pdf
    cd Figures
    pdflatex FinalNumericalResults.tex
    cd ../
  fi

  if [ ${PlotSystematics} -eq 1 ]
  then
    cp ${basedir}/latex/Systematics_chic$[nState-5].tex ./Systematics.tex
    pdflatex Systematics.tex
    mv Systematics.pdf Systematics${additionalName}.pdf
  fi
  rm *.aux
  rm *.log
  rm *.tex

  #rm -r Figures${additionalName}
  mkdir Figures${additionalName}
  mv Figures/* Figures${additionalName}/

done
