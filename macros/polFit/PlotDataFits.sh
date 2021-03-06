#!/bin/sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
# storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
#storagedir=$basedir/Psi/Data
# storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/results
# storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/testResults
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/MCclosure
datadir_Start=${basedir}/macros/DataFiles
# datadir_Start=/afs/hephy.at/user/t/tmadlener/CMSSW_5_3_11/src/ChiPol/macros/DataFiles

useRefittedChic=true

########## INPUTS ##########

#Take Care of Mean pT in ToyMC.h
NSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4; do

  cp ../../interface/rootIncludes.inc               rootIncludes.inc
  # cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
  # cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h
  # cp ../../interface/effsAndCuts_Psi$[nState-3]S.h  effsAndCuts.h
  cp ../../interface/commonVar.h commonVar.h
  cp ../../interface/ToyMC_chi.h ToyMC.h
  cp ../../interface/effsAndCuts_chi.h effsAndCuts.h

  touch polRapPtPlot.cc
  make

  #for JobID in Psi$[nState-3]S_${NSigma}Sigma_11Dec2012; do
  # for JobID in chic$[$nState-5]_30March2016_ML10_defaultSett chic$[$nState-5]_30March2016_ML30_defaultSett; do
  # for JobID in chic$[$nState-5]_11April2016_nonRefit_MCEff; do
  # for JobID in chic$[$nState-5]_11April2016_useRef_${useRefittedChic}_rejCBs; do
  # for JobID in chic$[$nState-5]_11April2016_useRef_${useRefittedChic}; do
  # for JobID in chic$[$nState-5]_30March2016_ML10_sigEff_postFix; do
  for JobID in jpsi_13May2016_mcClosure; do

    DataID=_jpsi_13May2016_MCclosure_rejCow_false_rejSea_false
    # DataID=_chic_30March2016_ML30 # fChi1MassLow = 0.3
    # DataID=_chic_$(echo $JobID | awk -F'_' -v OFS='_' '{print $2,$3}') # get the DataID from the JobID since all info is stored there
    # DataID=_chic_11April2016_nonRefit_useRef_${useRefittedChic}_rejCB
    # DataID=_chic_11April2016_nonRefit_useRef_${useRefittedChic}

    FidCuts=11
    if [ $nState -eq 4 ]
    then
      ptBinMin=1
      ptBinMax=5
    fi
    if [ $nState -eq 5 ]
    then
      ptBinMin=2
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


    MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

    additionalName=MPV${MPValgo}

    ############################
    # TreeID=chic$[nState-5]
    TreeID=Psi1S
    #datadir=${datadir_Start}/SetOfCuts${FidCuts}${DataID}/Psi$[nState-3]S/tmpFiles
    datadir=${datadir_Start}/SetOfCuts${FidCuts}${DataID}/tmpFiles

    frameSig=1
    polScenSig=3

    frameBkg=1
    polScenBkg=3

    nGenerations=10


    rapBinMin=1 #don't change
    if [ $nState -eq 4 ]
    then
      rapBinMax=1 #don't change
    fi
    if [ $nState -eq 5 ]
    then
      rapBinMax=3 #don't change
    fi
    if [ $nState -gt 5 ]
    then
      rapBinMax=1
    fi

    Jobdir=${storagedir}/${JobID}

    mkdir ${basedir}/macros/polFit/FiguresData
    mkdir ${Jobdir}
    mkdir ${Jobdir}/Figures
    mkdir ${Jobdir}/Figures/${TreeID}
    mkdir ${Jobdir}/Figures/${TreeID}/Figures
    mkdir ${basedir}/macros/polFit/FiguresData/${JobID}
    mkdir ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}

    #cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${Jobdir}/polGenRecFitPlot.cc
    #cp ${basedir}/macros/polFit/polRapPtPlot.cc ${Jobdir}/polRapPtPlot.cc
    #cp ${basedir}/macros/polFit/PlotFinalResults.cc ${Jobdir}/PlotFinalResults.cc
    #cp ${basedir}/macros/polFit/Makefile ${Jobdir}/Makefile
    #cp ${basedir}/macros/polFit/polGen.C ${Jobdir}/polGen.C
    #cp ${basedir}/macros/polFit/polRec.C ${Jobdir}/polRec.C
    #cp ${basedir}/macros/polFit/polFit.C ${Jobdir}/polFit.C
    cp ${basedir}/macros/polFit/polPlot.C ${Jobdir}/polPlot.C
    #
    #cp ../../interface/rootIncludes.inc ${Jobdir}/rootIncludes.inc
    #cp ../../interface/commonVar_Psi$[nState-3]S.h ${Jobdir}/commonVar.h
    #cp ../../interface/ToyMC_Psi$[nState-3]S.h ${Jobdir}/ToyMC.h
    #cp ../../interface/effsAndCuts_Psi$[nState-3]S.h ${Jobdir}/effsAndCuts.h

    #cd ${Jobdir}
    #touch polRapPtPlot.cc
    #make
    #cp ${basedir}/macros/polFit/polRapPtPlot polRapPtPlot_${TreeID}

    echo 'copy finished'

    for nSigma in 1;do #3 2 1

      cd ${Jobdir}
      cp ${basedir}/macros/polFit/polRapPtPlot polRapPtPlot_${TreeID}

      ./polRapPtPlot_${TreeID} ${nSigma}nSigma ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}frameSig ${polScenSig}polScen ${MPValgo}MPValgo ${nGenerations}nGenerations ${TreeID}=TreeID realdata ${Jobdir}=dirstruct ${nState}nState ${datadir}=realdatadir

      mv ${Jobdir}/TGraphResults_${TreeID}_temp.root ${Jobdir}/TGraphResults_${TreeID}.root
      cp ${Jobdir}/TGraphResults_${TreeID}.root ${Jobdir}/TGraphResults_${TreeID}_${nSigma}sigma.root

      cd ${Jobdir}/Figures/${TreeID}

      cp ${basedir}/latex/DataResults_vs_RapPt.tex .
      cp ${basedir}/latex/IndividualFitResults.tex ../../.
      mv ${Jobdir}/ToyNumericalResults.tex .

      pdflatex ToyNumericalResults.tex
      mv ToyNumericalResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/DataNumericalResults_${additionalName}.pdf
      rm *.aux
      rm *.log


      pdflatex DataResults_vs_RapPt.tex
      mv DataResults_vs_RapPt.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/DataResults_vs_RapPt_${additionalName}.pdf

      rm *.aux
      rm *.log
      rm DataResults_vs_RapPt.tex

      rap_=${rapBinMin}
      while [ $rap_ -le ${rapBinMax} ]
      do
        pT_=${ptBinMin}
        while [ $pT_ -le ${ptBinMax} ]
        do

          cd ${Jobdir}/Figures/${TreeID}

          filename=../lph_vs_lth_${TreeID}_rap${rap_}_pT${pT_}.pdf
          if test -s "$filename"
          then
            cd ../..
	    pdflatex "\newcommand\TreeBinID{${TreeID}_rap${rap_}_pT${pT_}}\input{IndividualFitResults.tex}"
	    mv IndividualFitResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/IndividualFitResults_rap${rap_}pt${pT_}_${additionalName}.pdf
          fi


          pT_=$[pT_+1]
        done
        rap_=$[rap_+1]
      done


    done


    cd ${Jobdir}


  done

  rm polRapPtPlot_${TreeID}
  rm IndividualFitResults.tex
  rm *.aux
  rm *.log
  cd ${basedir}/macros/polFit
  rm polRapPtPlot
done
