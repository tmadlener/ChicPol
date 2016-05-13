#!/bin/sh

homedir=$PWD
cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/MCclosure
datadir_Start=${basedir}/macros/DataFiles
# datadir_Start=/afs/hephy.at/user/t/tmadlener/CMSSW_5_3_11/src/ChiPol/macros/DataFiles # dir of old data

########## INPUTS ##########

#Batch submission system: 0/1
useBatch=0

#fracL=50 #in percent #MC closure: 25 for data sigmas, 50 for MC sigmas
#nSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4;do # chic1 = 6, chic2 = 7

  StatVarTotBGfraction=0     #apply statistical fluctuations on f_background
  StatVarTotBGmodel=0        #apply statistical fluctuations on Bg model
  StatVarRho=0               #apply statistical fluctuations on rho factor

  rapBinMin=1
  rapBinMax=1
  ptBinMin=5
  ptBinMax=5

  FidCuts=11

  nEff=1050				  #1050 parametrized truth #100001 parametrization with sigmoid
  #Jpsi #1101 MCtruthFineEta, 1080 MCTnPparam      #1030=soft-1060=tight-1070=mixed-111=soft-112=tight
  UseMCeff=false

  nDileptonEff=1
  UseMCDileptoneff=true

  nRhoFactor=1
  #nRhoFactor=325 ## old
  # nRhoFactor=326 ## newest

  useAmapApproach=false
  nAmap=1                    #frame/state/sigma/ID ( ID= 2 digits )
  nDenominatorAmap=1		     #the number here corresponds to the same notation as nEff
  cutDeltaREllDpt=false      # deltaR cut for Jpsi analysis

  nSample=20000
  nFits=30
  nSkipGen=0
  MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

  useRefittedChic=true
  # JobID=chic$[$nState-5]_30March2016_ML10_sigEff_postFix #fChi1MassLow = 0.1, sigmoid eff # NOTE: DO NOT USE "nEff" in this name!
  JobID=jpsi_13May2016_mcClosure
  # JobID=chic$[$nState-5]_30March2016_ML30_defaultSett #fChi1MassLow = 0.3
  # JobID=chic$[$nState-5]_11April2016_useRef_${useRefittedChic}_rejCBs
  # DataID=_chic_30March2016_ML10 # fChi1MassLow = 0.1
  DataID=_jpsi_13May2016_MCclosure
  # DataID=_chic_11April2016_nonRefit_useRef_${useRefittedChic}_rejCBs
  # DataID=_chic_30March2016_ML30 # fChi1MassLow = 0.3

  datadir=${datadir_Start}/SetOfCuts${FidCuts}${DataID}/tmpFiles
  # TreeID=chic$[nState-5]
  TreeID=Psi1S

  ########################################
  #useCentralFracL=0

  cd ${homedir}

  polScenSig=3
  polScenBkg=3
  frameSig=1
  frameBkg=1
  ConstEvents=15000
  UseConstEv=true
  nGenerations=${nFits}
  gen=false
  rec=false
  fit=true
  plot=false
  NewAccCalc=false

  ScenDir=Default_ScenDir
  mkdir -p ${storagedir}/${JobID}

  cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${storagedir}/${JobID}/polGenRecFitPlot.cc
  cp ${basedir}/macros/polFit/polRapPtPlot.cc ${storagedir}/${JobID}/polRapPtPlot.cc
  cp ${basedir}/macros/polFit/PlotFinalResults.cc ${storagedir}/${JobID}/PlotFinalResults.cc
  cp ${basedir}/macros/polFit/Makefile ${storagedir}/${JobID}/Makefile
  cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/polGen.C
  cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/polRec.C
  cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/polFit.C
  cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/polPlot.C

  cp ../../interface/rootIncludes.inc ${storagedir}/${JobID}/rootIncludes.inc
  cp ../../interface/commonVar.h ${storagedir}/${JobID}/commonVar.h
  cp ../../interface/ToyMC_chi.h ${storagedir}/${JobID}/ToyMC.h
  cp ../../interface/effsAndCuts_chi.h ${storagedir}/${JobID}/effsAndCuts.h

  cd ${storagedir}/${JobID}
  cp ${basedir}/macros/polFit/runDataFits.sh .

  touch polGenRecFitPlot.cc
  make

  rap_=${rapBinMin}
  while [ $rap_ -le ${rapBinMax} ]
  do
    pT_=${ptBinMin}
    while [ $pT_ -le ${ptBinMax} ]
    do

      nGen_=${nSkipGen}
      nGen_=$[nGen_+1]
      nMaxGen=$[nGenerations+nSkipGen]
      while [ $nGen_ -le $nMaxGen ]
      do

        if [ $useBatch -eq 0 ]
        then

          resultfilename=resultsMerged_chic$[nState-5]_rap${rap_}_pT${pT_}.root
          nActualGen=$[nGen_-nSkipGen]
          if [ $nSkipGen -ge 0 ]
          then
            if [ $nActualGen -eq 1 ]
            then
              cp results_chic$[nState-5]_rap${rap_}_pT${pT_}.root ${resultfilename}
            fi
          fi

        fi

        cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_}_Gen${nGen_}
        ./polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}=nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} nAmap=${nAmap} nDenominatorAmap=${nDenominatorAmap} useBatch=${useBatch} StatVarTotBGfraction=${StatVarTotBGfraction} StatVarTotBGmodel=${StatVarTotBGmodel} StatVarRho=${StatVarRho} nState=${nState}
        rm polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_}_Gen${nGen_}


        if [ $useBatch -eq 0 ]
        then

          mv results_chic$[nState-5]_rap${rap_}_pT${pT_}.root results_Fit${nGen_}_chic$[nState-5]_rap${rap_}_pT${pT_}.root

          if [ $nGen_ -eq 1 ]
          then
            cp results_Fit${nGen_}_chic$[nState-5]_rap${rap_}_pT${pT_}.root ${resultfilename}
          fi

          if [ $nGen_ -ge 2 ]
          then
            mv ${resultfilename} BUFFER_${resultfilename}
            hadd -f ${resultfilename} BUFFER_${resultfilename} results_Fit${nGen_}_chic$[nState-5]_rap${rap_}_pT${pT_}.root
            rm BUFFER_${resultfilename}
          fi

          #cp ${resultfilename} results_MergedUpToFit${nGen_}_$[nState-3]SUps_rap${rap_}_pT${pT_}.root

        fi

        nGen_=$[nGen_+1]
      done

      if [ $useBatch -eq 0 ]
      then
        mv ${resultfilename} results_chic$[nState-5]_rap${rap_}_pT${pT_}.root
        cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_}
        ./polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}=nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=true ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap ${nState}nState
        rm polGenRecFitPlot_chic$[nState-5]_rap${rap_}_pt${pT_}

      fi

      pT_=$[pT_+1]
    done
    rap_=$[rap_+1]
  done

done

rm *.so
rm *.d

mkdir ../tmp
